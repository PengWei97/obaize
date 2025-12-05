###############################################
# Phase-field shrinking of a circular void
# Only AC equation in the electrode:
#   - Double-well + gradient energy for xi
# No Li transport, no electrochemistry, no mechanics.
#
# Goal:
#   Test whether a void inside the Li electrode evolves
#   under AC curvature-driven dynamics, with a dummy
#   Li occupancy variable theta_m prepared for later use.
#
# Author: Wei Peng
# Date:   2025.12.05
###############################################

# Output naming
my_filename = s2_c1_pf_with_void1
my_interval = 20


###############################################################
#  Mesh
#  - 2D square: 100 μm × 100 μm
#  - Resolution: 250 × 250 QUAD4 elements
#  - Block 1: "electrode"   (0 – 50 μm, Li metal)
#  - Block 2: "electrolyte" (50 – 100 μm, SSE)
###############################################################
[Mesh]
  [./mesh]
    type      = GeneratedMeshGenerator
    dim       = 2
    nx        = 250
    ny        = 250
    xmin      = 0.0
    xmax      = 100e-6     # [m]
    ymin      = 0.0
    ymax      = 100e-6     # [m]
    elem_type = QUAD4
  [../]

  # Subdomain 1: Li electrode, left half (0 – 50 μm)
  [./electrode_block]
    type        = SubdomainBoundingBoxGenerator
    input       = mesh
    block_id    = 1
    block_name  = 'electrode'
    bottom_left = '0.0 0.0 0.0'
    top_right   = '50e-6 100e-6 0.0'
  [../]

  # Subdomain 2: solid electrolyte, right half (50 – 100 μm)
  [./electrolyte_block]
    type        = SubdomainBoundingBoxGenerator
    input       = electrode_block
    block_id    = 2
    block_name  = 'electrolyte'
    bottom_left = '50e-6 0.0 0.0'
    top_right   = '100e-6 100e-6 0.0'
  [../]

  type = MeshGeneratorMesh
[]


###############################################################
#  Variables
#
#  xi       : phase-field order parameter
#             1 = solid (Li), 0 = void
#
#  theta_m  : Li site occupancy (dimensionless), here used
#             only as a dummy variable with a time derivative
#             term (no diffusion, no source). It stays constant
#             in time in this test problem.
###############################################################
[Variables]
  [./xi]
    family = LAGRANGE
    order  = FIRST
  [../]

  [./theta_m]
    family = LAGRANGE
    order  = FIRST
  [../]
[]


###############################################################
#  Dampers
#  BoundingValueNodalDamper keeps xi in a safe range and
#  prevents Newton instability when gradients steepen.
###############################################################
[Dampers]
  [./limit_xi]
    type         = BoundingValueNodalDamper
    variable     = xi
    max_value    = 1.2
    min_value    = -0.2
    min_damping  = 1e-8   # allow very small damping factors
  [../]
[]


###############################################################
# Initial Conditions
#
#  xi:
#    - In the electrode: SmoothCircleIC with a circular void
#      (xi ≈ 0 inside, 1 outside).
#    - In the electrolyte: constant xi = 0 (void-like), but
#      no AC Kernels are applied there, so xi remains fixed.
#
#  theta_m:
#    - Initially uniform 0.9 in both electrode and electrolyte.
#      Only TimeDerivative kernel is applied (no flux/source),
#      so theta_m stays constant in time in this test.
###############################################################
[ICs]
  # Phase-field: diffuse circular void inside the electrode
  [./xi_ic1]
    type      = SmoothCircleIC
    variable  = xi
    x1        = 50e-6      # [m] void center x
    y1        = 50e-6      # [m] void center y
    radius    = 15e-6      # [m] void radius
    int_width = 1e-6       # [m] interface half width (diffuse layer)
    invalue   = 0.0        # void
    outvalue  = 1.0        # solid
    block     = 'electrode'
  [../]

  # xi in the electrolyte: fixed to 0.0 (here treated as void)
  # Note: no Kernels for xi in 'electrolyte', so xi there does
  #       not evolve in time.
  [./xi_ic2]
    type      = ConstantIC
    variable  = xi
    value     = 0.0
    block     = 'electrolyte'
  [../]

  # Li occupancy: initially uniform, here only a dummy variable
  [./theta_m_ic]
    type      = ConstantIC
    variable  = theta_m
    value     = 0.9
    block     = 'electrode electrolyte'
  [../]
[]


###############################################################
#  Kernels
#
#  1) theta_m equation (dummy):
#     - Only a TimeDerivative term, so theta_m stays constant.
#
#  2) xi equation (Allen–Cahn) in 'electrode' block:
#
#     ∂xi/∂t = -L * ( d g(xi)/d xi - kappa ∇² xi )
#
#     where:
#       g(xi) = w xi² (1 - xi)²
#
#     AC_time   → TimeDerivative term for xi
#     AC_bulk   → -L g'(xi)  (double-well chemical term)
#     AC_grad   → +L kappa ∇²xi  (gradient energy term)
###############################################################
[Kernels]

  # theta_m: pure time derivative (no diffusion / no source)
  # This simply defines a mass matrix, theta_m is constant in time.
  [./theta_time]
    type     = TimeDerivative
    variable = theta_m
    block    = 'electrode electrolyte'
  [../]

  # xi: time derivative ∂xi/∂t
  [./AC_time]
    type     = TimeDerivative
    variable = xi
    block    = 'electrode'
  [../]

  # xi: double-well derivative term -L g'(xi)
  [./AC_bulk]
    type     = AllenCahn
    variable = xi
    f_name   = g_xi        # material property g and g'
    mob_name = L           # mobility L
    block    = 'electrode'
  [../]

  # xi: gradient term L kappa ∇²xi
  [./AC_grad]
    type       = ACInterface
    variable   = xi
    kappa_name = kappa
    mob_name   = L
    block      = 'electrode'
  [../]

[]


###############################################################
#  Materials
#  L      : mobility in Allen–Cahn equation [m^2/(N·s)]
#  kappa  : gradient energy coefficient [N]
#  w      : double-well height [Pa]
#
#  g(xi) = w xi^2 (1 - xi)^2
#  DerivativeParsedMaterial provides g, g', g'' as g_xi.*
###############################################################
[Materials]

  # Constant phase-field parameters
  [./pf_constants]
    type        = GenericConstantMaterial
    prop_names  = 'L   kappa    w'
    prop_values = '1e-6 4.5e-7 3.5e6'
  [../]

  # Double-well potential g(xi) and its derivatives
  [./double_well]
    type                      = DerivativeParsedMaterial
    coupled_variables         = 'xi'
    material_property_names   = 'w'
    expression                = 'w * xi^2 * (1 - xi)^2'
    property_name             = g_xi
    derivative_order          = 2      # provide g, g', g''
    block                     = 'electrode'
    outputs                   = my_exodus
  [../]

[]


###############################################################
# Postprocessors (useful debug information)
###############################################################
[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]

  [./dt]
    type = TimestepSize
  [../]

  [./run_time]
    type         = PerfGraphData
    section_name = 'Root'
    data_type    = total
  [../]
[]


###############################################################
# Executioner
#  - Transient AC evolution
#  - Adaptive timestep to maintain Newton convergence
#  - PJFNK with hypre/BoomerAMG preconditioner
###############################################################
[Executioner]
  type        = Transient
  scheme      = bdf2
  solve_type  = PJFNK

  # Linear solver parameters
  petsc_options_iname  = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value  = 'hypre boomeramg 31 0.7'

  l_tol      = 1e-4
  l_max_its  = 30
  nl_max_its = 25
  nl_rel_tol = 1e-7

  end_time = 250.0
  dtmax    = 50.0
  dtmin    = 1e-8

  [./TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-1         # [s] initial timestep
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 8
  [../]
[]


###############################################################
# Outputs
###############################################################
[Outputs]
  [./my_checkpoint]
    type                  = Checkpoint
    file_base             = ./${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    additional_execute_on = 'INITIAL FINAL'
  [../]

  [./my_exodus]
    type                  = Nemesis
    file_base             = ./ex_${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]

  [./csv]
    type      = CSV
    file_base = ./csv_${my_filename}/out_${my_filename}
  [../]

  [./pgraph]
    type               = PerfGraphOutput
    time_step_interval = ${my_interval}
    level              = 2
    heaviest_branch    = true
    heaviest_sections  = 5
    execute_on         = 'TIMESTEP_END FINAL'
  [../]
[]
