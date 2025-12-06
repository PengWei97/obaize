###############################################
# Phase-field shrinking of a circular void
# Only AC equation: double-well + gradient energy
# No Li kinetics, no electrochemistry, no mechanics
#
# Goal:
#   Test whether a void evolves under AC curvature-driven dynamics.
#
# Author: Wei Peng
# Date: 2025.12.06
###############################################

# Output naming
my_filename = s2_c1_pf_with_void
my_interval = 20


###############################################################
#  Mesh
#  - 2D square: 100 μm × 100 μm
#  - Resolution: 250 × 250 QUAD4 elements
###############################################################
[Mesh]
  [./mesh]
    type   = GeneratedMeshGenerator
    dim    = 2
    nx     = 250
    ny     = 250
    xmin   = 0.0
    xmax   = 100e-6     # [m]
    ymin   = 0.0
    ymax   = 100e-6     # [m]
    elem_type = QUAD4
  [../]
[]
  

###############################################################
#  Variables
#  xi : phase-field order parameter
#       1 = solid (Li), 0 = void
###############################################################
[Variables]
  [./xi]
    family = LAGRANGE
    order  = FIRST
  [../]
[]

[AuxVariables]
  [./bnds]
  [../]
[]

[AuxKernels]
  [./bnds]
    type = BndsCalcAux
    variable = bnds
    v = 'xi'
  [../]
[]

###############################################################
#  Dampers
#  BoundingValueNodalDamper keeps xi in a safe range and
#  prevents Newton instability when gradients steepen.
###############################################################
[Dampers]
  [./limit_xi]
    type      = BoundingValueNodalDamper
    variable  = xi
    max_value = 1.2
    min_value = -0.2
    min_damping = 1e-8     # avoid rejecting steps too aggressively
  [../]
[]


###############################################################
# Initial Conditions
# - SmoothCircleIC creates a diffuse circular void.
# - int_width controls initial interface thickness.
###############################################################
[ICs]
  [./xi_ic]
    type      = SmoothCircleIC
    variable  = xi
    x1        = 50e-6      # [m] void center x
    y1        = 50e-6      # [m] void center y
    radius    = 15e-6      # [m] void radius
    int_width = 1e-6       # interface half width (diffuse layer)
    invalue   = 0.0        # void
    outvalue  = 1.0        # solid
  [../]
[]


###############################################################
#  Kernels: Allen–Cahn phase-field equation
#
#  ∂xi/∂t = -L * ( d(g)/d(xi)  - kappa ∇²xi )
#
#  Li_AC_time   → TimeDerivative term
#  Li_AC_bulk2  → -L g'(xi)  (double-well chemical term)
#  Li_AC_grad   → +L kappa ∇²xi  (gradient energy)
###############################################################
[Kernels]

  # Time derivative term ∂xi/∂t
  [./AC_time]
    type     = TimeDerivative
    variable = xi
  [../]

  # Double-well derivative term: -L g'(xi)
  [./AC_bulk]
    type     = AllenCahn
    variable = xi
    f_name   = g_xi        # material property g and g'
    mob_name = L           # mobility L
  [../]

  # Gradient term: L kappa ∇²xi
  [./AC_grad]
    type       = ACInterface
    variable   = xi
    kappa_name = kappa
    mob_name   = L
  [../]

[]


###############################################################
#  Materials
#  L      : mobility in Allen–Cahn equation
#  kappa  : gradient energy coefficient
#  w      : double-well height
#
#  g(xi) = w xi² (1 - xi)²
#  DerivativeParsedMaterial automatically generates g', g''
###############################################################
[Materials]

  # Constant parameters
  [./pf_constants]
    type        = GenericConstantMaterial
    prop_names  = 'L   kappa    w'
    prop_values = '1e-6 4.5e-7 3.5e6'
  [../]

  # Double-well potential g(xi)
  [./double_well]
    type                      = DerivativeParsedMaterial
    coupled_variables         = 'xi'
    material_property_names   = 'w'
    expression                = 'w * xi^2 * (1 - xi)^2'
    property_name             = g_xi
    derivative_order          = 2      # provide g, g', g''
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
    type        = PerfGraphData
    section_name = 'Root'
    data_type    = total
  [../]
[]

[Adaptivity]
  # initial_steps = 1
  # initial_marker = err_xi

  max_h_level = 1
  marker = err_xi

  [./Indicators]
    [./ind_xi]
      type = GradientJumpIndicator
      variable = xi
     [../]
    [./ind_bnds]
      type = GradientJumpIndicator
      variable = bnds
    [../]
  [../]
  [./Markers]
    [./err_xi]
      type = ErrorFractionMarker
      coarsen = 0.3
      refine = 0.95
      indicator = ind_xi
    [../]
    [./err_bnds]
      type = ErrorFractionMarker
      coarsen = 0.3
      refine = 0.95
      indicator = ind_bnds
    [../]
  [../]
[]

###############################################################
# Executioner
# - Adaptive timestep to maintain Newton convergence
# - PJFNK with hypre/BoolemerAMG
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
    dt                 = 1e-1         # initial timestep
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
