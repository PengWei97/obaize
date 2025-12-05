###############################################
# Phase-field evolution of a circular void
# in a Li metal / solid electrolyte system.
#
# Equations:
#   - xi      : Allen–Cahn phase-field equation in the electrode
#               with double-well + gradient energy + chemical
#               contribution from Li occupancy.
#   - theta_m : Li site occupancy, evolving via a simplified
#               diffusion equation and a coupling to xi.
#
# This file is intended as a test case for combined
# interface-driven and chemical-potential-driven evolution
# of a void within the Li electrode, with a simple
# dummy-like theta_m that will be extended later.
#
# Author: Wei Peng
# Date  : 2025-12-06
###############################################

# Output naming
my_filename = s3_c2_pf_df_with_void1
my_interval = 20
my_theta_m0 = 0.999999998

###############################################################
# Mesh
#   - 2D square: 100 μm × 100 μm
#   - Resolution: 250 × 250 QUAD4 elements
#   - Block "electrode"   : 0–50 μm  (Li metal)
#   - Block "electrolyte" : 50–100 μm (solid electrolyte)
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
# Variables
#
#  xi      : phase-field order parameter
#            (xi = 1: solid Li, xi = 0: void)
#
#  theta_m : Li site occupancy (dimensionless). In this
#            test it is primarily used to provide a
#            chemical driving term and a simple diffusion
#            model; it can be extended in future work.
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
# Dampers
#   BoundingValueNodalDamper keeps xi and theta_m within a
#   reasonable range and stabilizes the Newton iterations.
###############################################################
[Dampers]
  [./limit_xi]
    type         = BoundingValueNodalDamper
    variable     = xi
    max_value    = 1.2
    min_value    = -0.2
    min_damping  = 1e-8   # allow very small damping factors
  [../]

  [./limit_theta_m]
    type         = BoundingValueNodalDamper
    variable     = theta_m
    max_value    = 1.2
    min_value    = -0.2
    min_damping  = 1e-8
  [../]
[]


###############################################################
# Initial Conditions
#
#  xi:
#    - In the electrode: SmoothCircleIC with a circular void
#      (xi ≈ 0 inside, xi ≈ 1 outside).
#    - In the electrolyte: constant xi = 0 (void-like). No
#      xi kernels are applied there, so xi remains fixed.
#
#  theta_m:
#    - Initially uniform my_theta_m0 in both electrode and
#      electrolyte. In this file, theta_m obeys a simple
#      diffusion equation coupled weakly to xi.
###############################################################
[ICs]
  # xi: diffuse circular void inside the electrode
  [./ic_xi_void_electrode]
    type      = SmoothCircleIC
    variable  = xi
    x1        = 50e-6      # [m] void center x
    y1        = 50e-6      # [m] void center y
    radius    = 15e-6      # [m] void radius
    int_width = 1e-6       # [m] interface half width
    invalue   = 0.0        # void
    outvalue  = 1.0        # solid Li
    block     = 'electrode'
  [../]

  # xi in the electrolyte: fixed to 0.0 (void-like)
  # Note: no xi kernels act in 'electrolyte', so xi remains 0 there.
  [./ic_xi_electrolyte_void]
    type      = ConstantIC
    variable  = xi
    value     = 0.0
    block     = 'electrolyte'
  [../]

  # theta_m: initially uniform in both electrode and electrolyte
  [./ic_theta_uniform]
    type      = ConstantIC
    variable  = theta_m
    value     = ${my_theta_m0}
    block     = 'electrode electrolyte'
  [../]
[]


###############################################################
# Kernels
#
#  (A) xi equation (Allen–Cahn) in the 'electrode' block:
#
#      ∂xi/∂t = -L [ g'(xi) - kappa ∇²xi
#                    + (R T / Ω_L) h'(xi)
#                      ln((1 - theta_m)/(1 - theta_m0)) ]
#
#      Implemented as:
#        - xi_time       : TimeDerivative (∂xi/∂t)
#        - xi_chem_log   : chemical driving term from theta_m
#        - xi_bulk_dw    : -L g'(xi) (double-well potential)
#        - xi_grad       : +L kappa ∇²xi (gradient energy)
#
#  (B) theta_m equation (simple diffusion + coupling):
#
#      h(xi) ∂theta_m/∂t + theta_m h'(xi) ∂xi/∂t
#        = ∇·( D(theta_m, xi) ∇theta_m )
#
#      Implemented as:
#        - theta_mass          : h(xi) ∂theta_m/∂t
#        - theta_xi_coupled    : theta_m h'(xi) ∂xi/∂t
#        - theta_diffusion     : ∇·(conc_diffusivity ∇theta_m)
###############################################################
[Kernels]
  # ---- (A1) xi time derivative --------------------------------
  [./xi_time]
    type     = TimeDerivative
    variable = xi
    block    = 'electrode'
  [../]

  # ---- (A2) xi chemical driving term from Li occupancy --------
  #   Contribution to xi-equation:
  #     -(L R T / Ω_L) h'(xi) ln[(1 - theta_m)/(1 - theta_m0)]
  [./xi_chem_log]
    type     = LiVoidAllenCahn
    variable = xi
    mob_name = L                # mobility L [m^2/(N·s)]

    theta_m  = theta_m
    h_name   = h_xi             # interpolation function h(xi)

    R        = 8.314            # gas constant [J/(mol·K)]
    T        = 298.0            # temperature [K]
    Omega_L  = 13.1e-6          # reference molar volume Ω_L [m^3/mol]
    theta_m0 = ${my_theta_m0}   # reference occupancy
    block    = 'electrode'
  [../]

  # ---- (A3) xi double-well term: -L g'(xi) --------------------
  [./xi_bulk_dw]
    type     = AllenCahn
    variable = xi
    f_name   = g_xi             # material property g(xi) and g'(xi)
    mob_name = L                # mobility L
    block    = 'electrode'
  [../]

  # ---- (A4) xi gradient energy term: L kappa ∇²xi -------------
  [./xi_grad]
    type       = ACInterface
    variable   = xi
    kappa_name = kappa
    mob_name   = L
    block      = 'electrode'
  [../]

  # ---------- (B) theta_m equation: mass + coupling + diffusion ----------

  # (B1) Mass term: h(xi) ∂theta_m/∂t
  [./theta_mass]
    type           = MatTimeDerivative
    variable       = theta_m
    mat_prop_value = h_xi
    block          = 'electrode electrolyte'
  [../] 

  # (B2) Coupled time term: theta_m h'(xi) ∂xi/∂t
  [./theta_xi_coupled]
    type           = MatCoupledTimeDerivative
    variable       = theta_m
    v              = xi
    mat_prop_value = theta_m_times_hprime
    block          = 'electrode'
  [../]

  # (B3) Diffusion term: ∇·( conc_diffusivity ∇theta_m )
  [./theta_diffusion]
    type        = MatDiffusion
    variable    = theta_m
    diffusivity = conc_diffusivity        # [m^2/s]
    block       = 'electrode'
  [../]
[]


###############################################################
# Materials
#
#  Parameters:
#    - D_eff  : reference diffusivity [m^2/s]
#    - L      : Allen–Cahn mobility [m^2/(N·s)]
#    - kappa  : gradient-energy coefficient [N]
#    - w      : double-well height [Pa]
#
#  Double-well potential:
#    g(xi) = w xi^2 (1 - xi)^2
#
#  Interpolation function:
#    h(xi) = xi^2 (xi^2 - 3 xi + 3) + 1e-8
#    (the 1e-8 term regularizes h(xi) to avoid degeneracy)
#
#  Effective diffusivity:
#    D(theta_m, xi) = D_eff h(xi) / max(1 - theta_m, 1e-3)
###############################################################
[Materials]

  # Constant phase-field and diffusion parameters
  [./pf_constants]
    type        = GenericConstantMaterial
    prop_names  = 'D_eff     L        kappa    w'
    prop_values = '7.5e-13  1e-6     4.5e-7  3.5e6'
  [../]

  # Double-well potential g(xi) and its derivatives
  [./pf_double_well]
    type                    = DerivativeParsedMaterial
    coupled_variables       = 'xi'
    material_property_names = 'w'
    expression              = 'w * xi^2 * (1 - xi)^2'
    property_name           = g_xi
    derivative_order        = 2      # provides g, g', g''
    block                   = 'electrode'
    outputs                 = out_exodus
  [../]

  # Interpolation function h(xi) (regularized with +1e-8)
  [./pf_interp_h]
    type              = ParsedMaterial
    coupled_variables = 'xi'
    expression        = 'xi*xi*(xi*xi-3*xi+3) + 1.0e-8'
    property_name     = h_xi
    block             = 'electrode electrolyte'
    outputs           = out_exodus
  [../]
    
  # Effective diffusivity: conc_diffusivity(theta_m, xi)
  [./pf_conc_diffusivity]
    type                    = DerivativeParsedMaterial
    coupled_variables       = 'theta_m'
    material_property_names = 'D_eff h_xi'
    expression              = 'D_eff*h_xi/max(1-theta_m, 1e-3)'
    property_name           = conc_diffusivity
    derivative_order        = 2
    block                   = 'electrode'
    outputs                 = out_exodus
  [../]

  # theta_m * h'(xi), where h'(xi) = xi (4 xi^2 - 9 xi + 6)
  [./pf_theta_times_hprime]
    type              = DerivativeParsedMaterial
    property_name     = theta_m_times_hprime
    coupled_variables = 'theta_m xi'
    expression        = 'theta_m*xi*(4*xi^2-9*xi+6)'
    derivative_order  = 2
    outputs           = out_exodus
    block             = 'electrode'    # can be extended to 'electrolyte' if needed
  [../]
[]


###############################################################
# Postprocessors
#   - dofs     : total number of DOFs
#   - dt       : current time step
#   - run_time : total wall-clock run time
###############################################################
[Postprocessors]
  [./pp_num_dofs]
    type = NumDOFs
  [../]

  [./pp_dt]
    type = TimestepSize
  [../]

  [./pp_run_time]
    type         = PerfGraphData
    section_name = 'Root'
    data_type    = total
  [../]
[]


###############################################################
# Executioner
#   - Transient BDF2 scheme
#   - Fully coupled solution via PJFNK
#   - Adaptive time-stepping based on nonlinear iterations
###############################################################
[Executioner]
  type       = Transient
  scheme     = bdf2
  solve_type = PJFNK

  # Linear solver (PETSc / Hypre BoomerAMG)
  petsc_options_iname  = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value  = 'hypre    boomeramg      31                  0.7'

  l_tol      = 1e-4
  l_max_its  = 30
  nl_max_its = 25
  nl_rel_tol = 1e-7

  end_time = 240.0
  dtmax    = 50.0
  dtmin    = 1e-8

  [./TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-1       # [s] initial timestep
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 8
  [../]
[]


###############################################################
# Outputs
###############################################################
[Outputs]
  [./out_checkpoint]
    type                  = Checkpoint
    file_base             = ./${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    additional_execute_on = 'INITIAL FINAL'
  [../]

  [./out_exodus]
    type                  = Nemesis
    file_base             = ./ex_${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]

  [./out_csv]
    type      = CSV
    file_base = ./csv_${my_filename}/out_${my_filename}
  [../]

  [./out_perfgraph]
    type               = PerfGraphOutput
    time_step_interval = ${my_interval}
    level              = 2
    heaviest_branch    = true
    heaviest_sections  = 5
    execute_on         = 'TIMESTEP_END FINAL'
  [../]
[]
