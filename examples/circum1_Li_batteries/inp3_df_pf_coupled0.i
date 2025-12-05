###############################################
# Fully coupled Li metal / solid electrolyte phase-field model
#
# Primary variables:
#   - theta_m : dimensionless Li lattice-site occupancy (0–1)
#   - xi      : phase-field order parameter (xi = 1: Li metal, xi = 0: void)
#
# Author : Wei Peng
# Date   : 2025-12-06
###############################################

# File naming and basic parameters
my_filename  = s3_c1_pf_df_with_void0
my_interval  = 20
my_theta_m0  = 0.999999998    # reference occupancy (close to Li-metal saturation)

###############################################################
# Mesh
#   2D square domain: 100 μm × 100 μm
#   Structured QUAD4 mesh: 250 × 250 elements
###############################################################
[Mesh]
  [./mesh]
    type       = GeneratedMeshGenerator
    dim        = 2
    nx         = 250
    ny         = 250
    xmin       = 0.0
    xmax       = 100e-6       # [m]
    ymin       = 0.0
    ymax       = 100e-6       # [m]
    elem_type  = QUAD4
  [../]
[]

###############################################################
# Variables
#   xi      : phase-field order parameter (1: solid Li, 0: void)
#   theta_m : Li lattice-site occupancy
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
#   reasonable range to stabilize Newton iterations.
###############################################################
[Dampers]
  [./limit_xi]
    type        = BoundingValueNodalDamper
    variable    = xi
    max_value   = 1.2
    min_value   = -0.2
    min_damping = 1e-8        # avoid rejecting steps too aggressively
  [../]

  [./limit_theta_m]
    type        = BoundingValueNodalDamper
    variable    = theta_m
    max_value   = 1.2
    min_value   = -0.2
    min_damping = 1e-8
  [../]
[]

###############################################################
# Initial Conditions
#   - xi      : diffuse circular void embedded in Li metal
#   - theta_m : uniform initial Li occupancy
###############################################################
[ICs]
  # Phase-field: diffuse circular void
  [./xi_ic]
    type      = SmoothCircleIC
    variable  = xi
    x1        = 50e-6        # [m] center x-coordinate
    y1        = 50e-6        # [m] center y-coordinate
    radius    = 15e-6        # [m] void radius
    int_width = 10e-6        # [m] half width of diffuse interface
    invalue   = 0.0          # void
    outvalue  = 1.0          # solid Li
  [../]

  # Li occupancy: uniform initial value
  [./theta_m_ic]
    type     = ConstantIC
    variable = theta_m
    value    = ${my_theta_m0}
  [../]
[]

###############################################################
# Kernels
#
# (A) Allen–Cahn phase-field equation for xi:
#
#   ∂xi/∂t = -L [ g'(xi) - kappa ∇²xi
#                 + (R T / Ω_L) h'(xi) ln((1 - theta_m)/(1 - theta_m0)) ]
#
# Implemented as:
#   - TimeDerivative   : ∂xi/∂t
#   - LiVoidAllenCahn  : L (R T / Ω_L) h'(xi) ln(...)
#   - AllenCahn        : -L g'(xi) (double-well potential)
#   - ACInterface      : +L kappa ∇²xi (gradient-energy term)
#
# (B) theta_m equation (mass balance with interpolation h(xi)):
#
#   h(xi) ∂theta_m/∂t + theta_m h'(xi) ∂xi/∂t
#     = ∇·( D(theta_m, xi) ∇theta_m )
#
# Implemented as:
#   - MatTimeDerivative         : h(xi) ∂theta_m/∂t
#   - MatCoupledTimeDerivative  : theta_m h'(xi) ∂xi/∂t
#   - MatDiffusion              : ∇·( D_eff h(xi)/max(1 - theta_m, 1e-3) ∇theta_m )
###############################################################
[Kernels]

  # ---- (A1) xi time derivative --------------------------------
  [./AC_time]
    type     = TimeDerivative
    variable = xi
  [../]

  # ---- (A2) Chemical term from Li occupancy (no mechanics) ----
  #   Contribution to xi-equation:
  #     -(L R T / Ω_L) h'(xi) ln[(1 - theta_m)/(1 - theta_m0)]
  [./Li_AC_bulk1]
    type     = LiVoidAllenCahn
    variable = xi
    mob_name = L                  # mobility L [m^2/(N·s)]

    theta_m  = theta_m
    h_name   = h_xi               # interpolation function h(xi)

    R        = 8.314              # gas constant [J/(mol·K)]
    T        = 298.0              # temperature [K]
    Omega_L  = 13.1e-6            # reference molar volume Ω_L [m^3/mol]
    theta_m0 = ${my_theta_m0}     # reference Li occupancy
  [../]

  # ---- (A3) Double-well derivative: -L g'(xi) -----------------
  [./AC_bulk]
    type     = AllenCahn
    variable = xi
    f_name   = g_xi               # material property g(xi) and g'(xi)
    mob_name = L                  # mobility L
  [../]

  # ---- (A4) Gradient term: L kappa ∇²xi -----------------------
  [./AC_grad]
    type       = ACInterface
    variable   = xi
    kappa_name = kappa            # gradient-energy coefficient
    mob_name   = L
  [../]

  # ---------- (B) theta_m equation: time + coupling + diffusion ----------

  # (B1) Mass term: h(xi) ∂theta_m/∂t
  [./theta_time]
    type           = MatTimeDerivative
    variable       = theta_m
    mat_prop_value = h_xi
  [../]

  # (B2) Coupled time term: theta_m h'(xi) ∂xi/∂t
  [./theta_xi_coupled_time]
    type           = MatCoupledTimeDerivative
    variable       = theta_m
    v              = xi
    mat_prop_value = theta_m_times_hprime
  [../]

  # (B3) Diffusion term:
  #   ∇·( conc_diffusivity ∇theta_m )
  [./theta_diffusion]
    type        = MatDiffusion
    variable    = theta_m
    diffusivity = conc_diffusivity    # [m^2/s]
  [../]
[]

###############################################################
# Materials
#
# Constants:
#   - D_eff  : reference diffusivity in Li metal
#   - L      : Allen–Cahn mobility
#   - kappa  : gradient-energy coefficient
#   - w      : double-well barrier height
#
# Double-well potential:
#   g(xi) = w xi^2 (1 - xi)^2
#
# Interpolation function:
#   h(xi) = xi^2 (xi^2 - 3 xi + 3) + 1e-8
#   (the 1e-8 term regularizes h(xi) to avoid degeneracy)
#
# Effective diffusivity:
#   D(theta_m, xi) = D_eff h(xi) / max(1 - theta_m, 1e-3)
###############################################################
[Materials]

  # Constant material parameters
  [./pf_constants]
    type        = GenericConstantMaterial
    prop_names  = 'D_eff     L        kappa    w       const_test'
    prop_values = '7.5e-13  1e-6     4.5e-7  3.5e6   1.0'
  [../]

  # Double-well potential g(xi)
  [./double_well]
    type                    = DerivativeParsedMaterial
    coupled_variables       = 'xi'
    material_property_names = 'w'
    expression              = 'w * xi^2 * (1 - xi)^2'
    property_name           = g_xi
    derivative_order        = 2           # provides g, g', g''
    outputs                 = my_exodus
  [../]

  # Interpolation function h(xi), regularized by +1e-8
  [./interpore_func]
    type              = ParsedMaterial
    coupled_variables = 'xi'
    expression        = 'xi*xi*(xi*xi-3*xi+3) + 1.0e-8'
    property_name     = h_xi
    outputs           = my_exodus
  [../]

  # Effective diffusivity:
  #   conc_diffusivity = D_eff * h(xi) / max(1 - theta_m, 1e-3)
  [./conc_diffusivity]
    type                    = DerivativeParsedMaterial
    coupled_variables       = 'theta_m'
    material_property_names = 'D_eff h_xi'
    expression              = 'D_eff*h_xi/max(1-theta_m, 1e-3)'
    property_name           = conc_diffusivity
    derivative_order        = 2
    outputs                 = my_exodus
  [../]

  # theta_m * h'(xi),
  # where h'(xi) = xi (4 xi^2 - 9 xi + 6)
  [./theta_m_times_hprime]
    type              = DerivativeParsedMaterial
    property_name     = theta_m_times_hprime
    coupled_variables = 'theta_m xi'
    expression        = 'theta_m*xi*(4*xi^2-9*xi+6)'
    derivative_order  = 2
    outputs           = my_exodus
  [../]
[]

###############################################################
# Postprocessors
#   - NumDOFs       : total number of DOFs
#   - TimestepSize  : current time step
#   - PerfGraphData : wall-clock runtime profiling
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
#   - Transient BDF2 scheme
#   - Fully coupled solve via PJFNK
#   - Adaptive time-stepping based on nonlinear iterations
###############################################################
[Executioner]
  type       = Transient
  scheme     = bdf2
  solve_type = PJFNK

  # Linear solver (PETSc/Hypre BoomerAMG)
  petsc_options_iname  = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value  = 'hypre    boomeramg      31                  0.7'

  l_tol      = 1e-4
  l_max_its  = 30
  nl_max_its = 25
  nl_rel_tol = 1e-7

  end_time = 250.0
  dtmax    = 50.0
  dtmin    = 1e-8

  [./TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-5       # initial timestep
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 8
  [../]
[]

###############################################################
# Outputs
###############################################################
[Outputs]
  # Checkpoint files
  [./my_checkpoint]
    type                  = Checkpoint
    file_base             = ./${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    additional_execute_on = 'INITIAL FINAL'
  [../]

  # Exodus/Nemesis output for visualization
  [./my_exodus]
    type                  = Nemesis
    file_base             = ./ex_${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]

  # CSV output for postprocessors
  [./csv]
    type      = CSV
    file_base = ./csv_${my_filename}/out_${my_filename}
  [../]

  # Performance profiling
  [./pgraph]
    type               = PerfGraphOutput
    time_step_interval = ${my_interval}
    level              = 2
    heaviest_branch    = true
    heaviest_sections  = 5
    execute_on         = 'TIMESTEP_END FINAL'
  [../]
[]
