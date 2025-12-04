my_filename = s2_df_pf_Li
my_interval = 2

[Mesh]
  type      = GeneratedMesh
  dim       = 2
  nx        = 100
  ny        = 100
  xmax      = 100.0e-6   # 100 μm
  ymax      = 100.0e-6   # 100 μm
  elem_type = QUAD4
[]

[Variables]
  # Li lattice site fraction (0–1)
  [./theta_m]
    family = LAGRANGE
    order  = FIRST
  [../]

  # Phase-field order parameter (metal = 1, void = 0)
  [./xi]
    family = LAGRANGE
    order  = FIRST
  [../]
[]

# Optional bound dampers for stability (currently disabled)
[Dampers]
  [./limit_theta_m]
    type      = BoundingValueNodalDamper
    variable  = theta_m
    max_value = 0.9999999
    min_value = 1e-8
  [../]

  [./limit_xi]
    type      = BoundingValueNodalDamper
    variable  = xi
    max_value = 1.0
    min_value = 0.0
  [../]
[]

# ----------------------------------------------------
# Initial Conditions
#   Domain: 100 x 100 μm
#   Void:   circle, diameter 30 μm (radius 15 μm),
#           centered at (50 μm, 50 μm)
#
#   xi:
#     inside void (r < 15 μm): xi ≈ 0 (void)
#     outside                 : xi ≈ 1 (solid)
#
#   theta_m:
#     inside void             : 0.9 (Li-poor)
#     outside                 : 0.9 (Li-rich)
# ----------------------------------------------------
[ICs]
  [./xi_ic]
    type      = SmoothCircleIC
    variable  = xi
    x1        = 50.0e-6
    y1        = 50.0e-6
    radius    = 15.0e-6
    int_width = 1.0e-6
    invalue   = 0.0      # void
    outvalue  = 1.0      # solid
  [../]

  [./theta_m_ic]
    type      = SmoothCircleIC
    variable  = theta_m
    x1        = 50.0e-6
    y1        = 50.0e-6
    radius    = 15.0e-6
    int_width = 1.0e-6   # sharper interface for concentration
    invalue   = 0.999999998      # void: low Li
    outvalue  = 0.999999998      # solid: high Li
  [../]
[]

# ----------------------------------------------------
# Boundary Conditions
#
# For MatDiffusion and Allen–Cahn, the natural BC is zero flux:
#   -j · n = 0 for theta_m and xi on all boundaries.
# So we can omit explicit BC blocks here.
#
# If desired, explicit zero-flux BC could be added as:
# [BCs]
#   [./theta_no_flux]
#     type     = NeumannBC
#     variable = theta_m
#     boundary = 'left right top bottom'
#     value    = 0.0
#   [../]
# []

# ----------------------------------------------------
# Kernels
# ----------------------------------------------------
[Kernels]
  # --------------------
  # theta_m equation:
  #   h(xi) * d(theta_m)/dt
  # + theta_m * h'(xi) * d(xi)/dt
  # = div( D_eff(theta_m, xi) * grad(theta_m) )
  # --------------------

  # Time derivative term: h(xi) * theta_m_dot
  [./Li_df_time1]
    type            = MatTimeDerivative
    variable        = theta_m
    mat_prop_value  = h_xi
  [../]

  # Coupled time term: theta_m * h'(xi) * xi_dot
  [./Li_df_time2]
    type            = MatCoupledTimeDerivative
    variable        = theta_m
    v               = xi
    mat_prop_value  = theta_m_times_hprime
  [../]

  # Diffusion term: div( conc_diffusivity * grad(theta_m) )
  [./Li_diffusion]
    type        = MatDiffusion
    variable    = theta_m
    diffusivity = conc_diffusivity
  [../]

  # --------------------
  # xi equation (Allen–Cahn type):
  #   d(xi)/dt = -L * dF/dxi
  # including chemical free energy, double-well and gradient terms.
  # --------------------

  # Time derivative term: d(xi)/dt
  [./Li_AC_time]
    type     = TimeDerivative
    variable = xi
  [../]

  # Chemical contribution:
  #   -(L * R * T / Omega_L) * h'(xi) * ln[(1 - theta_m)/(1 - theta_m0)]
  [./Li_AC_bulk1]
    type     = LiVoidAllenCahn
    variable = xi
    mob_name = L

    theta_m  = theta_m
    h_name   = h_xi

    R        = 8.314       # J/(mol*K)
    T        = 298.0       # K
    Omega_L  = 13.1e-6     # m^3/mol
    theta_m0 = 0.999999998 # reference occupancy
  [../]

  # Double-well contribution: -L * g'(xi)
  [./Li_AC_bulk2]
    type     = AllenCahn
    variable = xi
    f_name   = g_xi
    mob_name = L
  [../]

  # Gradient contribution: L * kappa * Laplacian(xi)
  [./Li_AC_grad]
    type       = ACInterface
    variable   = xi
    kappa_name = kappa
    mob_name   = L
  [../]
[]

# ----------------------------------------------------
# Materials
# ----------------------------------------------------
[Materials]
  # Constant parameters: effective diffusivity, mobility, interface width, double-well height
  [./Material_prop_const]
    type        = GenericConstantMaterial
    prop_names  = 'D_eff     L      kappa    w'
    prop_values = '7.5e-13  1e-9   4.5e-7  3.5e6'
  [../]

  # Interpolation function h(xi) (phase-field switching)
  # h(xi) = xi^2 (xi^2 - 3 xi + 3)
  [./interpore_func]
    type              = DerivativeParsedMaterial
    block             = 0
    coupled_variables = 'xi'
    expression        = 'xi*xi*(xi*xi-3*xi+3)'
    property_name     = h_xi
    derivative_order  = 2   # h_xi, dh/dxi, d2h/dxi2
  [../]

  # Double-well potential:
  #   g(xi) = w * xi^2 * (1 - xi)^2
  [./double_well_func]
    type                      = DerivativeParsedMaterial
    block                     = 0
    coupled_variables         = 'xi'
    material_property_names   = 'w'
    expression                = 'w*xi^2*(1-xi)^2'
    property_name             = g_xi
    derivative_order          = 2   # g_xi, dg/dxi, d2g/dxi2
  [../]

  # Concentration-dependent effective diffusivity:
  #   D(theta_m, xi) = D_eff * h(xi) / max(1 - theta_m, 1e-3)
  [./conc_diffusivity]
    type                    = ParsedMaterial
    block                   = 0
    coupled_variables       = 'theta_m'
    material_property_names = 'D_eff h_xi'
    expression              = 'D_eff*h_xi/max(1-theta_m, 1e-3)'
    property_name           = conc_diffusivity
    outputs                 = my_exodus
  [../]

  # theta_m * h'(xi), where h'(xi) = xi * (4 xi^2 - 9 xi + 6)
  # A small offset is added to avoid exactly zero in some regions.
  [./theta_m_times_hprime]
    type              = ParsedMaterial
    property_name     = theta_m_times_hprime
    coupled_variables = 'theta_m xi'
    expression        = 'theta_m*xi*(4*xi^2-9*xi+6)'
  [../]
[]

# ----------------------------------------------------
# Preconditioning
# ----------------------------------------------------
[Preconditioning]
  [./pc]
    type = SMP
    full = true
  [../]
[]

# ----------------------------------------------------
# Executioner
# ----------------------------------------------------
[Executioner]
  type       = Transient
  solve_type = NEWTON
  scheme     = 'bdf2'
  line_search = basic

  # PETSc solver options: GMRES + hypre BoomerAMG
  petsc_options_iname  = '-pc_type -ksp_type'
  petsc_options_value  = 'hypre    gmres'

  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10

  l_max_its  = 80
  l_tol      = 1e-8

  # end_time = 300.0
  num_steps = 20
  dt       = 1.0e-3   # initial time step
  dtmax    = 1.0      # maximum time step

  [./TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1.0e-3
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 8
  [../]
[]

# ----------------------------------------------------
# Outputs
# ----------------------------------------------------
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