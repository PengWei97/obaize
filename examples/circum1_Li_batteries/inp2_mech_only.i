my_filename = s2_mech_only
my_interval = 2

my_theta_m0 = 0.999999998
my_bulk_modulus  = 6.81e9   # Pa
my_shear_modulus = 1.78e9   # Pa

[Mesh]
  type      = GeneratedMesh
  dim       = 2
  nx        = 100
  ny        = 100
  xmax      = 100.0e-6   # 100 μm = 1.0e-4 m
  ymax      = 100.0e-6   # 100 μm = 1.0e-4 m
  elem_type = QUAD4
[]

# ------------------------------------------------------------
# Global parameters: displacements for SolidMechanics
# ------------------------------------------------------------
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./sigma_h_xi]
    family = MONOMIAL   # elemental auxiliary variable
    order  = FIRST      # piecewise-linear inside each element
  [../]

  # Li lattice site fraction (0–1)
  [./theta_m]
    family = MONOMIAL
    order  = FIRST
  [../]

  # Phase-field order parameter (metal = 1, void = 0)
  [./xi]
    family = MONOMIAL
    order  = FIRST
  [../]
[]

# ------------------------------------------------------------
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
#     currently set to a uniform value my_theta_m0 everywhere
#     (change invalue/outvalue below if you want void/solid contrast)
# ------------------------------------------------------------
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
    invalue   = ${my_theta_m0}  # void region (currently same as solid)
    outvalue  = ${my_theta_m0}  # solid region
  [../]
[]

# ------------------------------------------------------------
# Functions: prescribed vertical displacement on the top boundary
# ------------------------------------------------------------
[Functions]
  [./left_press]
    type       = ParsedFunction
    expression = '2.0e-2*t*100.0e-6'   # linear vertical displacement over time (m)
  [../]
[]

# ------------------------------------------------------------
# Boundary conditions: solid mechanics
# ------------------------------------------------------------
# ------------------------------------------------------------
# [BCs]
#   [./y_displacement_left]
#     type     = FunctionDirichletBC
#     variable = disp_x
#     boundary = left
#     function = left_press
#   [../]

#   [./x_fixed_right]
#     type     = DirichletBC
#     variable = disp_x
#     boundary = right
#     value    = 0.0
#   [../]

#   [./y_fixed_bottom_top]
#     type     = DirichletBC
#     variable = disp_y
#     boundary = 'bottom top'
#     value    = 0.0
#   [../]
# []

[BCs]
  [./y_displacement_top]
    type     = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = left_press
  [../]

  [./x_fixed_left]
    type     = DirichletBC
    variable = disp_x
    boundary = left
    value    = 0.0
  [../]

  [./y_fixed_bottom]
    type     = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value    = 0.0
  [../]
[]

# ------------------------------------------------------------
# Solid mechanics physics block
# ------------------------------------------------------------
[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain          = SMALL
    incremental     = true
    add_variables   = true          # creates disp_x, disp_y

    # Register eigenstrain contribution from chemical expansion
    eigenstrain_names = 'chem_strain'

    generate_output = 'strain_xx strain_yy strain_xy 
                       stress_xx stress_yy stress_xy
                       vonmises_stress'
  [../]
[]

[AuxKernels]
  [./sigma_h_xi]
    type           = HydrostaticStressXiAux
    variable       = sigma_h_xi
    elastic_strain = elastic_strain     # from solid mechanics
    h_xi           = h_xi               # interpolation function
    bulk_modulus   = ${my_bulk_modulus} # K in Pa (from E, nu)
  [../]
[]

# ------------------------------------------------------------
# Materials
# ------------------------------------------------------------
[Materials]
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

  # # Double-well potential:
  # #   g(xi) = w * xi^2 * (1 - xi)^2
  # [./double_well_func]
  #   type                      = DerivativeParsedMaterial
  #   block                     = 0
  #   coupled_variables         = 'xi'
  #   material_property_names   = 'w'
  #   expression                = 'w*xi^2*(1-xi)^2'
  #   property_name             = g_xi
  #   derivative_order          = 2   # g_xi, dg/dxi, d2g/dxi2
  # [../]

  # Linear isotropic elasticity with phase-field interpolation
  [./elasticity_tensor]
    type          = ComputePFIsotropicElasticityTensor
    bulk_modulus  = ${my_bulk_modulus}  # Pa
    shear_modulus = ${my_shear_modulus} # Pa
    h_xi          = h_xi
  [../]

  # Chemical eigenstrain due to Li concentration (Ω in m^3/mol)
  [./chem_eigen]
    type             = ChemicalEigenstrain
    eigenstrain_name = chem_strain

    xi       = xi
    theta_m  = theta_m

    theta_m0  = ${my_theta_m0}
    Omega_L   = 13.1e-6
    Omega_v   = 6.0e-6
    h_xi_name = h_xi
  [../]   

  # Anand viscoplasticity
  # Units:
  #   A  [1/s], Q [J/mol], gas_constant [J/(mol·K)],
  #   stresses (S0, Sa0, H0) in Pa (consistent with elasticity),
  #   T_constant in K.
  [./anand_vp]
    type          = AnandViscoplasticityStressUpdate

    A  = 4.25e4                 # 1/s
    Q  = 37e3                   # J/mol
    gas_constant = 8.314        # J/(mol·K)
    m  = 0.15

    S0  = 2.0e6    # Pa, reference flow stress
    Sa0 = 1.1e6    # Pa, reference internal stress

    H0    = 1.0e7  # Pa, hardening modulus
    a_exp = 2.0
    n     = 0.05

    T_constant  = 298.0 # K
    yield_ratio = 0.1

    absolute_tolerance = 1e-8
  [../]

  # Stress update with inelastic (Anand) model
  [./radial_return_stress]
    type                     = ComputeMultipleInelasticStress
    tangent_operator         = elastic
    inelastic_models         = 'anand_vp'
    perform_finite_strain_rotations = false
  [../]
[]

# ---------------------------------------------------------------------
# Postprocessors:
#   - Average strain, stress and von Mises stress
#   - Average viscoplastic strain components
#   - Average effective viscoplastic strain (scalar, from StressUpdate)
# ---------------------------------------------------------------------
[Postprocessors]
  [./strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  [../]

  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]

  [./vonmises_stress]
    type = ElementAverageValue
    variable = vonmises_stress
  [../]
[]

# ------------------------------------------------------------
# Preconditioning
# ------------------------------------------------------------
[Preconditioning]
  [./pc]
    type = SMP
    full = true
  [../]
[]

# ------------------------------------------------------------
# Executioner
# ------------------------------------------------------------
[Executioner]
  type        = Transient
  solve_type  = NEWTON
  scheme      = 'bdf2'
  line_search = basic

  # PETSc solver options: GMRES + hypre BoomerAMG
  petsc_options_iname  = '-pc_type -ksp_type'
  petsc_options_value  = 'hypre    gmres'

  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10

  l_max_its  = 80
  l_tol      = 1e-8

  end_time = 0.1
  # num_steps = 20
  dt        = 1.0e-3   # s
  dtmax     = 1.0      # s

  [./TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1.0e-3
    growth_factor      = 1.3
    cutback_factor     = 0.7
    optimal_iterations = 8
  [../]
[]

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
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
