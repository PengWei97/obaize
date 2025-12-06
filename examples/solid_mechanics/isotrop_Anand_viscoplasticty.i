# =====================================================================
# 2D tensile test with isotropic Anand-type viscoplasticity
# implemented via AnandViscoplasticityStressUpdate.
#
# Features:
#   - Small-strain formulation with incremental strains
#   - J2 isotropic elasticity + Anand viscoplastic flow law
#   - Temperature-dependent viscosity term: exp(-Q/(R*T)), here T is constant
#   - Viscoplastic strain tensor stored as Material property
#     "viscoplastic_strain" by the StressUpdate class
# =====================================================================

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 10
  ymax = 10
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# ---------------------------------------------------------------------
# Boundary loading: prescribe vertical displacement on the top boundary.
# ---------------------------------------------------------------------
[Functions]
  [./top_pull]
    type = ParsedFunction
    expression = '0.1*t'         # Linear vertical displacement over time
  [../]
[]

# ---------------------------------------------------------------------
# Solid mechanics physics:
#   - Small strain
#   - Incremental formulation (required for radial return)
#   - Convenience "generate_output" to create AuxVariables for
#     selected strain/stress components and viscoplastic strain.
# ---------------------------------------------------------------------
[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain = SMALL
    incremental = true
    add_variables = true
    # Names here correspond to built-in outputs and Material tensor
    # components; "viscoplastic_strain_*" is from our StressUpdate class.
    generate_output = 'strain_yy stress_yy vonmises_stress' #  viscoplastic_strain_xx viscoplastic_strain_yy viscoplastic_strain_xy
  [../]
[]

# ---------------------------------------------------------------------
# Boundary conditions:
#   - Top: prescribed vertical displacement
#   - Left: fixed in x-direction (remove rigid-body translation)
#   - Bottom: fixed in y-direction (remove rigid-body translation)
# ---------------------------------------------------------------------
[BCs]
  [./y_pull_function]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = top_pull
  [../]

  [./x_fixed_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./y_fixed_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

# ---------------------------------------------------------------------
# Materials:
#   (1) Linear isotropic elasticity
#   (2) Anand viscoplastic model via AnandViscoplasticityStressUpdate
#   (3) ADComputeMultipleInelasticStress performs the radial return and
#       combines elastic and viscoplastic strain increments.
# ---------------------------------------------------------------------
[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 4.9e9      # 4.9 GPa for Li metal
    poissons_ratio = 0.38
  [../]

  [./anand_vp]
    type = AnandViscoplasticityStressUpdate

    # ---------- Visco-plastic parameters (from Table 2) ----------
    A  = 4.25e4                 # [1/s]
    Q  = 3.7e4                  # 37 kJ/mol -> 3.7e4 J/mol
    gas_constant = 8.314        # [J/(mol*K)]
    m  = 0.15

    S0  = 2.0e6                 # 2 MPa   -> 2.0e6 Pa (saturation coefficient)
    Sa0 = 1.1e6                 # Sa(t=0) -> 1.1e6 Pa (initial flow resistance)

    H0    = 1.0e7               # 10 MPa  -> 1.0e7 Pa
    a_exp = 2.0                 # hardening sensitivity a [-]
    n     = 0.05                # deformation resistance sensitivity n [-]

    # Constant absolute temperature T [K]
    T_constant = 298.0

    # Activation ratio: viscoplasticity active if sigma_eff > yield_ratio * Sa_old
    yield_ratio = 0.1

    absolute_tolerance = 1e-8
  [../]

  [./stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'anand_vp'
    tangent_operator = nonlinear      # 强制使用一致切线
    perform_finite_strain_rotations = false   # small-strain incremental formulation
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

  # [./vp_strain_yy]
  #   type = ElementAverageValue
  #   variable = viscoplastic_strain_yy
  # [../]

  # [./vp_strain_xy]
  #   type = ElementAverageValue
  #   variable = viscoplastic_strain_xy
  # [../]

  # # Effective viscoplastic strain scalar stored by RadialReturnStressUpdate
  # [./eff_vp_strain]
  #   type = ElementAverageMaterialProperty
  #   property = effective_viscoplastic_strain
  # [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

# ---------------------------------------------------------------------
# Time integration and nonlinear solution:
#   - Transient quasi-static solve with PJFNK
#   - Adaptive time stepping based on nonlinear iteration count
# ---------------------------------------------------------------------
[Executioner]
  # type = Transient
  # solve_type = 'NEWTON'       # instead of PJFNK

  # l_max_its = 80              # NEWTON 下通常线性迭代会少很多
  # l_tol     = 1e-5

  # nl_max_its = 15
  # nl_rel_tol = 1e-6           # 可以放松一点
  # nl_abs_tol = 1e-8

  type = Transient
  solve_type = PJFNK

  l_max_its = 80
  l_tol     = 1e-6

  nl_max_its = 20
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  # num_steps = 10
  end_time   = 1.0
  dtmin      = 1e-4

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5.0e-3
    growth_factor = 1.3
    cutback_factor = 0.5
    optimal_iterations = 6
  [../]
[]

# ---------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------
[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = true
  perf_graph = true
[]