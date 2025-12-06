###############################################
# TODO CASE: Fully-coupled diffusion–phase-field–electrochemistry–mechanics
#            for a circular void in Li metal / solid electrolyte
#
# Fields:
#   - xi      : Phase-field order parameter in the Li electrode
#               (xi = 1: solid Li, xi = 0: void/electrolyte).
#
#   - theta_m : Li site occupancy (dimensionless), providing a
#               chemical driving term for xi and obeying a
#               diffusion equation with mechano-chemical coupling.
#
#   - phi     : Electric potential, solving a conduction equation
#               with xi-dependent and material-dependent conductivity.
#
#   - disp_x, disp_y : Mechanical displacement components in 2D
#                      small-strain solid mechanics.
#
# Purpose:
#   This input aims at a fully coupled 4-field model:
#       diffusion (theta_m)
#       + phase-field (xi)
#       + electrochemistry (phi)
#       + mechanics (disp_x, disp_y)
#   for a Li metal / solid electrolyte system containing a circular void.
#
# IMPORTANT / TODO:
#   - When the right boundary applies a non-zero current via
#
#       [./phi_right_current]
#         type     = NeumannBC
#         variable = phi
#         boundary = 'right'
#         value    = -${my_i_app}
#       [../]
#
#     the strong conductivity jump between the Li metal electrode
#     and the solid electrolyte, combined with full coupling
#     (xi, theta_m, phi, mechanics), leads to convergence difficulties.
#   - This case has NOT been numerically validated yet; it is
#     a work-in-progress test that requires:
#       * mesh adaptivity near the interface and void,
#       * improved preconditioning / solver tuning,
#       * careful parameter calibration.
#
#   当前状态：本算例尚未通过系统性验证，在设置非零电流边界条件、
#   且存在电极/电解质导电率剧烈突变与四场强耦合时，容易导致不收敛。
#   后续需要通过网格自适应与数值算法优化来提高鲁棒性 —— 这是一个 TODO。
#
# Author: Wei Peng
# Date  : 2025-12-06
###############################################

# Output naming and basic parameters
my_filename = s5_fully_coupled
my_interval = 20
my_theta_m0 = 0.999999998

# Applied current density on the right boundary (counter electrode).
# NOTE:
#   - Non-zero my_i_app + large conductivity contrast
#     + strong multiphysics coupling => may cause convergence issues (TODO).
my_i_app = 1.0   # [A/m^2]

###############################################################
# Mesh
#
#   Geometry:
#     - 2D square: 100 μm × 100 μm
#     - Resolution: 250 × 250 QUAD4 elements
#
#   Subdomains:
#     - Block "electrode"   : 0 – 50 μm   (Li metal)
#     - Block "electrolyte" : 50 – 100 μm (solid electrolyte)
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

  # Subdomain 1: Li electrode (left half: 0 – 50 μm)
  [./electrode_block]
    type        = SubdomainBoundingBoxGenerator
    input       = mesh
    block_id    = 1
    block_name  = 'electrode'
    bottom_left = '0.0 0.0 0.0'
    top_right   = '50e-6 100e-6 0.0'
  [../]

  # Subdomain 2: solid electrolyte (right half: 50 – 100 μm)
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
# Global parameters
#
#   - displacements: names of mechanical displacement variables.
###############################################################
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

###############################################################
# Variables
#
#   xi       : Phase-field order parameter.
#   theta_m  : Li site occupancy (dimensionless).
#   phi      : Electric potential.
#   disp_x,y : Mechanical displacement (added automatically by
#              SolidMechanics/QuasiStatic with add_variables = true).
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

  [./phi]
    family = LAGRANGE
    order  = FIRST
  [../]
[]

###############################################################
# Auxiliary variables
#
#   i_mag      : |i| = |-sigma ∇phi|, current density magnitude.
#   sigma_h_xi : xi-weighted hydrostatic stress, used in
#                mechano-chemical coupling of theta_m.
###############################################################
[AuxVariables]
  [./i_mag]
    family = MONOMIAL
    order  = FIRST
  [../]

  [./sigma_h_xi]
    family = MONOMIAL
    order  = FIRST
  [../]
[]

###############################################################
# Dampers
#
#   BoundingValueNodalDamper keeps xi and theta_m within a
#   reasonable range during Newton iterations to improve
#   robustness.
###############################################################
[Dampers]
  [./limit_xi]
    type        = BoundingValueNodalDamper
    variable    = xi
    max_value   = 1.2
    min_value   = -0.2
    min_damping = 1e-8
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
#
#   xi:
#     - In the electrode: SmoothCircleIC defines a circular void
#       (xi ≈ 0 inside the void, xi ≈ 1 in surrounding Li).
#     - In the electrolyte: xi = 0 (void-like / non-metal),
#       and no xi kernels act there (xi remains fixed).
#
#   theta_m:
#     - Uniform initial value my_theta_m0 in both electrode and
#       electrolyte.
#
#   phi:
#     - 0 V everywhere at t = 0.
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

  # xi in the electrolyte: fixed to 0.0 (void-like / non-metal)
  # Note: no xi kernels act in 'electrolyte', so xi remains 0 there.
  [./ic_xi_electrolyte_void]
    type     = ConstantIC
    variable = xi
    value    = 0.0
    block    = 'electrolyte'
  [../]

  # theta_m: uniform in both electrode and electrolyte
  [./ic_theta_uniform]
    type     = ConstantIC
    variable = theta_m
    value    = ${my_theta_m0}
    block    = 'electrode electrolyte'
  [../]

  # phi: initial 0 V in the whole domain
  [./ic_phi_zero]
    type     = ConstantIC
    variable = phi
    value    = 0.0
  [../]
[]

###############################################################
# Boundary Conditions
#
#  (A) Electrostatic BCs for phi:
#
#    - Left boundary  (electrode / current collector):
#          phi = 0 (Dirichlet, grounded electrode).
#
#    - Right boundary (electrolyte / counter electrode):
#          n · (sigma ∇phi) = -i_app
#          (i_app > 0 => positive current leaving the domain).
#
#    - Top & bottom:
#          n · (sigma ∇phi) = 0  (electrically insulated).
#
#    NOTE / TODO:
#      When my_i_app ≠ 0, the conductivity jump and strong
#      coupling often cause convergence problems. This is part
#      of the known numerical challenge for this test.
#
#  (B) Interfacial BC for theta_m (Faradaic coupling):
#
#    - On the internal interface 'electrode_electrolyte',
#      ThetaFaradaicBC relates Li flux to normal current.
#
#  (C) Mechanical BCs:
#
#    - Top / bottom: no vertical displacement (u·n = 0, i.e. disp_y = 0).
#    - Right boundary: no horizontal displacement (u·n = 0, disp_x = 0).
#    - Left boundary: prescribed stack pressure (traction) in x-direction.
###############################################################
[BCs]
  # (A1) Left boundary: Dirichlet potential, grounded electrode
  [./phi_left]
    type     = DirichletBC
    variable = phi
    boundary = 'left'
    value    = 0.0        # [V]
  [../]

  # (A2) Right boundary: prescribed outward current density i_app
  #   Using Ohm's law i = -sigma ∇phi:
  #     i_n = i_app  =>  sigma ∇phi·n = -i_app.
  #   For MatDiffusion, NeumannBC "value" = sigma ∇phi·n.
  #   TODO: non-zero my_i_app + strong conductivity contrast
  #         + four-field coupling => convergence issues.
  [./phi_right_current]
    type     = NeumannBC
    variable = phi
    boundary = 'right'
    value    = -${my_i_app}     # = -i_app [A/m^2]
  [../]

  # (A3) Top and bottom: insulated (no normal current)
  [./phi_top_bottom_insulated]
    type     = NeumannBC
    variable = phi
    boundary = 'top bottom'
    value    = 0.0
  [../]

  # (B) theta_m on internal electrode/electrolyte interface:
  #     Faradaic coupling between Li flux and normal current density.
  [./theta_faradaic_interface]
    type         = ThetaFaradaicBC
    variable     = theta_m
    boundary     = 'electrode_electrolyte'

    # Coupled electric potential
    phi          = phi

    # Material property name for conductivity (sigma)
    conductivity = conductivity

    # Faraday constant [C/mol] and molar volumes
    F            = 96485.0
    Omega_L      = 13.1e-6   # [m^3/mol] reference molar volume of Li
    normal_sign  = 1.0       # +1: normal from electrode to electrolyte
  [../]

  # (C1) Mechanical BC: top/bottom – fixed in y-direction
  [./y_fixed_tb]
    type     = DirichletBC
    variable = disp_y
    boundary = 'top bottom'
    value    = 0.0
  [../]

  # (C2) Mechanical BC: right external boundary – fixed in x-direction
  [./x_fixed_right]
    type     = DirichletBC
    variable = disp_x
    boundary = 'right'
    value    = 0.0
  [../]

  # (C3) Mechanical BC: left external boundary – applied stack pressure
  #   n·σ·n = -p_applied along x-direction (compressive).
  [./left_stack_pressure]
    type      = Pressure
    variable  = disp_x
    boundary  = 'left'
    factor    = 6.0e5      # [Pa] = 0.6 MPa
    component = 0          # x-direction
  [../]
[]

###############################################################
# Solid mechanics physics
#
#   - Small-strain, quasi-static.
#   - Eigenstrain from chemical swelling (theta_m, xi).
#   - Anand viscoplasticity model + elastic response in electrode.
###############################################################
[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain            = SMALL
    incremental       = true
    eigenstrain_names = 'chem_strain'

    add_variables     = true
    generate_output   = 'strain_xx strain_yy strain_xy
                         stress_xx stress_yy stress_xy
                         vonmises_stress'
  [../]
[]

###############################################################
# Kernels
#
#  (A) xi: Allen–Cahn phase-field in 'electrode'
#
#      ∂xi/∂t = -L [ g'(xi) - kappa ∇²xi
#                    + (R T / Ω_L) h'(xi)
#                      ln((1 - theta_m)/(1 - theta_m0))
#                    + (Ω_v / Ω_L) h'(xi) ψ_e ]
#
#  (B) theta_m: diffusion + coupling + mechano-chemical term
#
#      h(xi) ∂theta_m/∂t + theta_m h'(xi) ∂xi/∂t
#        = ∇·( D(theta_m, xi) ∇theta_m )
#          + mechano-chemical contribution (LiMechCouplingKernel).
#
#  (C) phi: conduction
#
#      ∇·( conductivity ∇phi ) = 0.
###############################################################
[Kernels]
  # ---- (A1) xi time derivative --------------------------------
  [./xi_time]
    type     = TimeDerivative
    variable = xi
    block    = 'electrode'
  [../]

  # ---- (A2) xi chemical driving term from Li occupancy --------
  [./xi_chem_log]
    type     = LiVoidAllenCahn
    variable = xi
    mob_name = L                # mobility L [m^2/(N·s)]

    theta_m  = theta_m
    h_name   = h_xi             # interpolation function h(xi)

    R        = 8.314            # [J/(mol·K)]
    T        = 298.0            # [K]
    Omega_L  = 13.1e-6          # [m^3/mol]
    theta_m0 = ${my_theta_m0}   # reference occupancy
    block    = 'electrode'
  [../]

  # ---- (A3) xi double-well term: -L g'(xi) --------------------
  [./xi_bulk_dw]
    type     = AllenCahn
    variable = xi
    f_name   = g_xi             # material property g(xi) and g'(xi)
    mob_name = L
    block    = 'electrode'
  [../]

  # ---- (A4) xi elastic energy contribution --------------------
  #   Term proportional to h'(xi) * elastic_energy.
  [./xi_elastic_coupling]
    type           = LiVoidElasticAllenCahn
    variable       = xi
    block          = 'electrode'
    elastic_energy = elastic_energy
    h_name         = h_xi
    Omega_v        = 6.0e-6
    Omega_L        = 13.1e-6
  [../]

  # ---- (A5) xi gradient energy term: L kappa ∇²xi -------------
  [./xi_grad]
    type       = ACInterface
    variable   = xi
    kappa_name = kappa
    mob_name   = L
    block      = 'electrode'
  [../]

  # ---------- (B) theta_m: mass + coupling + diffusion ----------

  # (B1) Mass term: h(xi) ∂theta_m/∂t  (both blocks)
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

  # (B4) Mechano-chemical contribution using grad(sigma_h_xi)
  [./theta_mech_coupling]
    type        = LiMechCouplingKernel     # mech-coupled contribution with grad(sigma_h_xi)
    variable    = theta_m
    Deff        = D_eff
    h_xi        = h_xi
    sigma_h_xi  = sigma_h_xi

    Omega_Li    = 13.1e-6
    Omega_v     = 6.0e-6
    R           = 8.314
    T           = 298
    block       = 'electrode'
  [../]

  # ---------- (C) phi: conduction with conductivity ----------
  [./phi_diffusion]
    type        = MatDiffusion
    variable    = phi
    diffusivity = conductivity   # material property
  [../]
[]

###############################################################
# Auxiliary Kernels
#
#   - current_magnitude : i_mag = | -sigma ∇phi |.
#   - aux_sigma_h_xi    : xi-weighted hydrostatic stress for
#                         mechano-chemical coupling.
###############################################################
[AuxKernels]
  [./current_magnitude]
    type         = CurrentMagnitudeAux
    variable     = i_mag
    phi          = phi
    conductivity = conductivity
  [../]

  [./aux_sigma_h_xi]
    type           = HydrostaticStressXiAux
    variable       = sigma_h_xi
    elastic_strain = elastic_strain
    h_xi           = h_xi
    bulk_modulus   = 4.0833e3   # [Pa]; check magnitude vs. intended GPa scale
    block          = 'electrode'
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
#  Phase-field:
#    - g(xi)   = w xi^2 (1 - xi)^2
#    - h(xi)   = xi^2 (xi^2 - 3 xi + 3) + 1e-8  (regularized)
#
#  Diffusivity:
#    - D(theta_m, xi) = D_eff * h(xi) / max(1 - theta_m, 1e-3)
#
#  Conductivity:
#    - Electrode: sigma_s^xi(xi) (high in solid Li, vanishing in void).
#    - Electrolyte: constant sigma_el.
#
#  Mechanics:
#    - Phase-field-degraded isotropic elasticity tensor.
#    - Anand viscoplasticity (electrode).
#    - Chemical eigenstrain from Li swelling.
#    - Elastic energy density psi_e used to couple back to xi.
###############################################################
[Materials]
  # Constant phase-field and diffusion parameters
  [./pf_constants]
    type        = GenericConstantMaterial
    prop_names  = 'D_eff     L        kappa    w'
    prop_values = '7.5e-13  1e-9     4.5e-7  3.5e6'
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
    expression        = 'xi*xi*(xi*xi - 3*xi + 3) + 1.0e-8'
    property_name     = h_xi
    block             = 'electrode electrolyte'
    outputs           = out_exodus
  [../]

  # Effective diffusivity: conc_diffusivity(theta_m, xi)
  [./pf_conc_diffusivity]
    type                    = DerivativeParsedMaterial  # ParsedMaterial also possible
    coupled_variables       = 'theta_m'
    material_property_names = 'D_eff h_xi'
    expression              = 'D_eff * h_xi / max(1 - theta_m, 1e-3)'
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
    expression        = 'theta_m * xi * (4*xi^2 - 9*xi + 6)'
    derivative_order  = 2
    block             = 'electrode'    # can be extended to 'electrolyte' if needed
    outputs           = out_exodus
  [../]

  # Electrode conductivity: sigma_s^xi(xi)
  [./sigma_electrode]
    type                 = ParsedMaterial
    block                = 'electrode'
    property_name        = conductivity
    coupled_variables    = 'xi'
    constant_names       = 'sigma_s0'
    constant_expressions = '1.1e7'  # [S/m]
    expression           = 'sigma_s0 * max(0.0, min(1.0, pow(xi,15) * (pow(xi,4) - 3*xi*xi + 3)))'
    outputs              = out_exodus
  [../]

  # Electrolyte conductivity: constant sigma_el
  [./sigma_electrolyte]
    type        = GenericConstantMaterial
    block       = 'electrolyte'
    prop_names  = 'conductivity'
    prop_values = '5.5e-6'     # [S/m], adjust per reference
    outputs     = out_exodus
  [../]

  # Electrode elasticity: phase-field-degraded isotropic tensor
  [./elasticity_tensor]
    type          = ComputePFIsotropicElasticityTensor
    bulk_modulus  = 4.0833e3      # [Pa]; verify against intended elastic stiffness
    shear_modulus = 1.8846e3
    h_xi          = h_xi
    block         = 'electrode'
  [../]

  # Anand viscoplasticity model (electrode only)
  [./anand_vp]
    type          = AnandViscoplasticityStressUpdate
    block         = 'electrode'

    A            = 0.0       # To activate creep, set to reference value (e.g. paper Table 2)
    Q            = 3.7e4
    gas_constant = 8.314
    m            = 0.15

    S0           = 2.0e6
    Sa0          = 1.1e6

    H0           = 1.0e7
    a_exp        = 2.0
    n            = 0.05

    T_constant   = 298.0
    yield_ratio  = 0.1
    absolute_tolerance = 1e-8
  [../]

  # Combine elastic + inelastic (Anand) response
  [./radial_return_stress]
    type                     = ComputeMultipleInelasticStress
    block                    = 'electrode'
    tangent_operator         = elastic
    inelastic_models         = 'anand_vp'
    perform_finite_strain_rotations = false
  [../]

  # Chemical eigenstrain due to Li concentration (electrode only)
  [./chem_eigen]
    type             = ChemicalEigenstrain
    block            = 'electrode'
    eigenstrain_name = chem_strain

    xi       = xi
    theta_m  = theta_m

    theta_m0 = 0.999999998
    Omega_L  = 13.1e-6
    Omega_v  = 6.0e-6
  [../]

  # Elastic strain energy density psi_e (electrode only)
  [./pf_elastic_energy]
    type              = ComputePFElasticEnergy
    elastic_strain    = elastic_strain
    elasticity_tensor = elasticity_tensor
    elastic_energy    = elastic_energy
    block             = 'electrode'
  [../]
[]

###############################################################
# Postprocessors
#
#   - pp_num_dofs : total number of DOFs.
#   - pp_dt       : current time step size.
#   - pp_run_time : total wall-clock run time.
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
# Preconditioning
#
#   SMP (Split Matrix Preconditioning) with full coupling.
#   NOTE: For better robustness, further tuning or more advanced
#         preconditioners may be needed (TODO).
###############################################################
[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

###############################################################
# Executioner
#
#   - Transient BDF2 scheme.
#   - Fully coupled solution via PJFNK.
#   - Iteration-based adaptive time stepping.
#
#   NOTE / TODO:
#     For non-zero my_i_app and strong property jumps, expect
#     sensitivity to solver parameters; additional tuning and
#     mesh adaptivity will likely be required for robust runs.
###############################################################
[Executioner]
  type       = Transient
  scheme     = bdf2
  solve_type = PJFNK

  # Linear solver (PETSc / Hypre BoomerAMG)
  petsc_options_iname  = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value  = 'hypre    boomeramg      31                  0.7'

  line_search = bt

  l_tol      = 1e-4
  l_max_its  = 30
  nl_max_its = 25
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-18

  end_time = 240.0
  dtmax    = 50.0
  dtmin    = 1e-8

  [./TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1e-8       # [s] initial time step
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
