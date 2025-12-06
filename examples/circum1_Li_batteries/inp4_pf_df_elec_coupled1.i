###############################################
# Phase-field evolution of a circular void
# in a Li metal / solid electrolyte system
#
# Equations:
#   - xi      : Allen–Cahn phase-field equation in the Li electrode,
#               including a double-well bulk term, a gradient-energy
#               term, and a chemical contribution from Li occupancy
#               (theta_m).
#
#   - theta_m : Li site occupancy, evolving via a simplified diffusion
#               equation in the electrode and a weak coupling to xi.
#
#   - phi     : Electric potential, satisfying a steady conduction
#               equation with conductivity depending on xi and the
#               material (electrode / electrolyte).
#
# Goal:
#   Test the combined interface-driven (xi) and chemical-potential-
#   driven (theta_m) evolution of a void inside the Li electrode,
#   coupled to an electrostatic potential field phi.
#
# NOTE on current boundary conditions:
#   In this input file, we apply NO external current:
#     - All phi boundary conditions are homogeneous Neumann:
#           n · (sigma ∇phi) = 0
#       which corresponds to zero normal current density
#           i_n = -sigma (∇phi · n) = 0
#       on every external boundary (left/right/top/bottom).
#   即：这里所有与电流相关的边界条件均为 0，对应无外加电流。
#
# Author: Wei Peng
# Date  : 2025-12-06
###############################################

# Output naming and basic parameters
my_filename  = s4_c1_pf_df_elec_with_void1
my_interval  = 20
my_theta_m0  = 0.999999998
# my_i_app   = 1.0   # reserved for future non-zero current tests

###############################################################
# Mesh
#
#   - 2D square: 100 μm × 100 μm
#   - Resolution: 250 × 250 QUAD4 elements
#
#   Subdomains:
#     - Block "electrode"   : 0 – 50 μm  (Li metal)
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
# Variables
#
#   xi      : Phase-field order parameter
#             (xi = 1: solid Li, xi = 0: void/electrolyte).
#
#   theta_m : Li site occupancy (dimensionless), providing a
#             chemical driving force for xi and following a
#             simplified diffusion equation.
#
#   phi     : Electric potential, solving a conduction equation
#             with conductivity that depends on xi and the
#             material (electrode or electrolyte).
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
#   i_mag : Magnitude of current density |i| [A/m^2],
#           where i = -sigma * grad(phi).
###############################################################
[AuxVariables]
  [./i_mag]
    family = MONOMIAL
    order  = FIRST
  [../]
[]

###############################################################
# Dampers
#
#   BoundingValueNodalDamper keeps xi and theta_m within a
#   reasonable range during Newton iterations, improving
#   numerical stability.
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
#     - In the electrode: SmoothCircleIC defines a circular void,
#         xi ≈ 0 inside the void,
#         xi ≈ 1 in the surrounding Li metal.
#     - In the electrolyte: xi = 0 (void-like / non-metal region).
#       No xi kernels act in 'electrolyte', so xi remains fixed
#       there during the simulation.
#
#   theta_m:
#     - Initially uniform theta_m0 in both electrode and
#       electrolyte. Only in the electrode does theta_m diffuse
#       (via theta_diffusion); in the electrolyte it remains
#       essentially uniform.
#
#   phi:
#     - Initially 0 V throughout the domain. With zero current
#       boundary conditions and no sources, phi stays near 0.
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
  [./ic_xi_electrolyte_void]
    type     = ConstantIC
    variable = xi
    value    = 0.0
    block    = 'electrolyte'
  [../]

  # theta_m: initially uniform in both electrode and electrolyte
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
# Boundary Conditions for phi (electric potential)
#
#   In this test case, we impose homogeneous Neumann BC on
#   ALL external boundaries:
#
#       n · (sigma ∇phi) = 0   =>  i_n = 0
#
#   This corresponds to zero normal current density on each
#   boundary, i.e. NO externally applied current through
#   any boundary.
###############################################################
[BCs]
  # phi: zero normal current on all external boundaries
  [./phi_zero_current_all_boundaries]
    type     = NeumannBC
    variable = phi
    boundary = 'left right top bottom'
    value    = 0.0
  [../]

  # Note:
  #   - theta_m uses natural (zero-flux) Neumann BC by default,
  #     corresponding to no Li flux crossing the boundaries.
  #   - For future electrochemical coupling, Faradaic-type BCs
  #     can be added at the internal electrode/electrolyte
  #     interface, but are not used in this zero-current test.
[]
# ---------------------------------------------------------------------------
# (Reference BCs for non-zero applied current are kept below as comments
#  for future use.)
#
# [BCs]
#   # 1) Left boundary: Dirichlet potential, grounded electrode
#   [./phi_left]
#     type     = DirichletBC
#     variable = phi
#     boundary = 'left'
#     value    = 0.0        # [V]
#   [../]
#
#   # 2) Right boundary: prescribed outward current density i_app
#   #    Using Ohm's law i = -sigma grad(phi):
#   #      i_n = i_app  =>  sigma grad(phi)·n = -i_app.
#   [./phi_right_current]
#     type     = NeumannBC
#     variable = phi
#     boundary = 'right'
#     value    = -${my_i_app}     # = -i_app [A/m^2]
#   [../]
#
#   # 3) Top and bottom: insulated (no normal current)
#   [./phi_top_bottom_insulated]
#     type     = NeumannBC
#     variable = phi
#     boundary = 'top bottom'
#     value    = 0.0
#   [../]
#
#   # 4) Example Faradaic BC for theta_m at the internal interface
#   # [./theta_faradaic_interface]
#   #   type         = ThetaFaradaicBC
#   #   variable     = theta_m
#   #   boundary     = 'electrode_electrolyte'
#   #   phi          = phi
#   #   conductivity = conductivity
#   #   F            = 96485.0
#   #   Omega_L      = 13.1e-6
#   #   normal_sign  = 1.0
#   # [../]
# []
# ---------------------------------------------------------------------------

###############################################################
# Kernels
#
# (A) xi equation (Allen–Cahn) on 'electrode':
#
#   ∂xi/∂t = -L [ g'(xi) - kappa ∇²xi
#                 + (R T / Ω_L) h'(xi)
#                   ln((1 - theta_m)/(1 - theta_m0)) ]
#
#   Implemented via:
#     - xi_time     : TimeDerivative (∂xi/∂t)
#     - xi_chem_log : chemical driving term from theta_m
#     - xi_bulk_dw  : -L g'(xi)          (double-well)
#     - xi_grad     : +L kappa ∇²xi      (gradient energy)
#
# (B) theta_m equation (mass + coupling + diffusion):
#
#   h(xi) ∂theta_m/∂t + theta_m h'(xi) ∂xi/∂t
#     = ∇·( D(theta_m, xi) ∇theta_m )
#
#   Implemented via:
#     - theta_mass       : h(xi) ∂theta_m/∂t
#     - theta_xi_coupled : theta_m h'(xi) ∂xi/∂t
#     - theta_diffusion  : ∇·(conc_diffusivity ∇theta_m)
#
# (C) phi equation (steady conduction):
#
#   ∇·( conductivity ∇phi ) = 0
#
#   Implemented via:
#     - phi_diffusion : MatDiffusion with conductivity
#                       defined by materials.
###############################################################
[Kernels]
  # ---- (A1) xi time derivative --------------------------------
  [./xi_time]
    type     = TimeDerivative
    variable = xi
    block    = 'electrode'
  [../]

  # ---- (A2) xi chemical driving term from Li occupancy --------
  # Contribution:
  #   -(L R T / Ω_L) h'(xi) ln[(1 - theta_m)/(1 - theta_m0)]
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

  # ---------- (B) theta_m: mass + coupling + diffusion ----------

  # (B1) Mass term: h(xi) ∂theta_m/∂t (in both blocks)
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

  # ---------- (C) phi equation: conduction with conductivity ----------
  [./phi_diffusion]
    type        = MatDiffusion
    variable    = phi
    diffusivity = conductivity
  [../]
[]

###############################################################
# Auxiliary Kernels
#
#   i_mag = |i| = |-sigma ∇phi|
###############################################################
[AuxKernels]
  [./current_magnitude]
    type         = CurrentMagnitudeAux
    variable     = i_mag
    phi          = phi
    conductivity = conductivity
  [../]
[]

###############################################################
# Materials
#
#  Parameters:
#    - D_eff : reference diffusivity [m^2/s]
#    - L     : Allen–Cahn mobility [m^2/(N·s)]
#    - kappa : gradient-energy coefficient [N]
#    - w     : double-well height [Pa]
#
#  Double-well potential:
#    g(xi) = w xi^2 (1 - xi)^2
#
#  Interpolation function (regularized):
#    h(xi) = xi^2 (xi^2 - 3 xi + 3) + 1e-8
#
#  Effective diffusivity:
#    D(theta_m, xi) = D_eff * h(xi) / max(1 - theta_m, 1e-3)
#
#  Conductivity:
#    - In electrode: sigma_s^xi(xi), strongly xi-dependent
#    - In electrolyte: constant sigma_el
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
    expression        = 'xi*xi*(xi*xi - 3*xi + 3) + 1.0e-8'
    property_name     = h_xi
    block             = 'electrode electrolyte'
    outputs           = out_exodus
  [../]

  # Effective diffusivity: conc_diffusivity(theta_m, xi)
  [./pf_conc_diffusivity]
    type                    = DerivativeParsedMaterial
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

  # Electrode conductivity: sigma_s^xi(xi) (high in solid Li, ~0 in void)
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
    prop_values = '5.5e-6'     # [S/m], adjust as needed
    outputs     = out_exodus
  [../]
[]

###############################################################
# Postprocessors
#
#   - pp_num_dofs : total number of DOFs
#   - pp_dt       : current time step size
#   - pp_run_time : total wall-clock run time
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
#
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
    dt                 = 1e-1       # [s] initial time step
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
