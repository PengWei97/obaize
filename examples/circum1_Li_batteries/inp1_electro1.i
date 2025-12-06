###############################################
# Pure electric-field test problem
# - Only solve for electric potential phi
# - Phase-field xi is a fixed background field
#   used to define a phase-dependent conductivity
#   in the Li electrode (metal vs void).
#
# Domain:
#   Left half  (0–50 μm)  : Li metal electrode + void
#   Right half (50–100 μm): solid electrolyte
#
# Goal:
#   Debug electric equation, boundary conditions
#   and phase-field-dependent conductivity,
#   without coupling back to theta_m or xi evolution.
#
# Author: Wei Peng (simplified E-field version)
# Date  : 2025-12-06
###############################################

# File naming and basic parameters
my_filename = s1_phi_only_test
my_interval = 20

# Applied macroscopic current density [A/m^2]
my_i_app = 1.0
my_theta_m0 = 0.999999998
###############################################################
# Mesh
#   - 2D square: 100 μm × 100 μm
#   - Resolution: 250 × 250 QUAD4 elements
#   - Block "electrode"   : 0–50 μm  (Li metal)
#   - Block "electrolyte" : 50–100 μm (solid electrolyte)
###############################################################
[Mesh]
  [./base]
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

  # Subdomain 1: Li electrode, left half (0–50 μm)
  [./electrode_block]
    type        = SubdomainBoundingBoxGenerator
    input       = base
    block_id    = 1
    block_name  = 'electrode'
    bottom_left = '0.0 0.0 0.0'
    top_right   = '50e-6 100e-6 0.0'
  [../]

  # Subdomain 2: solid electrolyte, right half (50–100 μm)
  [./electrolyte_block]
    type        = SubdomainBoundingBoxGenerator
    input       = electrode_block
    block_id    = 2
    block_name  = 'electrolyte'
    bottom_left = '50e-6 0.0 0.0'
    top_right   = '100e-6 100e-6 0.0'
  [../]

  # Internal interface sideset between electrode and electrolyte
  #   normal points from electrode -> electrolyte
  [./electrode_electrolyte_interface]
    type          = SideSetsBetweenSubdomainsGenerator
    input         = electrolyte_block
    primary_block = 'electrode'
    paired_block  = 'electrolyte'
    new_boundary  = 'electrode_electrolyte'
    replace       = false
  [../]

  type = MeshGeneratorMesh
[]


###############################################################
# Variables
#
#  phi : electric potential [V], defined over both
#        electrode and electrolyte.
###############################################################
[Variables]
  [./phi]
    family = LAGRANGE
    order  = FIRST
  [../]

  [./theta_m]
    family = LAGRANGE
    order  = FIRST
    block  = 'electrode'
  [../]
[]

###############################################################
# Auxiliary variables
#
#  xi    : fixed phase-field describing void geometry.
#          xi ≈ 1 in solid Li, xi ≈ 0 in voids and electrolyte.
#
#  i_mag : magnitude of current density |i| [A/m^2]
#          computed from i = -sigma * grad(phi).
###############################################################
[AuxVariables]
  [./xi]
    family = LAGRANGE
    order  = FIRST
  [../]

  [./i_mag]
    family = MONOMIAL
    order  = FIRST
  [../]
[]

###############################################################
# Initial Conditions
#
#  xi:
#    - In the electrode: diffuse circular void
#      (xi ≈ 0 inside, xi ≈ 1 outside).
#    - In the electrolyte: xi = 0 (void-like / insulating).
#
#  phi:
#    - Initial guess: 0 V everywhere.
###############################################################
[ICs]
  # xi: diffuse circular void inside the electrode
  [./ic_xi_electrode_void]
    type      = SmoothCircleIC
    variable  = xi
    x1        = 50e-6      # [m] void center x (near interface)
    y1        = 50e-6      # [m] void center y
    radius    = 15e-6      # [m] void radius
    int_width = 1e-6       # [m] interface half width
    invalue   = 0.0        # void / insulating
    outvalue  = 1.0        # solid Li
    block     = 'electrode'
  [../]

  # xi in the electrolyte: set to 0.0 (no conduction from Li)
  [./ic_xi_electrolyte]
    type      = ConstantIC
    variable  = xi
    value     = 0.0
    block     = 'electrolyte'
  [../]

  # phi: initial 0 V in the whole domain
  [./ic_phi_zero]
    type     = ConstantIC
    variable = phi
    value    = 0.0
  [../]

  # theta_m: initially uniform in both electrode and electrolyte
  [./ic_theta_uniform]
    type      = ConstantIC
    variable  = theta_m
    value     = ${my_theta_m0}
    block     = 'electrode'
  [../]
[]


###############################################################
# Boundary Conditions for phi
#
#   Left boundary  (electrode / current collector): phi = 0.
#   Right boundary (electrolyte / counter electrode):
#       n · (sigma ∇phi) = -i_app   (i_app > 0 means
#       positive current leaving the domain).
#   Top & bottom: electrically insulated,
#       n · (sigma ∇phi) = 0.
###############################################################
[BCs]
  # 1) Left boundary: Dirichlet potential, grounded electrode
  [./phi_left]
    type     = DirichletBC
    variable = phi
    boundary = 'left'
    value    = 0.0        # [V]
  [../]

  # 2) Right boundary: prescribed outward current density i_app
  #    Using Ohm's law i = -sigma grad(phi):
  #      i_n = i_app  =>  sigma grad(phi)·n = -i_app.
  #    For MatDiffusion, NeumannBC "value" = sigma grad(phi)·n.
  [./phi_right_current]
    type     = NeumannBC
    variable = phi
    boundary = 'right'
    value    = -${my_i_app}     # = -i_app [A/m^2]
  [../]

  # 3) Top and bottom: insulated (no normal current)
  [./phi_top_bottom_insulated]
    type     = NeumannBC
    variable = phi
    boundary = 'top bottom'
    value    = 0.0
  [../]
[]


###############################################################
# Kernels
#
# Electric potential equation over the whole domain:
#
#   ∇·( σ(x, xi) ∇phi ) = 0
#
# implemented via MatDiffusion with a phase- and
# block-dependent conductivity.
###############################################################
[Kernels]
  [./phi_diffusion]
    type        = MatDiffusion
    variable    = phi
    diffusivity = conductivity   # material property
  [../]

  # (B1) Mass term: h(xi) ∂theta_m/∂t
  [./theta_mass]
    type           = MatTimeDerivative
    variable       = theta_m
    mat_prop_value = const_test
    block     = 'electrode'
  [../] 
[]


###############################################################
# Auxiliary Kernels
#
#  Compute magnitude of current density:
#    i = -sigma * grad(phi),
#    i_mag = |i|.
###############################################################
[AuxKernels]
  [./current_magnitude]
    type         = CurrentMagnitudeAux
    variable     = i_mag
    phi          = phi
    conductivity = conductivity   # same material property
  [../]
[]


###############################################################
# Materials
#
#  Conductivity in the electrode (metal + void):
#
#    sigma_s^xi = sigma_min + sigma_s0 * f(xi)
#
#  where
#    f(xi) = xi^15 (xi^4 - 3 xi^2 + 3)
#
#  so that:
#    f(0) = 0,  f(1) = 1,  f'(0) = 0,
#  and sigma_min provides a small lower bound to avoid
#  a fully degenerate PDE in the void.
#
#  Conductivity in the electrolyte:
#    sigma_el = const.
###############################################################
[Materials]
  # Constant phase-field and diffusion parameters
  [./pf_constants]
    type        = GenericConstantMaterial
    prop_names  = 'D_eff     L        kappa    w const_test'
    prop_values = '7.5e-13  1e-6     4.5e-7  3.5e6 1.0'
  [../]
  # Electrode conductivity: sigma_s^xi(xi)
  [./sigma_electrode]
    type                  = ParsedMaterial
    block                 = 'electrode'
    property_name         = conductivity
    coupled_variables     = 'xi'
    constant_names        = 'sigma_s0'
    constant_expressions  = '1.1e7'  # [S/m]
    expression            = 'sigma_s0 * max(0.0, min(1.0, pow(xi,15) * (pow(xi,4) - 3*xi*xi + 3)))'
    outputs                 = out_exodus
  [../]

  # Electrolyte conductivity: constant sigma_el
  [./sigma_electrolyte]
    type        = GenericConstantMaterial
    block       = 'electrolyte'
    prop_names  = 'conductivity'
    prop_values = '5.5e-6'     # [S/m], adjust per paper
    outputs                 = out_exodus
  [../]
  # Effective diffusivity: conc_diffusivity(theta_m, xi)
  [./pf_conc_diffusivity]
    type                    = DerivativeParsedMaterial
    coupled_variables       = 'theta_m'
    material_property_names = 'D_eff'
    expression              = 'D_eff/max(1-theta_m, 1e-3)'
    property_name           = conc_diffusivity
    block     = 'electrode'
    derivative_order        = 2
    outputs                 = out_exodus
  [../]
[]


###############################################################
# Postprocessors
#   - pp_num_dofs : total number of DOFs
#   - pp_run_time : total wall-clock runtime
###############################################################
[Postprocessors]
  [./pp_num_dofs]
    type = NumDOFs
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
  nl_abs_tol = 1e-12

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
  [./out_exodus]
    type                  = Exodus
    file_base             = ./ex_${my_filename}/out_${my_filename}
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
    execute_on         = 'FINAL'
  [../]
[]
