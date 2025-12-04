# ------------------------------------------------------------
# Step 1: 1D Li diffusion (occupancy θ) – no phase field, no mechanics
# Governing equation:
#   dθ/dt = ∇ · [ D_eff / (1 - θ) ∇θ ]
# Domain: 0 <= x <= 10 μm  (use SI units: 1e-5 m)
#
# Note:
#   Adaptive mesh refinement (AMR) is intentionally disabled here.
#   AMR will be a future optimization target, but the current AMR
#   strategy is not yet mature and may introduce unnecessary
#   complexity at this stage.
# ------------------------------------------------------------

my_filename = s1_diffusion_Li
my_interval = 2

[Mesh]
  type      = GeneratedMesh
  dim       = 1            # true 1D diffusion
  nx        = 100          # number of elements along x
  xmax      = 1.0e-5       # 10 μm in SI units
  elem_type = EDGE2
[]

[Variables]
  [./theta_m]
    family = LAGRANGE
    order  = FIRST
  [../]
[]

# ------------------------------------------------------------
# Initial Conditions
# Step profile:
#   θ = 0.9 for x < L/2, θ = 0.1 for x >= L/2
# ------------------------------------------------------------
[ICs]
  [./theta_m_ic]
    type     = BoundingBoxIC
    variable = theta_m

    x1 = 0.0
    x2 = 5.0e-6        # half of the domain length

    inside  = 0.9
    outside = 0.1
  [../]
[]

# ------------------------------------------------------------
# Boundary Conditions
# Natural Neumann (no-flux) BC is implied by default for diffusion,
# so no explicit BC block is required.
# Uncomment the following block if you want to enforce zero gradient.
# ------------------------------------------------------------
# [BCs]
#   [./concBCs]
#     type     = NeumannBC
#     variable = theta_m
#     boundary = 'left right'
#     value    = 0.0
#   [../]
# []

# ------------------------------------------------------------
# Kernels
#   TimeDerivative: (w, dθ/dt)
#   MatDiffusion:   (∇w, D(θ) ∇θ) with D = D_eff / (1 - θ)
# ------------------------------------------------------------
[Kernels]
  [./time]
    type     = TimeDerivative
    variable = theta_m
  [../]

  [./LiDiffusionPDE]
    type        = MatDiffusion
    variable    = theta_m
    diffusivity = conc_diffusivity    # material property name
  [../]
[]

# ------------------------------------------------------------
# Materials
#   D_eff: constant effective diffusion coefficient
#   conc_diffusivity: D(θ) = D_eff / (1 - θ)
# ------------------------------------------------------------
[Materials]
  [./material_prop_const]
    type        = GenericConstantMaterial
    prop_names  = 'D_eff'
    prop_values = '7.5e-13'   # m^2·s^-1 (example value)
  [../]

  [./conc_diffusivity]
    type                       = ParsedMaterial
    coupled_variables          = 'theta_m'
    material_property_names    = 'D_eff'
    property_name              = conc_diffusivity
    expression                 = 'D_eff/(1 - theta_m + 1e-6)'
    outputs                    = my_exodus
  [../]
[]

# ------------------------------------------------------------
# Preconditioning / Linear solver
#   Use SMP (LU factorization) with a single preconditioned solve.
# ------------------------------------------------------------
[Preconditioning]
  [./pc]
    type = SMP
    full = true

    petsc_options_iname  = '-pc_type -ksp_type'
    petsc_options_value  = 'lu       preonly'
  [../]
[]

# ------------------------------------------------------------
# Executioner
#   Transient diffusion solved by implicit Euler.
#   Fixed mesh is used here; AMR may be added in future work.
# ------------------------------------------------------------
[Executioner]
  type       = Transient
  scheme     = 'implicit-euler'
  solve_type = 'NEWTON'

  l_max_its  = 60
  l_tol      = 1.0e-4
  nl_max_its = 25
  nl_rel_tol = 1.0e-9

  num_steps  = 5
  dt         = 1.0e-2
  dtmax      = 1.0
[]

# ------------------------------------------------------------
# Outputs
#   - Checkpoint files
#   - Exodus/Nemesis output
#   - CSV output
#   - PerfGraph for performance analysis
# ------------------------------------------------------------
[Outputs]
  [./my_checkpoint]
    file_base             = ./${my_filename}/out_${my_filename}
    time_step_interval    = ${my_interval}
    type                  = Checkpoint
    additional_execute_on = 'INITIAL FINAL'
  [../]

  [./my_exodus]
    file_base             = ./ex_${my_filename}/out_${my_filename}
    type                  = Nemesis
    time_step_interval    = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]

  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type      = CSV
  [../]

  [./pgraph]
    type               = PerfGraphOutput
    time_step_interval = ${my_interval}
    level              = 2
    heaviest_branch    = true
    heaviest_sections  = 5
    execute_on         = 'TIMESTEP_END FINAL'
  [../]

  print_linear_residuals = false
[]
