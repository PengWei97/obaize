# ------------------------------------------------------------
# Step 1: 1D Li diffusion (occupancy θ) – no phase field, no mechanics
# Governing equation:
#   dθ/dt = ∇ · [ D_eff / (1 - θ) ∇θ ]
# Domain: 0 <= x <= 10 μm  (use SI units: 1e-5 m)
# ------------------------------------------------------------

my_filename = s2_diffusion_Li
my_interval = 10

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 5
  xmax = 10.0e-6
  ymax = 1.0e-6
  elem_type = QUAD4
[]

[Variables]
  [./theta_m]
    family = LAGRANGE
    order  = FIRST
  [../]
[]

# ------------------------
# Initial Conditions
#   Step profile:
#   θ = 0.9 for x < L/2, θ = 0.1 for x >= L/2
# ------------------------
[ICs]
  [./theta_m_ic]
    type = BoundingBoxIC
    variable = theta_m
    x1 = 0.0
    y1 = 0.0
    x2 = 5.0e-6
    y2 = 1.0e-6

    # int_width = 2.0e
    inside = 0.9
    outside = 0.1
  [../]
[]

# ------------------------
# Boundary Conditions
#   No-flux (natural Neumann) on both ends ⇒ nothing to specify
#   Optionally, you can enforce zero gradient explicitly:
# ------------------------
[BCs]
  [./concBCs]
    type = NeumannBC
    variable = theta_m
    boundary = 'left right top bottom'
    value = 0.0
  [../]
[]

# ------------------------
# Kernels
#   TimeDerivative: (w, dθ/dt)
#   Diffusion with nonlinear D = D_eff/(1 - θ)
# ------------------------
[Kernels]
  [./time]
    type = TimeDerivative
    variable = theta_m
  [../]
  [./LiDiffusionPDE]
    type = MatDiffusion
    variable = theta_m
    diffusivity = conc_diffusivity # D = D_eff/(1 - θ)
  [../]
[]

# ------------------------
# Materials
#   Constant effective diffusion coefficient D_eff
# ------------------------
[Materials]
  [./Material_prop_const]
    type = GenericConstantMaterial
    prop_names  = 'D_eff'
    prop_values = '7.5e-13'   # m^2/s (example value) '7.5e-13
  [../]
  [./conc_diffusivity]
    type = ParsedMaterial
    block = 0
    coupled_variables = 'theta_m'
    material_property_names = 'D_eff'

    expression = 'D_eff/(1-theta_m)'
    property_name = conc_diffusivity
    outputs = my_exodus
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'
  solve_type = 'PJFNK'

  l_max_its = 30
  l_tol = 1.0e-4
  nl_max_its = 25
  nl_rel_tol = 1.0e-9

  end_time = 200.0
  dtmax = 1.0
  dt = 1e-2
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-3       
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  # [./Adaptivity]
  #   initial_adaptivity = 0
  #   refine_fraction    = 0.4
  #   coarsen_fraction   = 0.05
  #   max_h_level        = 1
  # [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    time_step_interval = ${my_interval}
    type = Checkpoint
    additional_execute_on = 'INITIAL FINAL'
  [../]  
  [my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    type = Nemesis
    time_step_interval = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  [pgraph]
    type = PerfGraphOutput
    time_step_interval = ${my_interval}
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 5         # Default is 0
    execute_on = 'TIMESTEP_END FINAL'
  []
  print_linear_residuals = false
[]