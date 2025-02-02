E = 6.88e4
nu = 0.25

[GlobalParams]
  large_kinematics = true
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
  []
[]

[Kernels]
  [sdx]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  []
  [sdy]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  []
  [sdz]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
  [pull_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = 't'
    preset = false
  []
[]

[Materials]
  [elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [compute_strain]
    type = ComputeLagrangianStrain
    displacements = 'disp_x disp_y disp_z'
  []
  [flow_stress]
    type = DerivativeParsedMaterial
    f_name = flow_stress
    function = '320+688*effective_plastic_strain'
    material_property_names = 'effective_plastic_strain'
    additional_derivative_symbols = 'effective_plastic_strain'
    derivative_order = 2
    compute = false
  []
  [compute_stress]
    type = ComputeSimoHughesJ2PlasticityStress
    flow_stress_material = flow_stress
  []
[]

[Postprocessors]
  [sxx]
    type = ElementAverageValue
    variable = sxx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [exx]
    type = ElementAverageValue
    variable = exx
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[AuxVariables]
  [sxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 0
      index_j = 0
    []
  []
  [exx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 0
      index_j = 0
    []
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 5e-4
  end_time = 1e-1
[]

[Outputs]
  csv = true
[]
