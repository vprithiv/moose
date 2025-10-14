//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMolarAvgAdvection.h"

registerMooseObject("MooseApp", ADMolarAvgAdvection);

InputParameters
ADMolarAvgAdvection::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription("Advection kernel using a molar-average velocity material property.");
  params.addRequiredParam<MaterialPropertyName>(
      "molar_avg_velocity", "Name of the molar-average velocity material property.");
  return params;
}

ADMolarAvgAdvection::ADMolarAvgAdvection(const InputParameters & parameters)
  : ADKernelValue(parameters),
    _velocity(getADMaterialProperty<RealVectorValue>(getParam<MaterialPropertyName>("molar_avg_velocity")))
{
}

ADReal
ADMolarAvgAdvection::precomputeQpResidual()
{
  return _velocity[_qp] * _grad_u[_qp];
}

