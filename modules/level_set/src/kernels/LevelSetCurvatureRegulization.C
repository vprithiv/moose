//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetCurvatureRegulization.h"

registerMooseObject("LevelSetApp", LevelSetCurvatureRegulization);

InputParameters
LevelSetCurvatureRegulization::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("smoothes the curvature");
  params.addRequiredCoupledVar("grad_c", "The level set variable.");
  params.addRequiredParam<Real>("epsilon", "The epsilon coefficient.");
  return params;
}

LevelSetCurvatureRegulization::LevelSetCurvatureRegulization(
    const InputParameters & parameters)
  : ADKernel(parameters),
    _grad_c(adCoupledVectorValue("grad_c")),
    _epsilon(getParam<Real>("epsilon"))
{
}


ADReal
LevelSetCurvatureRegulization::computeQpResidual()
{
  ADReal s = std::numeric_limits<ADReal>::epsilon();
  ADRealVectorValue n(0, 0, 0);
  if (_grad_c[_qp].norm_sq() > 1.0e-8)
  {
    s = _grad_c[_qp].norm();
    n = _grad_c[_qp] / s;
  }
  return _test[_i][_qp] * _u[_qp] - _grad_test[_i][_qp] * (n - _epsilon * _grad_u[_qp]);
}
