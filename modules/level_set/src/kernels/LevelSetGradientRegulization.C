//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetGradientRegulization.h"

registerMooseObject("LevelSetApp", LevelSetGradientRegulization);

InputParameters
LevelSetGradientRegulization::validParams()
{
  InputParameters params = ADVectorKernel::validParams();
  params.addRequiredCoupledVar("c", "Level set variable");
  params.addParam<bool>("normal_regulization", false, "Level set variable");
  return params;
}

LevelSetGradientRegulization::LevelSetGradientRegulization(
    const InputParameters & parameters)
  : ADVectorKernel(parameters),
    _grad_c(adCoupledGradient("c")),
    _normal_regulization(getParam<bool>("normal_regulization"))
{
}


ADReal
LevelSetGradientRegulization::computeQpResidual()
{
  if (!_normal_regulization)
    return _test[_i][_qp] * (_u[_qp] - _grad_c[_qp]);
  else
  {
    if (_grad_c[_qp].norm_sq() > 1.0e-8)
      return _test[_i][_qp] * (_u[_qp] - _grad_c[_qp] / _grad_c[_qp].norm());
    else
      return _test[_i][_qp] * _u[_qp];
  }
}

