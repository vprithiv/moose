//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetAdvectionSUPG.h"

registerMooseObject("LevelSetApp", LevelSetAdvectionSUPG);

InputParameters
LevelSetAdvectionSUPG::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription(
      "SUPG stablization term for the advection portion of the level set equation.");
  params.addRequiredCoupledVar("velocity", "Velocity vector variable.");
  return params;
}

LevelSetAdvectionSUPG::LevelSetAdvectionSUPG(const InputParameters & parameters)
  : ADKernelGrad(parameters), _velocity(adCoupledVectorValue("velocity"))
{
}

ADRealVectorValue
LevelSetAdvectionSUPG::precomputeQpResidual()
{
  computeQpVelocity();
  ADReal tau;
  if (_velocity.norm() > 1.0e-10)
    tau = _current_elem->hmin() / (2 * _velocity.norm());
  else
    tau = 0.0;
  return (tau * _velocity) * (_velocity * _grad_u[_qp]);
}
