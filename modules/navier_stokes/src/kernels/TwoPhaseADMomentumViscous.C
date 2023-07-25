//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseADMomentumViscous.h"
#include "MooseMesh.h"
#include "INSADObjectTracker.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", TwoPhaseADMomentumViscous);

InputParameters
TwoPhaseADMomentumViscous::validParams()
{
  InputParameters params = ADVectorKernelGrad::validParams();
  params.addRequiredCoupledVar("level_set_gradient",
                               "Regularized gradient of the level set variable");
  params.addParam<MaterialPropertyName>("mu_name", "mu", "Dynamic viscosity variable");
  params.addRequiredParam<Real>("surface_tension", "surface tension coefficient.");
  return params;
}

TwoPhaseADMomentumViscous::TwoPhaseADMomentumViscous(
    const InputParameters & parameters)
  : ADVectorKernelGrad(parameters),
    _mu(getADMaterialProperty<Real>("mu_name")),
    _grad_c(adCoupledVectorValue("level_set_gradient")),
    _sigma(getParam<Real>("surface_tension"))
{
}


ADRealTensorValue
TwoPhaseADMomentumViscous::precomputeQpResidual()
{
  ADRealVectorValue tolerance_val = RealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE);
  ADReal grad_c_norm = (_grad_c[_qp] + tolerance_val).norm();
  ADRealVectorValue n_phi = _grad_c[_qp] /grad_c_norm;
  ADReal delta = grad_c_norm;

  RankTwoTensor iden(RankTwoTensor::initIdentity);
  ADRankTwoTensor proj;
  proj.vectorOuterProduct(n_phi, n_phi);
  proj = iden - proj;

  ADRankTwoTensor tau = _mu[_qp] * (_grad_u[_qp] + _grad_u[_qp].transpose());
  ADRankTwoTensor tau_phi =_sigma*delta*proj;
//   return tau+tau_phi;
  return tau_phi;

}
