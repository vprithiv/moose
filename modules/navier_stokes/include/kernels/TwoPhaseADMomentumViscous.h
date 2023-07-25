//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include "ADKernelGrad.h"

class TwoPhaseADMomentumViscous : public ADVectorKernelGrad
{
public:
  static InputParameters validParams();  
  TwoPhaseADMomentumViscous(const InputParameters & parameters);

protected:
  virtual ADRealTensorValue precomputeQpResidual() override;
  const ADMaterialProperty<Real> & _mu;
  const ADVectorVariableValue & _grad_c;
  const Real _sigma;
};
