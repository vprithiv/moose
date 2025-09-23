//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernelValue.h"

/**
 * Kernel applying advection using a molar-average velocity material property.
 */
class ADMolarAvgAdvection : public ADKernelValue
{
public:
  static InputParameters validParams();

  ADMolarAvgAdvection(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  /// Molar-average velocity material property
  const ADMaterialProperty<RealVectorValue> & _velocity;
};

