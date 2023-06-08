//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralVariablePostprocessor.h"

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class LevelSetAverageValue : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  LevelSetAverageValue(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  /// The accumulated volume to return as a PostprocessorValue
  Real _volume;

  /// The level set contour to consider for computing inside vs. outside of the volume
  const Real & _threshold;

  /// Flag for triggering the internal volume calculation
  const bool _inside;

  /// level set variable used as coupled variable for thresholding volume
  const VariableValue & _levelset_var;

};
