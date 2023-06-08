//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetAverageValue.h"
#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", LevelSetAverageValue);

InputParameters
LevelSetAverageValue::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addClassDescription("Computes the volumetric average of a variable inside or outside of levelset contour");
  params.addParam<Real>(
      "threshold", 0.0, "The level set threshold to consider for computing area/volume.");
  params.addRequiredCoupledVar("levelset_var",
                               "levelset variable to threshold volume");

  MooseEnum loc("inside=0 outside=1", "inside");
  params.addParam<MooseEnum>("location", loc, "The location of the area/volume to be computed.");  
  return params;
}

LevelSetAverageValue::LevelSetAverageValue(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters),
    _threshold(getParam<Real>("threshold")),
    _inside(getParam<MooseEnum>("location") == "inside"),
    _levelset_var(coupledValue("levelset_var"))
{
}

void
LevelSetAverageValue::initialize()
{
  ElementIntegralVariablePostprocessor::initialize();
  _volume = 0;
}

void
LevelSetAverageValue::execute()
{
  Real sum = 0;
  Real cnt = 0;
  Real n = _u.size();

  // Perform the check for inside/outside outside the qp loop for speed
  if (_inside)
  {
    for (_qp = 0; _qp < n; ++_qp)
      if (_levelset_var[_qp] <= _threshold)
        {
          sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral();
          cnt++;
        }

  }
  else
  {
    for (_qp = 0; _qp < n; ++_qp)
      if (_levelset_var[_qp] > _threshold)
        {
          sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral();
          cnt++;
        }
  }
  _volume += cnt / n * _current_elem_volume;
  _integral_value+= sum;

}

Real
LevelSetAverageValue::getValue()
{
  return _integral_value / _volume;
}



void
LevelSetAverageValue::finalize()
{
  ElementIntegralVariablePostprocessor::finalize();
  gatherSum(_volume);
}

void
LevelSetAverageValue::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const LevelSetAverageValue & pps = static_cast<const LevelSetAverageValue &>(y);
  _volume += pps._volume;
}
