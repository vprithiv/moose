//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetOlssonCircularArc.h"

registerMooseObject("LevelSetApp", LevelSetOlssonCircularArc);

InputParameters
LevelSetOlssonCircularArc::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Implementation of 'CircularArc' ranging from 0 to 1.");
  params.addParam<RealVectorValue>(
      "center", RealVectorValue(0.5, 0.5, 0), "The center of the CircularArc.");
  params.addParam<Real>("radius", 0.15, "The radius of the CircularArc.");
  params.addParam<Real>("epsilon", 0.01, "The interface thickness.");
  MooseEnum concaveflag("false=0 true=1", "true");
  params.addParam<MooseEnum>("concave_facing_right", concaveflag, "The location of the area/volume to be computed.");
  return params;
}

LevelSetOlssonCircularArc::LevelSetOlssonCircularArc(const InputParameters & parameters)
  : Function(parameters),
    _center(getParam<RealVectorValue>("center")),
    _radius(getParam<Real>("radius")),
    _epsilon(getParam<Real>("epsilon")),
    _concave_facing_right(getParam<MooseEnum>("concave_facing_right") == "true")
{
}

Real
LevelSetOlssonCircularArc::value(Real /*t*/, const Point & p) const
{
  if (_concave_facing_right)
  {
    if (p(0)<=_center(0))
    {
    const Real x = ((p - _center).norm() - _radius) / _epsilon;
    return 1.0 / (1 + std::exp(x));
    }
    else
    {
    return 1.0;
    }
  }
  else
  {
    if (p(0)>_center(0))
    {
    const Real x = ((p - _center).norm() - _radius) / _epsilon;
    return 1.0 / (1 + std::exp(x));
    }
    else
    {
    return 1.0;
    }

  }

}

ADReal
LevelSetOlssonCircularArc::value(const ADReal & /*t*/, const ADPoint & p) const
{
  if (_concave_facing_right)
  {
    if (p(0)<=_center(0))
    {
    const auto x = ((p - _center).norm() - _radius) / _epsilon;
    return 1.0 / (1 + std::exp(x));
    }
    else
    {
    return 1.0;
    }
  }
  else
  {
    if (p(0)>_center(0))
    {
    const auto x = ((p - _center).norm() - _radius) / _epsilon;
    return 1.0 / (1 + std::exp(x));
    }
    else
    {
    return 0.0;
    }

  }
}

RealGradient
LevelSetOlssonCircularArc::gradient(Real /*t*/, const Point & p) const
{
  Real norm = (p - _center).norm();
  Real g = (norm - _radius) / _epsilon;
  RealGradient output;

  Real g_prime;

  if (_concave_facing_right)
  {
    if (p(0)<=_center(0))
    {
        for (const auto i : make_range(Moose::dim))
        {
            g_prime = (p(i) - _center(i)) / (_epsilon * norm);
            output(i) = -(g_prime * std::exp(g)) / ((std::exp(g) + 1) * (std::exp(g) + 1));
        }
    }
    else
    {
    return 0.0;
    }
  }
  else
  {
    if (p(0)>_center(0))
    {
        for (const auto i : make_range(Moose::dim))
        {
            g_prime = (p(i) - _center(i)) / (_epsilon * norm);
            output(i) = -(g_prime * std::exp(g)) / ((std::exp(g) + 1) * (std::exp(g) + 1));
        }
    }
    else
    {
    return 0.0;
    }

  }

  return output;
}
 