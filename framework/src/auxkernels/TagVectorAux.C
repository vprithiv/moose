//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TagVectorAux.h"

registerMooseObject("MooseApp", TagVectorAux);

template <typename T>
InputParameters
TagVectorAuxTempl<T>::validParams()
{
  InputParameters params = TagAuxBase<AuxKernel>::validParams();
  params.addRequiredParam<TagName>("vector_tag", "Tag Name this Aux works on");
  params.addClassDescription("Couple a tag vector, and return its nodal value");
  return params;
}

template <typename T>
TagVectorAuxTempl<T>::TagVectorAuxTempl(const InputParameters & parameters)
  : TagAuxBase<AuxKernelTempl<T>>(parameters),
    _v(this->template vectorTagValueHelper<T>("v", "vector_tag")),
    _v_var(*this->getFieldVar("v", 0))
{
  this->checkCoupledVariable(&_v_var, &_var);
}

template <typename T>
T
TagVectorAuxTempl<T>::computeValue()
{
  return _v[_qp] / (_scaled ? 1.0 : _v_var.scalingFactor());
}

template class TagVectorAuxTempl<Real>;
template class TagVectorAuxTempl<RealVectorValue>;
