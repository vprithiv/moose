//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "TagAuxBase.h"

/**
 * The value of a tagged vector for a given node and a given variable is coupled to
 * the current AuxVariable. TagVectorAux returns the coupled nodal value.
 */
template <typename T>
class TagVectorAuxTempl : public TagAuxBase<AuxKernelTempl<T>>
{
public:
  static InputParameters validParams();

  TagVectorAuxTempl(const InputParameters & parameters);

protected:
  virtual T computeValue() override;

  const typename OutputTools<T>::VariableValue & _v;
  const MooseVariableBase & _v_var;

  using TagAuxBase<AuxKernelTempl<T>>::_qp;
  using TagAuxBase<AuxKernelTempl<T>>::_var;
  using TagAuxBase<AuxKernelTempl<T>>::_scaled;
};

typedef TagVectorAuxTempl<Real> TagVectorAux;
typedef TagVectorAuxTempl<RealVectorValue> VectorTagVectorAux;