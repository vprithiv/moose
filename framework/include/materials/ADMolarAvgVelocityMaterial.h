//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

#include "MooseTypes.h"

#include <vector>

/**
 * Material computing the molar-average velocity based on species fluxes and a bulk velocity.
 *
 * The material expects fluxes for the first N-1 species and eliminates the last species using the
 * constraint \sum_i J_i = 0. In an (N-1)-species basis the molar-average velocity is constructed as
 *
 *   \vec{u} = \vec{v} + \frac{1}{\rho} \sum_{i=1}^{N-1} (M_i - M_N) \vec{J}_i,
 *
 * where \vec{v} is the coupled bulk velocity, M_i the molecular weights, and \rho the mixture
 * density.
 */
class ADMolarAvgVelocityMaterial : public Material
{
public:
  static InputParameters validParams();

  ADMolarAvgVelocityMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Molecular weights for each species (size N)
  const std::vector<Real> _molecular_weights;
  /// Flux material properties for the first N-1 species
  std::vector<const ADMaterialProperty<RealVectorValue> *> _fluxes;
  /// Mixture density material property
  const ADMaterialProperty<Real> & _density;
  /// Coupled bulk velocity
  const ADVectorVariableValue & _velocity;
  /// Output molar-average velocity material property
  ADMaterialProperty<RealVectorValue> & _molar_avg_velocity;
};

