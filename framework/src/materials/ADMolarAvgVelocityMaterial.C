//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMolarAvgVelocityMaterial.h"

registerMooseObject("MooseApp", ADMolarAvgVelocityMaterial);

InputParameters
ADMolarAvgVelocityMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes the molar-average velocity from species fluxes, molecular "
                             "weights, a mixture density, and a coupled bulk velocity.");

  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "fluxes", "Material property names for the species fluxes except for the eliminated species.");
  params.addRequiredParam<std::vector<Real>>("molecular_weights", "Molecular weights for all species.");
  params.addRequiredParam<MaterialPropertyName>("density", "Name of the mixture density material property.");
  params.addRequiredParam<MaterialPropertyName>(
      "molar_avg_velocity", "Name of the molar-average velocity material property to declare.");
  params.addRequiredCoupledVar("velocity", "Bulk velocity vector variable.");

  return params;
}

ADMolarAvgVelocityMaterial::ADMolarAvgVelocityMaterial(const InputParameters & parameters)
  : Material(parameters),
    _molecular_weights(getParam<std::vector<Real>>("molecular_weights")),
    _density(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("density"))),
    _velocity(adCoupledVectorValue("velocity")),
    _molar_avg_velocity(
        declareADProperty<RealVectorValue>(getParam<MaterialPropertyName>("molar_avg_velocity")))
{
  const auto & flux_names = getParam<std::vector<MaterialPropertyName>>("fluxes");

  if (_molecular_weights.size() < 2)
    paramError("molecular_weights", "At least two molecular weights must be provided.");

  if (flux_names.size() + 1 != _molecular_weights.size())
    paramError("fluxes", "Number of fluxes must be exactly one less than the number of molecular weights.");

  _fluxes.reserve(flux_names.size());
  for (const auto & name : flux_names)
    _fluxes.push_back(&getADMaterialProperty<RealVectorValue>(name));
}

void
ADMolarAvgVelocityMaterial::computeQpProperties()
{
  ADRealVectorValue weighted_flux_sum;
  weighted_flux_sum.zero();

  const Real last_molecular_weight = _molecular_weights.back();

  for (std::size_t i = 0; i < _fluxes.size(); ++i)
  {
    const ADRealVectorValue & flux_i = (*_fluxes[i])[_qp];
    weighted_flux_sum += (_molecular_weights[i] - last_molecular_weight) * flux_i;
  }

  const ADReal inv_density = 1.0 / _density[_qp];
  const ADRealVectorValue correction = inv_density * weighted_flux_sum;

  _molar_avg_velocity[_qp] = _velocity[_qp] + correction;
}

