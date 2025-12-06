//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "ComputePFIsotropicElasticityTensor.h"

#include <algorithm>
#include <cmath>

registerMooseObject("obaizeApp", ComputePFIsotropicElasticityTensor);

InputParameters
ComputePFIsotropicElasticityTensor::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();
  params.addClassDescription(
      "Isotropic elasticity tensor degraded by a phase-field interpolation h(xi).");

  params.addRequiredParam<Real>("bulk_modulus",
                                "Reference bulk modulus K of the material.");
  params.addRequiredParam<Real>("shear_modulus",
                                "Reference shear modulus G of the material.");
  params.addRequiredParam<MaterialPropertyName>(
      "h_xi", "Name of the interpolation function material property h(xi).");
  params.addParam<Real>(
      "kappa0",
      1.0e-4,
      "Residual stiffness fraction in the void: h_eff = kappa0 + (1-kappa0)*h(xi).");
  params.declareControllable("bulk_modulus shear_modulus");
  return params;
}

ComputePFIsotropicElasticityTensor::ComputePFIsotropicElasticityTensor(
    const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _h_xi(getMaterialProperty<Real>(parameters.get<MaterialPropertyName>("h_xi"))),
    _lambda0(0.0),
    _E0(0.0),
    _nu0(0.0),
    _effective_stiffness_local(0.0),
    _kappa0(getParam<Real>("kappa0"))
{
  if (_bulk_modulus <= 0.0)
    mooseError("Bulk modulus must be positive in material '" + name() + "'.");

  if (_shear_modulus <= 0.0)
    mooseError("Shear modulus must be positive in material '" + name() + "'.");

  issueGuarantee(_elasticity_tensor_name, Guarantee::ISOTROPIC);
  issueGuarantee("effective_stiffness", Guarantee::ISOTROPIC);
}

void
ComputePFIsotropicElasticityTensor::residualSetup()
{
  // Base Lamé constant from K and G
  _lambda0 = _bulk_modulus - 2.0 / 3.0 * _shear_modulus;

  // Base Young's modulus and Poisson's ratio
  const Real K = _bulk_modulus;
  const Real G = _shear_modulus;

  _E0 = (9.0 * K * G) / (3.0 * K + G);
  _nu0 = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G));

  if (_nu0 <= -1.0 || _nu0 >= 0.5)
    mooseError("Poisson's ratio inferred from K and G must lie in (-1, 0.5) "
               "in material '" +
               name() + "'.");
}

void
ComputePFIsotropicElasticityTensor::computeQpElasticityTensor()
{
  const Real h_raw = _h_xi[_qp];
  const Real h_clamped = std::max(h_raw, 0.0); // 避免负值
  const Real h_eff = _kappa0 + (1.0 - _kappa0) * h_clamped;

  // 不再有刚度完全为 0 的情况
  const Real G_eff     = h_eff * _shear_modulus;
  const Real K_eff     = h_eff * _bulk_modulus;
  const Real lambda_eff = h_eff * _lambda0;

  std::vector<Real> iso_const(2);
  iso_const[0] = lambda_eff;
  iso_const[1] = G_eff;

  _Cijkl.fillFromInputVector(iso_const, RankFourTensor::symmetric_isotropic);
  _elasticity_tensor[_qp] = _Cijkl;

  const Real E_eff  = h_eff * _E0;
  const Real nu_eff = _nu0;

  const Real bulk_like =
      std::sqrt((E_eff * (1.0 - nu_eff)) / ((1.0 + nu_eff) * (1.0 - 2.0 * nu_eff)));
  const Real shear_like = std::sqrt(G_eff);

  _effective_stiffness_local = std::max(bulk_like, shear_like);
  _effective_stiffness[_qp]  = _effective_stiffness_local;
}

