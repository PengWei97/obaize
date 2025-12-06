//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "ComputePFElasticEnergy.h"

registerMooseObject("obaizeApp", ComputePFElasticEnergy);
registerMooseObject("obaizeApp", ADComputePFElasticEnergy);

template <bool is_ad>
InputParameters
ComputePFElasticEnergyTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute elastic strain energy density using the degraded elasticity tensor: "
      "psi_e^xi = 0.5 * epsilon_e : C^xi : epsilon_e.");

  params.addRequiredParam<MaterialPropertyName>(
      "elastic_strain",
      "Name of the elastic strain material property (RankTwoTensor).");

  params.addRequiredParam<MaterialPropertyName>(
      "elasticity_tensor",
      "Name of the (possibly phase-field-degraded) elasticity tensor (RankFourTensor).");

  params.addParam<MaterialPropertyName>(
      "elastic_energy",
      "elastic_energy",
      "Name of the output elastic energy density material property.");

  return params;
}

template <bool is_ad>
ComputePFElasticEnergyTempl<is_ad>::ComputePFElasticEnergyTempl(
    const InputParameters & parameters)
  : Material(parameters),
    _elastic_strain(
        this->template getGenericMaterialProperty<RankTwoTensor, is_ad>(
            this->template getParam<MaterialPropertyName>("elastic_strain"))),
    _elasticity_tensor(
        this->template getGenericMaterialProperty<RankFourTensor, is_ad>(
            this->template getParam<MaterialPropertyName>("elasticity_tensor"))),
    _elastic_energy(
        this->template declareGenericProperty<Real, is_ad>(
            this->template getParam<MaterialPropertyName>("elastic_energy")))
{
}

template <bool is_ad>
void
ComputePFElasticEnergyTempl<is_ad>::computeQpProperties()
{
  const GenericRankTwoTensor<is_ad> & eps = _elastic_strain[_qp];
  const GenericRankFourTensor<is_ad> & C  = _elasticity_tensor[_qp];

  // sigma^xi = C^xi : epsilon_e
  const GenericRankTwoTensor<is_ad> sigma = C * eps;

  // psi_e^xi = 0.5 * sigma^xi : epsilon_e
  _elastic_energy[_qp] = 0.5 * sigma.doubleContraction(eps);
}

// 显式实例化模板
template class ComputePFElasticEnergyTempl<false>;
template class ComputePFElasticEnergyTempl<true>;
