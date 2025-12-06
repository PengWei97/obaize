//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "Material.h"

/**
 * Compute elastic strain energy density using the (possibly degraded)
 * elasticity tensor:
 *
 *   psi_e^xi = 0.5 * epsilon_e : C^xi : epsilon_e
 *
 * 模板参数 is_ad 控制是否使用 AD 类型。
 */

template <bool is_ad>
class ComputePFElasticEnergyTempl : public Material
{
public:
  static InputParameters validParams();

  ComputePFElasticEnergyTempl(const InputParameters & parameters);

  using Material::_qp;

  virtual void computeQpProperties() override;

protected:
  /// Elastic strain tensor epsilon_e
  const GenericMaterialProperty<RankTwoTensor, is_ad> & _elastic_strain;

  /// (Possibly degraded) elasticity tensor C^xi
  const GenericMaterialProperty<RankFourTensor, is_ad> & _elasticity_tensor;

  /// Elastic energy density psi_e^xi
  GenericMaterialProperty<Real, is_ad> & _elastic_energy;
};

// 非 AD 和 AD typedef，方便在输入文件中使用
typedef ComputePFElasticEnergyTempl<false> ComputePFElasticEnergy;
typedef ComputePFElasticEnergyTempl<true>  ADComputePFElasticEnergy;
