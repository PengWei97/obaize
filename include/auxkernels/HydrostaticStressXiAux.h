//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "AuxKernel.h"
#include "MooseTypes.h"
#include "ADReal.h"

/**
 * 模板版本的 HydrostaticStressXiAux：
 *
 *   sigma_h^xi = h_eff(xi) * K * tr(eps_e)
 *   h_eff = kappa0 + (1-kappa0) * max(h(xi), 0)
 *
 * - is_ad = false : 使用普通 MaterialProperty<T>
 * - is_ad = true  : 使用 ADMaterialProperty<T>，通过 raw_value() 抽取 Real 值
 */
template <bool is_ad>
class HydrostaticStressXiAuxTempl : public AuxKernelTempl<Real>
{
public:
  static InputParameters validParams();

  HydrostaticStressXiAuxTempl(const InputParameters & parameters);

protected:
  /// 计算 aux 变量值（总是返回 Real）
  virtual Real computeValue() override;

  /// 弹性应变 eps_e (RankTwoTensor<Real> 或 RankTwoTensor<ADReal>)
  const GenericMaterialProperty<RankTwoTensor, is_ad> & _elastic_strain;

  /// 插值函数 h(xi) (Real 或 ADReal)
  const GenericMaterialProperty<Real, is_ad> & _h_xi;

  /// 参考体积模量 K
  const Real _bulk_modulus;

  /// void 中残余刚度比例
  const Real _kappa0;
};

// 方便使用的别名：非 AD + AD 两个具体类型
using HydrostaticStressXiAux    = HydrostaticStressXiAuxTempl<false>;
using ADHydrostaticStressXiAux = HydrostaticStressXiAuxTempl<true>;
