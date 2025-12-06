//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HydrostaticStressXiAux.h"

#include "MooseError.h"
#include <algorithm>          // std::max
#include "metaphysicl/raw_type.h"  // MetaPhysicL::raw_value

// 显式实例化模板，避免多个 TU 里重复编译
template class HydrostaticStressXiAuxTempl<false>;
template class HydrostaticStressXiAuxTempl<true>;

// 注册两个 Moose 对象
registerMooseObject("obaizeApp", HydrostaticStressXiAux);
registerMooseObject("obaizeApp", ADHydrostaticStressXiAux);

//
// validParams
//
template <bool is_ad>
InputParameters
HydrostaticStressXiAuxTempl<is_ad>::validParams()
{
  InputParameters params = AuxKernelTempl<Real>::validParams();

  params.addClassDescription(
      "Compute phase-field-degraded hydrostatic stress "
      "sigma_h^xi = h_eff(xi) * K * tr(epsilon_e), "
      "with h_eff = kappa0 + (1-kappa0)*max(h(xi), 0).");

  params.addRequiredParam<MaterialPropertyName>(
      "elastic_strain",
      "Name of the elastic strain material property (RankTwoTensor).");

  params.addRequiredParam<MaterialPropertyName>(
      "h_xi",
      "Name of the interpolation function material property h(xi).");

  params.addRequiredParam<Real>(
      "bulk_modulus",
      "Reference bulk modulus K of the solid phase.");

  params.addParam<Real>(
      "kappa0",
      1.0e-4,
      "Residual stiffness fraction in the void: "
      "h_eff = kappa0 + (1-kappa0)*h(xi).");

  return params;
}

//
// 构造函数
//
template <bool is_ad>
HydrostaticStressXiAuxTempl<is_ad>::HydrostaticStressXiAuxTempl(
    const InputParameters & parameters)
  : AuxKernelTempl<Real>(parameters),
    _elastic_strain(
        this->template getGenericMaterialProperty<RankTwoTensor, is_ad>(
            parameters.get<MaterialPropertyName>("elastic_strain"))),
    _h_xi(
        this->template getGenericMaterialProperty<Real, is_ad>(
            parameters.get<MaterialPropertyName>("h_xi"))),
    _bulk_modulus(this->template getParam<Real>("bulk_modulus")),
    _kappa0(this->template getParam<Real>("kappa0"))
{
  if (_bulk_modulus <= 0.0)
    mooseError("bulk_modulus must be positive in HydrostaticStressXiAux.");
}

//
// computeValue
//
template <bool is_ad>
Real
HydrostaticStressXiAuxTempl<is_ad>::computeValue()
{
  using std::max;
  using MetaPhysicL::raw_value;

  // 取出 h(xi)（可能是 Real 或 ADReal），统一抽成 Real
  const auto & h_val   = _h_xi[this->_qp];
  const Real   h_raw   = raw_value(h_val);
  const Real   h_eff   = _kappa0 + (1.0 - _kappa0) * max(h_raw, 0.0);
  const Real   K_xi    = h_eff * _bulk_modulus;

  // 取出 eps_e 并计算 trace，再抽 raw_value
  const auto & eps_e      = _elastic_strain[this->_qp];
  const auto   tr_eps_eAD = eps_e.trace();        // GenericReal<is_ad>
  const Real   tr_eps_e   = raw_value(tr_eps_eAD);

  return K_xi * tr_eps_e;
}
