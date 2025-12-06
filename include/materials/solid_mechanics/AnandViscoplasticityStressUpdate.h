//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RadialReturnStressUpdate.h"

/**
 * J2 viscoplasticity model based on Anand's flow rule.
 * Templated on is_ad to support both Real and ADReal.
 */

template <bool is_ad>
class AnandViscoplasticityStressUpdateTempl : public RadialReturnStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  AnandViscoplasticityStressUpdateTempl(const InputParameters & parameters);

  using Material::_qp;
  using RadialReturnStressUpdateTempl<is_ad>::_base_name;
  using RadialReturnStressUpdateTempl<is_ad>::_three_shear_modulus;
  using RadialReturnStressUpdateTempl<is_ad>::_dt;
  using RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn;

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;

  virtual GenericReal<is_ad>
  computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                  const GenericReal<is_ad> & scalar) override;

  virtual GenericReal<is_ad>
  computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                    const GenericReal<is_ad> & scalar) override;

  virtual void iterationFinalize(const GenericReal<is_ad> & scalar) override;

  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plasticStrainIncrement) override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  /// Anand hardening law: Sa^{n+1} = Sa_old + dt * Sa_dot(F_cr_old, Sa_old, T)
  GenericReal<is_ad>
  computeUpdatedFlowResistance(const GenericReal<is_ad> & F_cr_old,
                               const GenericReal<is_ad> & T,
                               const GenericReal<is_ad> & Sa_old) const;

  /// Anand flow rate F_cr = A exp(-Q/RT) [sinh(sigma_eff / Sa)]^(1/m)
  GenericReal<is_ad>
  computeFlowRate(const GenericReal<is_ad> & sigma_eff,
                  const GenericReal<is_ad> & Sa,
                  const GenericReal<is_ad> & T) const;

  // ---------- Anand material parameters (Real, not AD) ----------
  const Real _A;
  const Real _Q;
  const Real _R;
  const Real _m;
  const Real _H0;
  const Real _a_exp;
  const Real _n;
  const Real _S0;
  const Real _Sa0;

  /// activation condition factor: viscoplastic flow active if sigma_eff > yield_ratio * Sa_old
  const Real _yield_ratio;

  // ---------- Temperature handling ----------
  const GenericVariableValue<is_ad> * const _temperature;
  const bool _use_temperature;
  const Real _T_constant;

  // ---------- Stateful material properties ----------
  /// Anand flow resistance Sa
  GenericMaterialProperty<Real, is_ad> & _Sa;
  const MaterialProperty<Real> & _Sa_old;

  /// Viscoplastic strain tensor
  GenericMaterialProperty<RankTwoTensor, is_ad> & _vp_strain;
  const MaterialProperty<RankTwoTensor> & _vp_strain_old;

  // ---------- Internal scalar state for the scalar solve ----------
  GenericReal<is_ad> _activation_condition;
  GenericReal<is_ad> _F_cr;
  GenericReal<is_ad> _dFcr_dscalar;
};

// Non-AD and AD typedefs
typedef AnandViscoplasticityStressUpdateTempl<false> AnandViscoplasticityStressUpdate;
typedef AnandViscoplasticityStressUpdateTempl<true> ADAnandViscoplasticityStressUpdate;
