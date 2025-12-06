//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AnandViscoplasticityStressUpdate.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

#include <cmath>

// 注册非 AD 和 AD 两个版本
registerMooseObject("obaizeApp", AnandViscoplasticityStressUpdate);
registerMooseObject("obaizeApp", ADAnandViscoplasticityStressUpdate);

template <bool is_ad>
InputParameters
AnandViscoplasticityStressUpdateTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("J2 viscoplasticity model based on Anand's flow rule. "
                             "The class solves for the effective viscoplastic strain increment "
                             "using an implicit integration scheme and updates the flow "
                             "resistance according to Anand's hardening law.");

  // ---------- Anand material parameters (same notation as in the paper) ----------

  params.addRequiredParam<Real>("A", "Pre-exponential factor A [1/s]");
  params.addRequiredParam<Real>("Q", "Activation energy Q [J/mol]");
  params.addParam<Real>("gas_constant", 8.314, "Gas constant R [J/(mol*K)]");
  params.addRequiredParam<Real>("m", "Strain-rate sensitivity exponent m (0 < m <= 1)");
  params.addRequiredParam<Real>("H0", "Hardening constant H0 [Pa]");
  params.addRequiredParam<Real>("a_exp", "Hardening sensitivity a [-]");
  params.addRequiredParam<Real>("n", "Deformation resistance sensitivity n [-]");
  params.addRequiredParam<Real>("S0", "Deformation resistance saturation coefficient S0 [Pa]");
  params.addRequiredParam<Real>("Sa0", "Initial value of the flow resistance Sa(t=0) [Pa]");

  // Activation condition for viscoplastic flow
  params.addParam<Real>("yield_ratio",
                        0.1,
                        "Viscoplastic flow is active if sigma_eff > yield_ratio * Sa_old");

  // Temperature: either a constant value or a coupled variable
  params.addParam<Real>("T_constant",
                        298.15,
                        "Constant temperature if no temperature variable is coupled [K]");
  params.addCoupledVar("temperature",
                       "Optional temperature variable; "
                       "if not provided, T_constant will be used");

  // Effective inelastic strain name (used by RadialReturnStressUpdate)
  params.set<std::string>("effective_inelastic_strain_name") = "effective_viscoplastic_strain";

  return params;
}

template <bool is_ad>
AnandViscoplasticityStressUpdateTempl<is_ad>::AnandViscoplasticityStressUpdateTempl(
    const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),

    // Anand parameters
    _A(this->template getParam<Real>("A")),
    _Q(this->template getParam<Real>("Q")),
    _R(this->template getParam<Real>("gas_constant")),
    _m(this->template getParam<Real>("m")),
    _H0(this->template getParam<Real>("H0")),
    _a_exp(this->template getParam<Real>("a_exp")),
    _n(this->template getParam<Real>("n")),
    _S0(this->template getParam<Real>("S0")),
    _Sa0(this->template getParam<Real>("Sa0")),
    _yield_ratio(this->template getParam<Real>("yield_ratio")),

    // Temperature handling
    _temperature(this->isParamValid("temperature")
                     ? &this->template coupledGenericValue<is_ad>("temperature")
                     : nullptr),
    _use_temperature(this->isParamValid("temperature")),
    _T_constant(this->template getParam<Real>("T_constant")),

    // Stateful material properties
    _Sa(this->template declareGenericProperty<Real, is_ad>("anand_flow_resistance")),
    _Sa_old(this->template getMaterialPropertyOld<Real>("anand_flow_resistance")),
    _vp_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>(
        this->_base_name + "viscoplastic_strain")),
    _vp_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>(
        this->_base_name + "viscoplastic_strain")),

    _activation_condition(-1.0),
    _F_cr(0.0),
    _dFcr_dscalar(0.0)
{
  // --- Stabilization: enforce a numerically safe range for m
  mooseAssert(_m > 0.05 && _m <= 1.0,
              "Anand parameter m is out of a numerically safe range (0.05, 1].");
}

template <bool is_ad>
void
AnandViscoplasticityStressUpdateTempl<is_ad>::initQpStatefulProperties()
{
  // Initial flow resistance Sa(t=0) and zero viscoplastic strain
  _Sa[_qp] = _Sa0;

  // --- Stabilization: floor Sa at a small fraction of S0
  const Real Sa_floor = 0.01 * _S0;
  if (_Sa[_qp] < Sa_floor)
    _Sa[_qp] = Sa_floor;

  _vp_strain[_qp].zero();
}

template <bool is_ad>
void
AnandViscoplasticityStressUpdateTempl<is_ad>::propagateQpStatefulProperties()
{
  // Copy old state to current step
  _Sa[_qp] = _Sa_old[_qp];
  _vp_strain[_qp] = _vp_strain_old[_qp];

  // --- Stabilization: keep Sa above a floor
  const Real Sa_floor = 0.01 * _S0;
  if (_Sa[_qp] < Sa_floor)
    _Sa[_qp] = Sa_floor;

  // Let the base class propagate its own state (effective inelastic strain, etc.)
  this->propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
AnandViscoplasticityStressUpdateTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  // Base class does trial elastic state and related bookkeeping
  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);

  // Use old Sa as initial guess for the current step
  _Sa[_qp] = _Sa_old[_qp];

  // --- Stabilization: floor Sa
  const Real Sa_floor = 0.01 * _S0;
  if (_Sa[_qp] < Sa_floor)
    _Sa[_qp] = Sa_floor;

  _vp_strain[_qp] = _vp_strain_old[_qp];

  // Simple activation condition: only attempt viscoplastic update
  // if effective_trial_stress is sufficiently larger than Sa_old
  _activation_condition = effective_trial_stress - _yield_ratio * _Sa_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
AnandViscoplasticityStressUpdateTempl<is_ad>::computeResidual(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> residual = 0.0;

  mooseAssert(_activation_condition != -1.0,
              "Activation condition was not updated by computeStressInitialize");

  // If the activation condition is not met, no viscoplastic flow is considered
  if (_activation_condition <= 0.0)
    return 0.0;

  // Effective stress in the scalar space:
  // sigma_eff = sigma_trial - 3G * Δp
  GenericReal<is_ad> sigma_eff = effective_trial_stress - _three_shear_modulus * scalar;

  // --- Stabilization: do not allow negative effective stress
  if (sigma_eff < 0.0)
    sigma_eff = 0.0;

  // Temperature at this point
  const GenericReal<is_ad> T =
      _use_temperature ? (*_temperature)[_qp] : GenericReal<is_ad>(_T_constant);

  // Use Sa from the previous time step inside the flow rule (explicit coupling)
  GenericReal<is_ad> Sa_used = _Sa_old[_qp];

  // --- Stabilization: Sa floor
  const Real Sa_floor = 0.01 * _S0;
  if (Sa_used < Sa_floor)
    Sa_used = Sa_floor;

  // Compute Anand flow rate F_cr and store it for later use
  _F_cr = computeFlowRate(sigma_eff, Sa_used, T);

  // Derivative dF_cr/d(Δp).
  //
  // F_cr = A exp(-Q/RT) [sinh(x)]^(1/m),   x = sigma_eff / Sa
  // dF/dsigma = F_cr * (1/m) * coth(x) / Sa
  // dsigma/d(Δp) = -3G
  // => dF/d(Δp) = -3G * F_cr * (1/m) * coth(x) / Sa
  GenericReal<is_ad> x = sigma_eff / Sa_used;

  // --- Stabilization: limit x to avoid overflow in sinh/cosh
  const Real x_max = 20.0;
  if (x > x_max)
    x = x_max;
  else if (x < -x_max)
    x = -x_max;

  const GenericReal<is_ad> sinh_x = std::sinh(x);
  const GenericReal<is_ad> cosh_x = std::cosh(x);

  GenericReal<is_ad> coth_x = 0.0;
  if (std::abs(sinh_x) > 1e-16)
    coth_x = cosh_x / sinh_x;
  else
    // safeguard: very small sinh -> use a large but finite value
    coth_x = (x >= 0.0 ? GenericReal<is_ad>(1e16) : GenericReal<is_ad>(-1e16));

  _dFcr_dscalar = -_three_shear_modulus * _F_cr * (1.0 / _m) * coth_x / Sa_used;

  // Residual: R = F_cr * dt - Δp
  residual = _F_cr * _dt - scalar;

  return residual;
}

template <bool is_ad>
GenericReal<is_ad>
AnandViscoplasticityStressUpdateTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & /*effective_trial_stress*/,
    const GenericReal<is_ad> & /*scalar*/)
{
  // If no viscoplastic flow is activated, the residual is effectively R = -Δp,
  // but the radial return framework will recognize that we are elastic.
  if (_activation_condition <= 0.0)
    return 1.0;

  // dR/d(Δp) = dF_cr/d(Δp) * dt - 1
  const GenericReal<is_ad> derivative = _dFcr_dscalar * _dt - 1.0;
  return derivative;
}

template <bool is_ad>
void
AnandViscoplasticityStressUpdateTempl<is_ad>::iterationFinalize(
    const GenericReal<is_ad> & /*scalar*/)
{
  // We do not update Sa inside the scalar Newton iteration.
  // Sa is updated explicitly in computeStressFinalize using the
  // last converged value of F_cr.
}

template <bool is_ad>
GenericReal<is_ad>
AnandViscoplasticityStressUpdateTempl<is_ad>::computeUpdatedFlowResistance(
    const GenericReal<is_ad> & F_cr_old,
    const GenericReal<is_ad> & T,
    const GenericReal<is_ad> & Sa_old) const
{
  // Helper: A exp(-Q / (R T))
  const GenericReal<is_ad> Ae = _A * std::exp(-_Q / (_R * T));

  // Sa_star = S0 * (F_cr / Ae)^n
  GenericReal<is_ad> Sa_star = _S0;
  if (F_cr_old > 0.0 && Ae > 0.0)
    Sa_star = _S0 * std::pow(F_cr_old / Ae, _n);

  // --- Stabilization: Sa_old floor
  const Real Sa_floor = 0.01 * _S0;
  GenericReal<is_ad> Sa_old_safe = Sa_old;
  if (Sa_old_safe < Sa_floor)
    Sa_old_safe = Sa_floor;

  const GenericReal<is_ad> ratio = 1.0 - Sa_old_safe / Sa_star;
  const GenericReal<is_ad> abs_ratio = std::abs(ratio);
  const GenericReal<is_ad> sign_ratio = (ratio >= 0.0 ? 1.0 : -1.0);

  // Sa_dot = H0 * |1 - Sa/Sa_star|^a * sign(1 - Sa/Sa_star) * F_cr
  const GenericReal<is_ad> Sa_dot =
      _H0 * std::pow(abs_ratio, _a_exp) * sign_ratio * F_cr_old;

  GenericReal<is_ad> Sa_new = Sa_old_safe + _dt * Sa_dot;

  // --- Stabilization: limit single-step relative change of Sa
  const Real max_rel_change = 0.5; // at most ±50% per time step
  const GenericReal<is_ad> Sa_max = (1.0 + max_rel_change) * Sa_old_safe;
  const GenericReal<is_ad> Sa_min_rel = (1.0 - max_rel_change) * Sa_old_safe;

  if (Sa_new > Sa_max)
    Sa_new = Sa_max;
  else if (Sa_new < Sa_min_rel)
    Sa_new = Sa_min_rel;

  // --- Stabilization: enforce absolute floor
  if (Sa_new < Sa_floor)
    Sa_new = Sa_floor;

  return Sa_new;
}

template <bool is_ad>
GenericReal<is_ad>
AnandViscoplasticityStressUpdateTempl<is_ad>::computeFlowRate(
    const GenericReal<is_ad> & sigma_eff,
    const GenericReal<is_ad> & Sa,
    const GenericReal<is_ad> & T) const
{
  const GenericReal<is_ad> Ae = _A * std::exp(-_Q / (_R * T));

  // --- Stabilization: floor Sa based on S0
  const Real Sa_floor = 0.01 * _S0;
  GenericReal<is_ad> Sa_safe = Sa;
  if (Sa_safe < Sa_floor)
    Sa_safe = Sa_floor;

  // x = sigma_eff / Sa_safe
  GenericReal<is_ad> x = sigma_eff / Sa_safe;

  // --- Stabilization: limit x to avoid overflow in sinh(x)
  const Real x_max = 20.0;
  if (x > x_max)
    x = x_max;
  else if (x < 0.0)
    x = 0.0; // effective stress should not drive negative viscoplastic flow

  GenericReal<is_ad> sinh_x = std::sinh(x);

  // --- Stabilization: avoid log(0) later
  if (sinh_x < 1e-16)
    sinh_x = 1e-16;

  // F_cr = Ae * [sinh(x)]^(1/m)
  // Use log form to avoid overflow: log(F_cr) = log(Ae) + (1/m)*log(sinh(x))
  GenericReal<is_ad> F_cr = 0.0;
  if (Ae > 0.0)
  {
    GenericReal<is_ad> log_F = std::log(Ae) + (1.0 / _m) * std::log(sinh_x);

    // --- Stabilization: cap log(F_cr) to keep exp(log_F) within double range
    const Real log_F_max = 700.0; // exp(700) ~ 1e304
    if (log_F > log_F_max)
      log_F = log_F_max;

    F_cr = std::exp(log_F);
  }
  else
    F_cr = 0.0;

  return F_cr;
}

template <bool is_ad>
void
AnandViscoplasticityStressUpdateTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plasticStrainIncrement)
{
  // Update the viscoplastic strain tensor using the increment provided
  // by the radial return mapping
  _vp_strain[_qp] += plasticStrainIncrement;

  // Update the flow resistance Sa using the converged F_cr
  const GenericReal<is_ad> T =
      _use_temperature ? (*_temperature)[_qp] : GenericReal<is_ad>(_T_constant);
  _Sa[_qp] = computeUpdatedFlowResistance(_F_cr, T, _Sa_old[_qp]);
}

// 显式实例化模板
template class AnandViscoplasticityStressUpdateTempl<false>;
template class AnandViscoplasticityStressUpdateTempl<true>;
