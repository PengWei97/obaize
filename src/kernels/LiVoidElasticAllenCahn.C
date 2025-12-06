#include "LiVoidElasticAllenCahn.h"

registerMooseObject("obaizeApp", LiVoidElasticAllenCahn);

InputParameters
LiVoidElasticAllenCahn::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();

  params.addClassDescription(
      "Allen–Cahn elastic driving term for Li–void phase-field model, "
      "using DerivativeParsedMaterial for h(xi) and a material psi_e.");

  // Elastic energy density psi_e
  params.addRequiredParam<MaterialPropertyName>(
      "elastic_energy", "Material property name for elastic energy density psi_e.");

  // Base name of interpolation function h(xi)
  params.addRequiredParam<MaterialPropertyName>(
      "h_name", "Base material property name for interpolation function h(xi).");

  // Constants
  params.addRequiredParam<Real>("Omega_v", "Partial molar volume Omega_v.");
  params.addRequiredParam<Real>("Omega_L", "Reference molar volume Omega_L.");

  return params;
}

LiVoidElasticAllenCahn::LiVoidElasticAllenCahn(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _psi_e(getMaterialProperty<Real>("elastic_energy")),
    _hprime(getMaterialPropertyDerivative<Real>("h_name", _var.name())),
    _h2prime(getMaterialPropertyDerivative<Real>("h_name", _var.name(), _var.name())),
    _Omega_v(getParam<Real>("Omega_v")),
    _Omega_L(getParam<Real>("Omega_L"))
{
}

void
LiVoidElasticAllenCahn::initialSetup()
{
  // Let base class do its checks for mobility L
  ACBulk<Real>::initialSetup();

  // Additionally check h(xi) has required derivatives
  validateDerivativeMaterialPropertyBase<Real>("h_name");
}

/**
 * Compute dF/dxi or d/dxi (dF/dxi)
 *
 * dF/dxi = (Omega_v / Omega_L) * h'(xi) * psi_e
 *
 * For the Jacobian we approximate:
 *   d/dxi(dF/dxi) ≈ (Omega_v / Omega_L) * h''(xi) * psi_e
 * treating psi_e as independent of xi in this kernel.
 */
Real
LiVoidElasticAllenCahn::computeDFDOP(PFFunctionType type)
{
  const Real A = _Omega_v / _Omega_L;

  switch (type)
  {
    case Residual:
      // dF/dxi
      return A * _hprime[_qp] * _psi_e[_qp];

    case Jacobian:
      // d/dxi(dF/dxi) * phi_j
      return A * _h2prime[_qp] * _psi_e[_qp] * _phi[_j][_qp];

    default:
      mooseError("Invalid PFFunctionType in LiVoidElasticAllenCahn");
  }
}
