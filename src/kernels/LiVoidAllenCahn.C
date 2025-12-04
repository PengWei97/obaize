#include "LiVoidAllenCahn.h"

registerMooseObject("obaizeApp", LiVoidAllenCahn);

InputParameters
LiVoidAllenCahn::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();

  params.addClassDescription(
      "Allen–Cahn chemical driving term for Li–void phase-field model, "
      "using DerivativeParsedMaterial for h(xi).");

  // Couple theta_m
  params.addRequiredCoupledVar("theta_m", "Mobile lithium fraction (theta_m)");

  // Base name of interpolation function h(xi)
  params.addRequiredParam<MaterialPropertyName>(
      "h_name", "Base material property name for interpolation function h(xi)");

  // Constants
  params.addRequiredParam<Real>("R", "Gas constant");
  params.addRequiredParam<Real>("T", "Absolute temperature");
  params.addRequiredParam<Real>("Omega_L", "Reference molar volume Ω_L");
  params.addRequiredParam<Real>("theta_m0", "Reference mobile Li fraction theta_m0");

  return params;
}

LiVoidAllenCahn::LiVoidAllenCahn(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _theta_m(coupledValue("theta_m")),
    _hprime(getMaterialPropertyDerivative<Real>("h_name", _var.name())),
    _h2prime(getMaterialPropertyDerivative<Real>("h_name", _var.name(), _var.name())),

    // Constants
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),
    _Omega_L(getParam<Real>("Omega_L")),
    _theta_m0(getParam<Real>("theta_m0"))
{
}

void
LiVoidAllenCahn::initialSetup()
{
  // Let base class do its checks for mobility L
  ACBulk<Real>::initialSetup();

  // Additionally check h(xi) has required derivatives
  validateDerivativeMaterialPropertyBase<Real>("h_name");
}

/**
 * Compute dF/dxi or d/dxi (dF/dxi)
 *
 * dF/dxi = (R*T/Ω_L) * h'(xi) * ln((1 - theta_m) / (1 - theta_m0))
 */
Real
LiVoidAllenCahn::computeDFDOP(PFFunctionType type)
{
  const Real A    = _R * _T / _Omega_L;
  const Real eps  = 1e-8;
  const Real tiny = 1e-14;

  // 1 - theta_m0 must stay positive
  const Real denom0 = 1.0 - _theta_m0;
  if (denom0 <= 0.0)
    mooseError("In LiVoidAllenCahn: 1 - theta_m0 must be positive.");

  Real th = _theta_m[_qp];

  // If theta_m is already NaN/Inf, do not contribute to residual/Jacobian
  if (!std::isfinite(th))
    return 0.0;

  // Clamp theta_m to avoid 1 - theta_m <= 0
  if (th >= 1.0 - eps)
    th = 1.0 - eps;

  const Real denom = 1.0 - th;
  const Real ratio = std::max(denom / denom0, tiny);
  const Real log_term = std::log(ratio);

  switch (type)
  {
    case Residual:
      // dF/dxi
      return A * _hprime[_qp] * log_term;

    case Jacobian:
      // d/dxi(dF/dxi) ≈ A * h''(xi) * ln(...)
      return A * _h2prime[_qp] * log_term;

    default:
      mooseError("Invalid PFFunctionType in LiVoidAllenCahn");
  }
}
  
/**
 * Off-diagonal Jacobian contribution from dependency on theta_m.
 *
 * For F = (R*T/Ω_L) h'(xi) ln[(1 - theta_m)/(1 - theta_m0)],
 *
 *   ∂F/∂theta_m = (R*T/Ω_L) h'(xi) * ( -1 / (1 - theta_m) )
 *
 * The contribution to the residual is L * F, so:
 *
 *   dR_i / dtheta_j = ∫ L * ∂F/∂theta_m * phi_j * test_i dV
 */
Real
LiVoidAllenCahn::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Start with base-class contribution (e.g. for other coupled variables)
  Real val = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  // Only add extra coupling for theta_m
  if (jvar != coupled("theta_m"))
    return val;

  const Real A   = _R * _T / _Omega_L;
  const Real eps = 1e-8;

  Real th = _theta_m[_qp];

  // If theta_m is already NaN/Inf, skip off-diagonal contribution
  if (!std::isfinite(th))
    return val;

  // Clamp theta_m to avoid division by zero in 1/(1 - theta_m)
  if (th >= 1.0 - eps)
    th = 1.0 - eps;

  const Real denom = 1.0 - th;

  // ∂F/∂theta_m
  const Real dFd_dtheta = A * _hprime[_qp] * (-1.0 / denom);

  // Add contribution: L * ∂F/∂theta_m * phi_j * test_i
  val += _L[_qp] * dFd_dtheta * _phi[_j][_qp] * _test[_i][_qp];

  return val;
}