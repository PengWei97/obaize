#include "LiMechCouplingKernel.h"

registerMooseObject("obaizeApp", LiMechCouplingKernel);

InputParameters
LiMechCouplingKernel::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<MaterialPropertyName>(
      "Deff", "Effective diffusion coefficient D_eff");
  params.addRequiredParam<MaterialPropertyName>(
      "h_xi", "Interpolation function h(xi)");
  params.addRequiredCoupledVar(
      "sigma_h_xi", "Hydrostatic stress field sigma_h^xi");

  params.addRequiredParam<Real>(
      "Omega_Li", "Partial molar volume of Li");
  params.addRequiredParam<Real>(
      "Omega_v", "Partial molar volume of vacancy");
  params.addRequiredParam<Real>(
      "R", "Gas constant");
  params.addRequiredParam<Real>(
      "T", "Absolute temperature");

  return params;
}

LiMechCouplingKernel::LiMechCouplingKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _Deff(getMaterialProperty<Real>("Deff")),
    _h_xi(getMaterialProperty<Real>("h_xi")),
    _grad_sigma_h_xi(coupledGradient("sigma_h_xi")),
    _sigma_var(coupled("sigma_h_xi")),
    _Omega_Li(getParam<Real>("Omega_Li")),
    _Omega_v(getParam<Real>("Omega_v")),
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T"))
{
}

Real
LiMechCouplingKernel::computeQpResidual()
{
  // coeff = - D_eff * h(xi) * theta_m * (Omega_Li - Omega_v) / (R T)
  const Real coeff =
      -_Deff[_qp] * _h_xi[_qp] * _u[_qp] * (_Omega_Li - _Omega_v) / (_R * _T);

  // Residual: -∫ D h theta (Ω_Li-Ω_v)/(RT) ∇sigma · ∇w
  return coeff * (_grad_sigma_h_xi[_qp] * _grad_test[_i][_qp]);
}

Real
LiMechCouplingKernel::computeQpJacobian()
{
  // dR/d theta_m: replace theta_m by basis function phi_j
  const Real coeff =
      -_Deff[_qp] * _h_xi[_qp] * (_Omega_Li - _Omega_v) / (_R * _T);

  return coeff * _phi[_j][_qp] * (_grad_sigma_h_xi[_qp] * _grad_test[_i][_qp]);
}

Real
LiMechCouplingKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Off-diagonal Jacobian with respect to sigma_h_xi
  if (jvar == _sigma_var)
  {
    // d/d sigma: grad_sigma -> grad_phi_j
    const Real coeff =
        -_Deff[_qp] * _h_xi[_qp] * _u[_qp] * (_Omega_Li - _Omega_v) / (_R * _T);

    return coeff * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
  }

  return 0.0;
}
