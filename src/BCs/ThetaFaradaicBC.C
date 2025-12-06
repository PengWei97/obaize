#include "ThetaFaradaicBC.h"

registerMooseObject("obaizeApp", ThetaFaradaicBC);

InputParameters
ThetaFaradaicBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();

  // Coupled electric potential
  params.addRequiredCoupledVar("phi", "Electric potential");

  // Electrical conductivity sigma (can be electrode or electrolyte,
  // depending on how the BC is applied)
  params.addRequiredParam<MaterialPropertyName>(
      "conductivity", "Name of conductivity material property");

  // Faraday constant [C/mol]
  params.addParam<Real>("F", 96485.0, "Faraday constant (C/mol)");

  // Li molar volume Ω_L [m^3/mol]; used to convert from molar flux
  // boundary condition to occupancy-rate boundary condition for theta_m
  params.addRequiredParam<Real>("Omega_L",
                                "Molar volume of Li, Omega_L [m^3/mol]");

  // Optional sign to flip the normal, if needed
  params.addParam<Real>("normal_sign", 1.0,
                        "Sign to flip the normal if needed (usually +1)");

  return params;
}

ThetaFaradaicBC::ThetaFaradaicBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _grad_phi(coupledGradient("phi")),
    _phi_var(coupled("phi")),
    _sigma(getMaterialProperty<Real>(
        getParam<MaterialPropertyName>("conductivity"))),
    _F(getParam<Real>("F")),
    _Omega_L(getParam<Real>("Omega_L")),
    _normal_sign(getParam<Real>("normal_sign"))
{
}

Real
ThetaFaradaicBC::computeQpResidual()
{
  // Outward normal from the primary side (electrode) into the electrolyte
  const RealVectorValue & n = _normals[_qp];

  // grad(phi) · n
  const Real gradphi_dot_n = _grad_phi[_qp] * n;

  // Electric current density normal component:
  // i_n = n · i = -sigma * (grad(phi) · n)
  const Real i_normal = -_sigma[_qp] * gradphi_dot_n;

  // According to the chemical balance after dividing by c_L^m:
  //  -j_tilde_Li · n = (i * Omega_L) / (z F), with z = 1 for Li.
  // Here we take normal_sign to allow flipping the overall sign
  // if the mesh normal orientation is opposite to the convention.
  const Real j_theta = _normal_sign * i_normal * _Omega_L / _F;

  // Weak form boundary term: \int_{Gamma_r} (i * Omega_L / F) * test dS
  return j_theta * _test[_i][_qp];
}

Real
ThetaFaradaicBC::computeQpJacobian()
{
  // For now, j_theta does not explicitly depend on theta_m
  return 0.0;
}

Real
ThetaFaradaicBC::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  // A more accurate implementation would add off-diagonal Jacobian
  // w.r.t. phi (since i_normal depends on grad(phi)), but we keep it
  // zero for now and rely on the non-linear iterations.
  return 0.0;
}
