#include "ChemicalEigenstrain.h"

registerMooseObject("obaizeApp", ChemicalEigenstrain);

InputParameters
ChemicalEigenstrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();

  params.addClassDescription(
      "Computes isotropic chemical eigenstrain for Li metal "
      "based on phase-field xi and Li occupancy theta_m.");

  params.addRequiredCoupledVar("xi", "Phase-field variable (0=void, 1=metal)");
  params.addRequiredCoupledVar("theta_m", "Li lattice occupancy");

  params.addRequiredParam<Real>("theta_m0",
        "Reference Li occupancy for stress-free state");

  params.addRequiredParam<Real>("Omega_L",
        "Molar volume of Li lattice sites [m^3/mol]");

  params.addRequiredParam<Real>("Omega_v",
        "Molar volume of vacancies [m^3/mol]");

  params.addRequiredParam<MaterialPropertyName>("h_xi_name",
        "Material property name of h(xi)");

  return params;
}

ChemicalEigenstrain::ChemicalEigenstrain(const InputParameters & parameters)
  : ComputeEigenstrainBase(parameters),

    // ⭐ Correct way: materials use coupledValue(), NOT getVar()
    _xi(coupledValue("xi")),
    _theta_m(coupledValue("theta_m")),

    _theta_m0(getParam<Real>("theta_m0")),
    _Omega_L(getParam<Real>("Omega_L")),
    _Omega_v(getParam<Real>("Omega_v")),

    _c_Lm(1.0 / _Omega_L),
    _Omega_diff(_Omega_L - _Omega_v),

    // ⭐ Correct way to retrieve a material property (non-AD version)
    _h_xi(getMaterialProperty<Real>(getParam<MaterialPropertyName>("h_xi_name")))
{
}

void
ChemicalEigenstrain::computeQpEigenstrain()
{
  const Real h_val       = _h_xi[_qp];
  const Real theta_diff  = _theta_m[_qp] - _theta_m0;

  // Scalar chemical strain, Eqs. (23)-(24)
  const Real eps_c_scalar =
      (1.0/3.0) * h_val * _c_Lm * theta_diff * _Omega_diff;

  this->_eigenstrain[_qp].zero();

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    this->_eigenstrain[_qp](i, i) = eps_c_scalar;
}
