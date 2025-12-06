#include "CurrentMagnitudeAux.h"

// 把 "ObaizeApp" 换成你真实的 app 名（和其它 registerMooseObject 一样）
registerMooseObject("obaizeApp", CurrentMagnitudeAux);

InputParameters
CurrentMagnitudeAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("phi", "Electric potential");
  params.addRequiredParam<MaterialPropertyName>("conductivity",
                                                "Name of conductivity material property");
  return params;
}

CurrentMagnitudeAux::CurrentMagnitudeAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_phi(coupledGradient("phi")),
    _sigma(getMaterialProperty<Real>(getParam<MaterialPropertyName>("conductivity")))
{
}

Real
CurrentMagnitudeAux::computeValue()
{
  // grad(phi) at this quadrature point
  const RealVectorValue & grad_phi = _grad_phi[_qp];

  // current density vector: i = -sigma * grad(phi)
  RealVectorValue i_vec = -_sigma[_qp] * grad_phi;

  // return its magnitude
  return i_vec.norm();
}
