//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MatCoupledTimeDerivative.h"

registerMooseObject("obaizeApp", MatCoupledTimeDerivative);

InputParameters
MatCoupledTimeDerivative::validParams()
{
  InputParameters params = CoupledTimeDerivative::validParams();
  params.addRequiredParam<MaterialPropertyName>("mat_prop_value",
                                                "Material "
                                                "property to multiply by time "
                                                "derivative");
  return params;
}

MatCoupledTimeDerivative::MatCoupledTimeDerivative(const InputParameters & parameters)
  : CoupledTimeDerivative(parameters), _mat_prop_value(getMaterialProperty<Real>("mat_prop_value"))
{
}

Real
MatCoupledTimeDerivative::computeQpResidual()
{
  return _mat_prop_value[_qp] * CoupledTimeDerivative::computeQpResidual();
}

Real
MatCoupledTimeDerivative::computeQpJacobian()
{
  return _mat_prop_value[_qp] * CoupledTimeDerivative::computeQpJacobian();
}

Real
MatCoupledTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_var)
    return _mat_prop_value[_qp] *
           CoupledTimeDerivative::computeQpOffDiagJacobian(jvar);

  return 0.0;
}