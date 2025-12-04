#pragma once

#include "CoupledTimeDerivative.h"

class MatCoupledTimeDerivative : public CoupledTimeDerivative
{
public:
  static InputParameters validParams();
  MatCoupledTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const MaterialProperty<Real> & _mat_prop_value;
};
