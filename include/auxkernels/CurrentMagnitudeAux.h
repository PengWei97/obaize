#pragma once

#include "AuxKernel.h"

/**
 * Compute the magnitude of the current density:
 *
 *   i = -sigma * grad(phi)
 *   |i| = sqrt( i_x^2 + i_y^2 + i_z^2 )
 *
 * and store it in a scalar AuxVariable.
 */
class CurrentMagnitudeAux : public AuxKernel
{
public:
  static InputParameters validParams();
  CurrentMagnitudeAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Coupled electric potential
  const VariableGradient & _grad_phi;

  /// Conductivity sigma(x, xi) from Material
  const MaterialProperty<Real> & _sigma;
};
