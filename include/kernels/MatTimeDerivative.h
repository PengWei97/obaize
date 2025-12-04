/**
 * This kernel extends the standard MOOSE TimeDerivative kernel by multiplying the
 * time-derivative term \f$ \dot{u} \f$ with a MaterialProperty.
 *
 * In weak form, the residual contribution is:
 *
 *   \f[
 *     R = \int_{\Omega} \eta \, M(\mathbf{x}) \, \dot{u} \, d\Omega,
 *   \f]
 *
 * where:
 * - \f$ u \f$ is the primary variable,
 * - \f$ \eta \f$ is the test function,
 * - \f$ M(\mathbf{x}) \f$ is a user-specified material property.
 *
 * This class originates from the XFEM module's TestMatTimeDerivative kernel
 * and is a general-purpose utility kernel whenever a coefficient is needed
 * in front of a time derivative term.
 *
 * Typical usage example:
 *
 *   - modeling \f$ h(\xi) \, \dot{\theta}_m \f$ in Li diffusion equations,
 *   - introducing a spatial/phase-field/temperature dependent scaling on \f$ \dot{u} \f$,
 *   - implementing position-dependent capacitance-like terms.
 *
 * The material property name is supplied through input parameters, and the kernel
 * automatically retrieves it as `_mat_prop_value`.
 */

#pragma once

#include "TimeDerivative.h"

class MatTimeDerivative : public TimeDerivative
{
public:
  static InputParameters validParams();

  MatTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  const MaterialProperty<Real> & _mat_prop_value;
};
