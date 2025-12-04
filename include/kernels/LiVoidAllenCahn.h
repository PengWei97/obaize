#pragma once

#include "ACBulk.h"

/**
 * LiVoidAllenCahn
 *
 * Implements the chemical driving term for the void phase-field evolution:
 *
 *   (eta_xi, (R*T/Ω_L) * h'(xi) * ln((1 - theta_m) / (1 - theta_m0)) )
 *
 * The full Allen–Cahn contribution handled by ACBulk is:
 *
 *   Residual:        L * dF/dxi
 *   Jacobian (diag): L * d^2F/dxi^2 + (dL/dxi) * dF/dxi
 *
 * Here we define:
 *
 *   dF/dxi = (R*T/Ω_L) * h'(xi) * ln((1 - theta_m)/(1 - theta_m0))
 *
 * h(xi) is provided by a DerivativeParsedMaterial with base property name "h_name".
 * Its first and second derivatives h'(xi), h''(xi) are automatically obtained via
 * Material derivatives.
 */
class LiVoidAllenCahn : public ACBulk<Real>
{
public:
  static InputParameters validParams();
  LiVoidAllenCahn(const InputParameters & parameters);

  virtual void initialSetup();
protected:
  /// Compute dF/dxi or its derivative w.r.t xi (for residual / Jacobian)
  virtual Real computeDFDOP(PFFunctionType type) override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  /// Coupled mobile Li fraction θ_m
  const VariableValue & _theta_m;

  /// h'(xi) from DerivativeMaterialInterface: d(h_name)/d(xi)
  const MaterialProperty<Real> & _hprime;

  /// h''(xi) from DerivativeMaterialInterface: d²(h_name)/d(xi)d(xi)
  const MaterialProperty<Real> & _h2prime;

  /// Constants from input
  const Real _R;
  const Real _T;
  const Real _Omega_L;
  const Real _theta_m0;
};
