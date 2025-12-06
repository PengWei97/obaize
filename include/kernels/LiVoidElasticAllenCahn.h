#pragma once

#include "ACBulk.h"

/**
 * LiVoidElasticAllenCahn
 *
 * Elastic part of the chemical potential for the void phase-field evolution:
 *
 *   dF/dxi = (Omega_v / Omega_L) * h'(xi) * psi_e
 *
 * The full Allen–Cahn contribution handled by ACBulk is:
 *
 *   Residual:        L * dF/dxi
 *   Jacobian (diag): L * d^2F/dxi^2 + (dL/dxi) * dF/dxi
 *
 * Here we define:
 *
 *   dF/dxi = (Omega_v / Omega_L) * h'(xi) * psi_e
 *
 * h(xi) is provided by a DerivativeParsedMaterial with base property name "h_name".
 * psi_e is the elastic strain energy density, typically from ComputePFElasticEnergy.
 */
class LiVoidElasticAllenCahn : public ACBulk<Real>
{
public:
  static InputParameters validParams();
  LiVoidElasticAllenCahn(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  /// Compute dF/dxi or its derivative w.r.t xi (for residual / Jacobian)
  virtual Real computeDFDOP(PFFunctionType type) override;

  /// Elastic strain energy density psi_e
  const MaterialProperty<Real> & _psi_e;

  /// h'(xi) from DerivativeMaterialInterface: d(h_name)/d(xi)
  const MaterialProperty<Real> & _hprime;

  /// h''(xi) from DerivativeMaterialInterface: d²(h_name)/d(xi)d(xi)
  const MaterialProperty<Real> & _h2prime;

  /// Constants
  const Real _Omega_v;
  const Real _Omega_L;
};
