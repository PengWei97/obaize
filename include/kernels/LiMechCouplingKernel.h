#pragma once

#include "Kernel.h"

/**
 * LiMechCouplingKernel
 *
 * Contribution of the mechano-chemical coupling term:
 *   -∫ D_eff h(xi) theta_m (Omega_Li - Omega_v) / (R T) ∇sigma_h^xi · ∇w dΩ
 *
 * Primary variable: theta_m
 * Coupled variable: sigma_h_xi
 */
class LiMechCouplingKernel : public Kernel
{
public:
  static InputParameters validParams();
  LiMechCouplingKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Effective diffusion coefficient D_eff
  const MaterialProperty<Real> & _Deff;

  /// Interpolation function h(xi)
  const MaterialProperty<Real> & _h_xi;

  /// Gradient of coupled hydrostatic stress sigma_h^xi
  const VariableGradient & _grad_sigma_h_xi;

  /// Variable id for sigma_h_xi (for off-diagonal Jacobian)
  const unsigned int _sigma_var;

  /// Molar volumes and thermal constants
  const Real _Omega_Li;
  const Real _Omega_v;
  const Real _R;
  const Real _T;
};
