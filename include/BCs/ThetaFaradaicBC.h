#pragma once

#include "IntegratedBC.h"

/**
 * Faradaic flux boundary condition for theta_m (Li occupancy)
 * on the electrode–electrolyte interface.
 *
 * The chemical balance after dividing by c_L^m gives the boundary flux:
 *
 *   - j_tilde_Li · n = (i * Omega_L) / (z F),     (z = 1 for Li)
 *
 * where
 *   i = -sigma * grad(phi)
 *
 * Therefore the boundary contribution to the weak form is:
 *
 *   R_i =  (normal_sign * i * Omega_L / F) * test_theta_i
 *
 * This BC enforces the conversion between electric current crossing
 * the interface and the change of Li lattice-site occupancy theta_m.
 *
 * NOTE:
 *   Currently computeQpJacobian() and computeQpOffDiagJacobian()
 *   return 0. While a full Jacobian would include terms wrt phi,
 *   this simplified version works well for moderate coupling.
 */
class ThetaFaradaicBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  ThetaFaradaicBC(const InputParameters & parameters);

protected:
  /// Residual contribution at quadrature point
  virtual Real computeQpResidual() override;

  /// Diagonal Jacobian wrt theta_m (currently zero)
  virtual Real computeQpJacobian() override;

  /// Off-diagonal Jacobian wrt other variables (phi, etc.) (currently zero)
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // ----------------------------------------------------
  // Coupled fields & material properties
  // ----------------------------------------------------

  /// ∇phi (electric potential gradient)
  const VariableGradient & _grad_phi;

  /// Variable index of phi (unused here but kept for consistency)
  const unsigned int _phi_var;

  /// Electrical conductivity sigma(x, xi)
  const MaterialProperty<Real> & _sigma;

  // ----------------------------------------------------
  // Physical constants
  // ----------------------------------------------------

  /// Faraday constant [C/mol]
  const Real _F;

  /// Li molar volume Omega_L [m^3/mol], needed for boundary flux conversion
  const Real _Omega_L;

  /// Optional sign to flip the normal direction (+1 or -1)
  const Real _normal_sign;
};
