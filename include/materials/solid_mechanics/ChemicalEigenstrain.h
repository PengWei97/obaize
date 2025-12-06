//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1

#pragma once

#include "ComputeEigenstrainBase.h"

/**
 * ChemicalEigenstrain
 *
 * Computes isotropic chemical eigenstrain for Li metal, following Eqs. (23)-(24):
 *
 *   eps_c = (1/3)*h(xi)*c_Lm*(theta_m - theta_m0)*(Omega_L - Omega_v) * I
 *
 * with
 *    c_Lm       = 1 / Omega_L
 *    h(xi)      = material property (DerivativeParsedMaterial)
 *    theta_m    = Li lattice occupancy
 *    theta_m0   = reference occupancy
 *    Omega_L    = Li site molar volume
 *    Omega_v    = vacancy molar volume
 *
 * Eigenstrain is volumetric (isotropic). It vanishes automatically in voids
 * because h(0) = 0.
 */
class ChemicalEigenstrain : public ComputeEigenstrainBase
{
public:
  static InputParameters validParams();

  ChemicalEigenstrain(const InputParameters & parameters);

protected:
  /// Compute eigenstrain at the current quadrature point
  virtual void computeQpEigenstrain() override;

  // -------- Coupled variables ----------
  const VariableValue & _xi;
  const VariableValue & _theta_m;

  // -------- Parameters ----------
  const Real _theta_m0;
  const Real _Omega_L;
  const Real _Omega_v;

  const Real _c_Lm;        // 1 / Omega_L
  const Real _Omega_diff;  // Omega_L - Omega_v

  // -------- Material property h(xi) ----------
  const MaterialProperty<Real> & _h_xi;
};
