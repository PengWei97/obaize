#pragma once

#include "ComputeElasticityTensorBase.h"

/**
 * Isotropic elasticity tensor degraded by a phase-field interpolation h(xi):
 *   G^xi = h(xi) * G,  K^xi = h(xi) * K.
 */
class ComputePFIsotropicElasticityTensor : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();

  ComputePFIsotropicElasticityTensor(const InputParameters & parameters);

  virtual void residualSetup() override;
  virtual void computeQpElasticityTensor() override;

protected:
  /// Reference bulk and shear moduli (intact material)
  const Real & _bulk_modulus;
  const Real & _shear_modulus;

  /// Interpolation function h(xi)
  const MaterialProperty<Real> & _h_xi;

  Real _lambda0;
  Real _E0;
  Real _nu0;
  Real _effective_stiffness_local;

  Real _kappa0;   // <<--- new: residual stiffness fraction

  RankFourTensor _Cijkl;

  using ComputeElasticityTensorBase::name;
  using ComputeElasticityTensorBase::_elasticity_tensor_name;
  using ComputeElasticityTensorBase::issueGuarantee;
  using ComputeElasticityTensorBase::_elasticity_tensor;
  using ComputeElasticityTensorBase::_effective_stiffness;
  using ComputeElasticityTensorBase::_qp;
};
