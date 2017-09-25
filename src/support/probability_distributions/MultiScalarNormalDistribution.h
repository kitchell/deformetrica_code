/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _MultiScalarNormalDistribution_h
#define _MultiScalarNormalDistribution_h

/// Class file.
#include "AbstractNormalDistribution.h"

/**
 *  \brief      MultiScalarNormalDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw independent samples from normal distributions of equal variance.
 */
template<class ScalarType>
class MultiScalarNormalDistribution : public AbstractNormalDistribution<ScalarType> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract probability distribution type.
  typedef AbstractNormalDistribution<ScalarType> Superclass;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructors.
  MultiScalarNormalDistribution();

  /// Copy constructor.
  MultiScalarNormalDistribution(const MultiScalarNormalDistribution &other);

  /// Makes a copy of the object.
  virtual MultiScalarNormalDistribution *Clone() { return new MultiScalarNormalDistribution(*this); }

  /// Destructor.
  virtual ~MultiScalarNormalDistribution();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the covariance matrix of the normal distribution.
  virtual MatrixType GetCovariance() const {
    return diagonal_matrix<ScalarType>(Superclass::m_Mean.size(), m_Variance);
  }
  /// Sets the covariance matrix of the normal distribution.
  virtual void SetCovariance(MatrixType const &cov) { m_Variance = cov(0, 0); }

  /// Returns the cholesky squared root of the covariance matrix.
  virtual MatrixType GetCovarianceSqrt() const {
    return diagonal_matrix<ScalarType>(Superclass::m_Mean.size(), m_VarianceSqrt);
  }
  /// Sets the cholesky squared root of the covariance matrix.
  virtual void SetCovarianceSqrt(MatrixType const &covSqrt) { m_VarianceSqrt = covSqrt(0, 0); }

  /// Returns the inverse of the covariance matrix.
  virtual MatrixType GetCovarianceInverse() const {
    return diagonal_matrix<ScalarType>(Superclass::m_Mean.size(), m_VarianceInverse);
  }
  /// Sets the the inverse of the covariance matrix.
  virtual void SetCovarianceInverse(MatrixType const &covInv) { m_VarianceInverse = covInv(0, 0); }

  /// Returns the variance of the normal distribution.
  ScalarType GetVariance() const { return m_Variance; }
  /// Sets the covariance matrix of the normal distribution.
  void SetVariance(ScalarType const &var) {
    m_Variance = var;
    m_VarianceSqrt = std::sqrt(var);
    m_VarianceInverse = 1.0 / var;
    m_VarianceLog = std::log(var);
  }

  /// Returns the squared root of the variance.
  ScalarType GetVarianceSqrt() const { return m_VarianceSqrt; }
  /// Sets the cholesky squared root of the covariance matrix.
  void SetVarianceSqrt(ScalarType const &varSqrt) {
    m_Variance = pow(varSqrt, 2);
    m_VarianceSqrt = varSqrt;
    m_VarianceInverse = 1.0 / m_Variance;
    m_VarianceLog = std::log(m_Variance);
  }

  /// Returns the inverse of the variance.
  ScalarType GetVarianceInverse() const { return m_VarianceInverse; }
  /// Sets the the inverse of the covariance matrix.
  void SetVarianceInverse(ScalarType const &varInv) {
    m_Variance = 1.0 / varInv;
    m_VarianceSqrt = std::sqrt(m_Variance);
    m_VarianceInverse = varInv;
    m_VarianceLog = std::log(m_Variance);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples from the normal distribution to generate a realization.
  virtual VectorType Sample() const;

  /// Computes the log-likelihood of an observation.
  virtual ScalarType ComputeLogLikelihood(VectorType const &obs) const {
    return ComputeLogLikelihood(obs, 1.0);
  }
  /// Variation which vectorizes the linear variable in a first step.
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs) const {
    return ComputeLogLikelihood(obs, 1.0);
  }
  /// Variation featuring a temperature parameter (for the McmcSaem algorithm). Vectorizes \e obs in a first step.
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs,
                                          ScalarType const &temperature) const {
    assert(Superclass::m_Mean.size() == obs.n_elem());
    VectorType aux = obs.vectorize();
    return ComputeLogLikelihood(aux, temperature);
  }
  /// Variation featuring a temperature parameter (for the McmcSaem algorithm).
  virtual ScalarType ComputeLogLikelihood(VectorType const &obs,
                                          ScalarType const &temperature) const;

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Variance of the independent normal distributions.
  ScalarType m_Variance;
  /// Standard deviation of the idependent normal distributions.
  ScalarType m_VarianceSqrt;
  /// Inverse of the variance of the independent normal distributions.
  ScalarType m_VarianceInverse;
  /// Logarithm of the variance of the independent normal distributions.
  ScalarType m_VarianceLog;

}; /* class MultiScalarNormalDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "MultiScalarNormalDistribution.txx"
//#endif

#endif /* _MultiScalarNormalDistribution_h */
