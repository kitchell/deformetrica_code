/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NormalDistribution_h
#define _NormalDistribution_h

/// Class file.
#include "AbstractNormalDistribution.h"

/**
 *  \brief      NormalDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a normal distribution, determined by its mean and variance.
 */
template<class ScalarType>
class NormalDistribution : public AbstractNormalDistribution<ScalarType> {
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
  NormalDistribution();

  /// Copy constructor.
  NormalDistribution(const NormalDistribution &other);

  /// Makes a copy of the object.
  virtual NormalDistribution *Clone() { return new NormalDistribution(*this); }

  /// Destructor.
  virtual ~NormalDistribution();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the covariance matrix of the normal distribution.
  virtual MatrixType GetCovariance() const { return m_Covariance; }
  /// Sets the covariance matrix of the normal distribution.
  virtual void SetCovariance(MatrixType const &cov) {
    m_Covariance = cov;
    m_CovarianceSqrt = chol(cov);
    m_CovarianceInverse = inverse_sympd(cov);
    m_CovarianceLogDeterminant = log_det(cov);
  }

  /// Returns the cholesky squared root of the covariance matrix.
  virtual MatrixType GetCovarianceSqrt() const { return m_CovarianceSqrt; }
  /// Sets the cholesky squared root of the covariance matrix.
  virtual void SetCovarianceSqrt(MatrixType const &covSqrt) {
    m_CovarianceSqrt = covSqrt;
    m_Covariance = covSqrt * covSqrt.transpose();
    m_CovarianceInverse = inverse_sympd(m_Covariance);
    m_CovarianceLogDeterminant = log_det(m_Covariance);
  }

  /// Returns the inverse of the covariance matrix.
  virtual MatrixType GetCovarianceInverse() const { return m_CovarianceInverse; }
  /// Sets the the inverse of the covariance matrix.
  virtual void SetCovarianceInverse(MatrixType const &covInv) {
    m_CovarianceInverse = covInv;
    m_Covariance = inverse_sympd(covInv);
    m_CovarianceSqrt = chol(m_Covariance);
    m_CovarianceLogDeterminant = log_det(m_Covariance);
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

  /// Covariance matrix of the normal distribution.
  MatrixType m_Covariance;
  /// Square root of the covariance matrix (Cholesky square root).
  MatrixType m_CovarianceSqrt;
  /// Inverse of the covariance matrix.
  MatrixType m_CovarianceInverse;
  /// Log-determinant of the covariance matrix.
  ScalarType m_CovarianceLogDeterminant;

}; /* class NormalDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "NormalDistribution.txx"
//#endif

#endif /* _NormalDistribution_h */
