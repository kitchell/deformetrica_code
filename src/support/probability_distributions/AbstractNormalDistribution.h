/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AbstractNormalDistribution_h
#define _AbstractNormalDistribution_h

/// Class file.
#include "AbstractProbabilityDistribution.h"

/**
 *  \brief      AbstractNormalDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a normal distribution, determined by its mean and variance.
 */
template<class ScalarType>
class AbstractNormalDistribution : public AbstractProbabilityDistribution<ScalarType> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract probability distribution type.
  typedef AbstractProbabilityDistribution<ScalarType> Superclass;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructors.
  AbstractNormalDistribution(); // Standard normal distribution.

  /// Copy constructor.
  AbstractNormalDistribution(const AbstractNormalDistribution &other);

  /// Makes a copy of the object.
  virtual AbstractNormalDistribution *Clone() = 0;

  /// Destructor.
  virtual ~AbstractNormalDistribution();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the mean vector of the normal distribution.
  VectorType GetMean() const { return m_Mean; }
  /// Sets the mean vector of the normal distribution, for different input data structures.
  virtual void SetMean(ScalarType const &mean) { m_Mean = VectorType(1, mean); }
  virtual void SetMean(VectorType const &mean) { m_Mean = mean; }
  virtual void SetMean(MatrixType const &mean) { m_Mean = mean.vectorize(); }
  virtual void SetMean(MatrixListType const &mean) { m_Mean = mean.vectorize(); }
  virtual void SetMean(LinearVariableType const &mean) { m_Mean = mean.vectorize(); }

  /// Returns the covariance matrix of the normal distribution.
  virtual MatrixType GetCovariance() const = 0;
  /// Sets the covariance matrix of the normal distribution.
  virtual void SetCovariance(MatrixType const &cov) = 0;

  /// Returns the cholesky squared root of the covariance matrix.
  virtual MatrixType GetCovarianceSqrt() const = 0;
  /// Sets the cholesky squared root of the covariance matrix.
  virtual void SetCovarianceSqrt(MatrixType const &covSqrt) = 0;

  /// Returns the inverse of the covariance matrix.
  virtual MatrixType GetCovarianceInverse() const = 0;
  /// Sets the the inverse of the covariance matrix.
  virtual void SetCovarianceInverse(MatrixType const &covInv) = 0;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples from the normal distribution to generate a realization.
  virtual VectorType Sample() const = 0;

  /// Computes the log-likelihood of an observation. Vectorizes the linear variable in a first step.
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs) const = 0;
  /// Computes the log-likelihood of an observation.
  virtual ScalarType ComputeLogLikelihood(VectorType const &obs) const = 0;

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Mean of the normal distribution.
  VectorType m_Mean;

}; /* class AbstractNormalDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "AbstractNormalDistribution.txx"
//#endif

#endif /* _AbstractNormalDistribution_h */
