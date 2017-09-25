/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DisplacementFieldNormalDistribution_h
#define _DisplacementFieldNormalDistribution_h

/// Class file.
#include "AbstractNormalDistribution.h"

/// Support files.
#include "KernelFactory.h"
#include "MultiScalarNormalDistribution.h"

/**
 *  \brief      DisplacementFieldNormalDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw independent samples from normal distributions of equal variance.
 */
template<class ScalarType>
class DisplacementFieldNormalDistribution : public AbstractNormalDistribution<ScalarType> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract probability distribution type.
  typedef AbstractNormalDistribution<ScalarType> Superclass;


  /// Multi scalar normal distribution type.
  typedef MultiScalarNormalDistribution<ScalarType> MultiScalarNormalDistributionType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructors.
  DisplacementFieldNormalDistribution();

  /// Copy constructor.
  DisplacementFieldNormalDistribution(const DisplacementFieldNormalDistribution &other);

  /// Makes a copy of the object.
  virtual DisplacementFieldNormalDistribution *Clone() { return new DisplacementFieldNormalDistribution(*this); }

  /// Destructor.
  virtual ~DisplacementFieldNormalDistribution();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the mean vector of the normal distribution, for different input data structures.
  virtual void SetMean(ScalarType const &mean) {
    Superclass::SetMean(mean);
    m_MomentaDistribution.SetMean(VectorType(Superclass::m_Mean.size(), 0.0));
  }
  virtual void SetMean(VectorType const &mean) {
    Superclass::SetMean(mean);
    m_MomentaDistribution.SetMean(VectorType(Superclass::m_Mean.size(), 0.0));
  }
  virtual void SetMean(MatrixType const &mean) {
    Superclass::SetMean(mean);
    m_MomentaDistribution.SetMean(VectorType(Superclass::m_Mean.size(), 0.0));
  }
  virtual void SetMean(MatrixListType const &mean) {
    Superclass::SetMean(mean);
    m_MomentaDistribution.SetMean(VectorType(Superclass::m_Mean.size(), 0.0));
  }
  virtual void SetMean(LinearVariableType const &mean) {
    Superclass::SetMean(mean);
    m_MomentaDistribution.SetMean(VectorType(Superclass::m_Mean.size(), 0.0));
  }

  /// Returns the covariance matrix of the normal distribution.
  virtual MatrixType GetCovariance() const { return m_MomentaDistribution.GetCovariance(); }
  /// Sets the covariance matrix of the normal distribution.
  virtual void SetCovariance(MatrixType const &cov) { m_MomentaDistribution.SetCovariance(cov); }

  /// Returns the cholesky squared root of the covariance matrix.
  virtual MatrixType GetCovarianceSqrt() const { return m_MomentaDistribution.GetCovarianceSqrt(); }
  /// Sets the cholesky squared root of the covariance matrix.
  virtual void SetCovarianceSqrt(MatrixType const &covSqrt) { m_MomentaDistribution.SetCovarianceSqrt(covSqrt); }

  /// Returns the inverse of the covariance matrix.
  virtual MatrixType GetCovarianceInverse() const { return m_MomentaDistribution.GetCovarianceInverse(); }
  /// Sets the the inverse of the covariance matrix.
  virtual void SetCovarianceInverse(MatrixType const &covInv) { m_MomentaDistribution.SetCovarianceInverse(covInv); }

  /// Returns the variance of the normal distribution.
  ScalarType GetVariance() const { return m_MomentaDistribution.GetVariance(); }
  /// Sets the covariance matrix of the normal distribution.
  void SetVariance(ScalarType const &var) { m_MomentaDistribution.SetVariance(var); }

  /// Returns the squared root of the variance.
  ScalarType GetVarianceSqrt() const { return m_MomentaDistribution.GetVarianceSqrt(); }
  /// Sets the cholesky squared root of the covariance matrix.
  void SetVarianceSqrt(ScalarType const &varSqrt) { m_MomentaDistribution.SetVarianceSqrt(varSqrt); }

  /// Returns the inverse of the variance.
  ScalarType GetVarianceInverse() const { return m_MomentaDistribution.GetVarianceInverse(); }
  /// Sets the the inverse of the covariance matrix.
  void SetVarianceInverse(ScalarType const &varInv) { m_MomentaDistribution.SetVarianceInverse(varInv); }

  /// Sets the shape dimension.
  void SetKernelType(const KernelEnumType type) { m_KernelType = type; }
  /// Sets the shape dimension.
  void SetKernelWidth(const ScalarType kw) { m_KernelWidth = kw; }
  /// Sets the shape dimension.
  void SetSubsamplingStepSize(const unsigned int sss) { m_SubsamplingStepSize = sss; }
  /// Sets the shape dimension.
  void SetShapeDimension(const unsigned int sd) { m_ShapeDimension = sd; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples from the normal distribution to generate a realization.
  virtual VectorType Sample() const;

  /// Computes the log-likelihood of an observation. Vectorizes the linear variable in a first step.
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs) const {
    std::cerr << "Error : DisplacementFieldNormalDistribution::ComputeLogLikelihood method is not available yet."
              << std::endl;
    return 0;
  }
  /// Computes the log-likelihood of an observation.
  virtual ScalarType ComputeLogLikelihood(VectorType const &obs) const {
    std::cerr << "Error : DisplacementFieldNormalDistribution::ComputeLogLikelihood method is not available yet."
              << std::endl;
    return 0;
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Distribution of the momenta.
  MultiScalarNormalDistributionType m_MomentaDistribution;

  /// Type of the kernel.
  KernelEnumType m_KernelType;
  /// Size of the kernel.
  ScalarType m_KernelWidth;
  /// Subsampling step size.
  unsigned int m_SubsamplingStepSize;
  /// Dimension of the shape.
  unsigned int m_ShapeDimension;

}; /* class DisplacementFieldNormalDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "DisplacementFieldNormalDistribution.txx"
//#endif

#endif /* _DisplacementFieldNormalDistribution_h */
