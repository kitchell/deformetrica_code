/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah. All rights reserved. This file is     *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

/// Class file.
#include "AbstractSampler.h"

/// Support files.
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"

using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      Anisotropic Metropolis Adjusted Langevin Algorithm.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A sampler simulates a random variable according to a target distribution.
 */

template<class ScalarType, unsigned int Dimension>
class AmalaSampler : public AbstractSampler<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract sampler type.
  typedef AbstractSampler<ScalarType, Dimension> Superclass;

  /// Abstract statistical model type.
  typedef typename Superclass::StatisticalModelType StatisticalModelType;
  /// Longitudinal data set type.
  typedef typename Superclass::LongitudinalDataSetType LongitudinalDataSetType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  AmalaSampler();

  /// Copy constructor.
  AmalaSampler(const AmalaSampler &other);

  /// Makes a copy of the object.
  virtual AmalaSampler *Clone() { return new AmalaSampler(*this); }

  /// Destructor.
  virtual ~AmalaSampler();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the truncation threshold.
  void Sethreshold(const ScalarType &s) { m_Threshold = s; }

  /// Sets the stochastic step sizes for the mean.
  void SetMeanScales(const std::map<std::string, ScalarType> &s) { m_MeanScales = s; }
  /// Sets the stochastic step sizes for the covariance.
  void SetCovarianceScales(const std::map<std::string, ScalarType> &s) { m_CovarianceScales = s; }

  /// Sets the regularization parameter.
  void SetRegularizations(const std::map<std::string, ScalarType> &s) { m_Regularizations = s; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples.
  virtual void Sample(LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER,
                      VectorType &acceptanceRates);

  /// Initializes the sampler.
  virtual void Initialize(LinearVariableMapType const &popRER,
                          LinearVariablesMapType const &indRER);

  /// Adapts the proposal distributions based on the detected acceptance rates and the iteration number.
  virtual void AdaptProposalDistributions(const VectorType &detectedAcceptanceRates,
                                          const unsigned int &iterationNumber,
                                          const bool verbose = false);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the Amala proposal distribution at point z = (popRER, indRER).
  void ComputeAmalaProposalDistribution(LinearVariableMapType const &popRER,
                                        LinearVariablesMapType const &indRER,
                                        VectorType &totRER,
                                        NormalDistributionType &prop) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Truncation threshold. Provides more stability to the Amala sampler.
  ScalarType m_Threshold;
  /// Stochastic step sizes, which control the mean of the proposal distribution.
  std::map<std::string, ScalarType> m_MeanScales;
  /// Stochastic step sizes, which control the covariance of the proposal distribution.
  std::map<std::string, ScalarType> m_CovarianceScales;
  /// Regularization parameters, which ensures that the proposal distribution covariance matrix is sympd.
  std::map<std::string, ScalarType> m_Regularizations;

  /// Memory of the current random effects realization vector.
  VectorType m_CurrentTotalRER;
  /// Memory of the current proposal distribution.
  NormalDistributionType m_CurrentProposalRED;

  /// Size parameters of the random effects realizations.
  std::vector<std::vector<unsigned int>> m_SizeParametersRER;
  /// Total number of scalars in the random effects realizations.
  unsigned int m_TotalSizeRER;

  /// For bug-tracking.
  static unsigned int m_Count;

}; /* class AmalaSampler */

