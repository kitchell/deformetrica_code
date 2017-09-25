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
#include "ProbabilityDistributions.h"


using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      Symmetric random walk Metropolis Hasting within Gibbs sampler object class.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A sampler simulates a random variable according to a target distribution.
 */

template<class ScalarType, unsigned int Dimension>
class SrwMhwgSampler : public AbstractSampler<ScalarType, Dimension> {
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
  SrwMhwgSampler();

  /// Copy constructor.
  SrwMhwgSampler(const SrwMhwgSampler &other);

  /// Makes a copy of the object.
  virtual SrwMhwgSampler *Clone() { return new SrwMhwgSampler(*this); }

  /// Destructor.
  virtual ~SrwMhwgSampler();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the population proposal distributions.
  void SetPopulationProposalDistribution(const AbstractNormalDistributionMapType &map) {
    m_PopulationProposalDistribution = map;
  }
  void SetPopulationProposalDistribution(const MultiScalarNormalDistributionMapType &map) {
    m_PopulationProposalDistribution.clear();
    for (auto it = map.begin(); it != map.end(); ++it) {
      m_PopulationProposalDistribution[it->first] =
          std::static_pointer_cast<AbstractNormalDistributionType>(it->second);
    }
  }

  /// Sets the individual proposal distributions.
  void SetIndividualProposalDistribution(const AbstractNormalDistributionMapType &map) {
    m_IndividualProposalDistribution = map;
  }
  void SetIndividualProposalDistribution(const MultiScalarNormalDistributionMapType &map) {
    m_IndividualProposalDistribution.clear();
    for (auto it = map.begin(); it != map.end(); ++it) {
      m_IndividualProposalDistribution[it->first] =
          std::static_pointer_cast<AbstractNormalDistributionType>(it->second);
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples. ("RER" Random Effects Realization).
  virtual void Sample(LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER,
                      VectorType &acceptanceRates) {
    Sample(popRER, indRER, acceptanceRates, 1.0);
  }
  /// Variation with a temperature parameter.
  virtual void Sample(LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER,
                      VectorType &acceptanceRates,
                      const ScalarType &temperature);

  /// Initializes the sampler. Empty function for the Srwhwg sampler.
  virtual void Initialize(LinearVariableMapType const &popRER,
                          LinearVariablesMapType const &indRER) {};

  /// Adapts the proposal distributions based on the detected acceptance rates and the iteration number.
  virtual void AdaptProposalDistributions(const VectorType &detectedAcceptanceRates,
                                          const unsigned int &iterationNumber,
                                          const bool verbose = false);

  /// Saves the current state of the sampler.
  virtual void SaveState(def::utils::DeformationState &state) const;
  /// Recovers a previously saved state.
  virtual void RecoverState(def::utils::DeformationState &state);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Population proposal distributions. ("PD" Proposal Distribution).
  AbstractNormalDistributionMapType m_PopulationProposalDistribution;
  /// Individual proposal distributions. ("PD" Proposal Distribution).
  AbstractNormalDistributionMapType m_IndividualProposalDistribution;

}; /* class SrwMhwgSampler */

