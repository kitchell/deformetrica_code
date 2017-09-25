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

/// Core files.
#include "AbstractStatisticalModel.h"
#include "LongitudinalDataSet.h"

/// Support files.
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"
#include <src/support/utilities/GeneralSettings.h>

/// Standard files.
#include <lib/ThreadPool/ThreadPool.h>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <src/support/utilities/SerializeDeformationState.h>

using namespace def::algebra;
using namespace def::proba;
/**
 *	\brief      AbstractSampler object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A sampler simulates a random variable according to a target distribution.
 */
template<class ScalarType, unsigned int Dimension>
class AbstractSampler {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Typedefs :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Possible types of samplers.
  typedef enum {
    null,         // Null value.
    SrwMhwg,     // Symmetric random walk - Metropolis Hasting within Gibbs.
    Mala,        // Metropolis-adjusted Langevin algorithm.
    Amala,       // Anisotropic Metropolis-adjusted Langevin algorithm.
  } SamplerType;

  /// Abstract statistical model type.
  typedef AbstractStatisticalModel<ScalarType, Dimension> StatisticalModelType;

  /// Longitudinal dataset type.
  typedef LongitudinalDataSet<ScalarType, Dimension> LongitudinalDataSetType;

  /// Deformable multi object type
  typedef typename StatisticalModelType::DeformableMultiObjectType DeformableMultiObjectType;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  AbstractSampler();

  /// Copy constructor.
  AbstractSampler(const AbstractSampler &other);

  /// Makes a copy of the object.
  virtual AbstractSampler *Clone() = 0;

  /// Destructor
  virtual ~AbstractSampler();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns true if the sampler is a symmetric random walk.
  inline bool IsSrwMhwg() const { return (m_Type == SrwMhwg); }
  /// Sets the symmetric random walk type.
  inline void SetSrwMhwgType() { m_Type = SrwMhwg; }

  /// Returns true if the sampler is a mala.
  inline bool IsMala() const { return (m_Type == Mala); }
  /// Sets the amala type.
  inline void SetMalaType() { m_Type = Mala; }

  /// Returns true if the sampler is a mala.
  inline bool IsAmala() const { return (m_Type == Amala); }
  /// Sets the amala type.
  inline void SetAmalaType() { m_Type = Amala; }

  /// Returns the statistical model.
  std::shared_ptr<StatisticalModelType> GetStatisticalModel() const { return m_StatisticalModel; }
  /// Sets the statistical model to \e model.
  void SetStatisticalModel(std::shared_ptr<StatisticalModelType> const model) { m_StatisticalModel = model; }

  /// Returns the data set.
//  StatisticalModelType *GetDataSet() const { return m_DataSet; }
  /// Sets the data set to \e dataSet
  void SetDataSet(LongitudinalDataSetType *const dataSet) {
    m_DataSet = dataSet;
    m_NumberOfSubjects = dataSet->GetNumberOfSubjects();
  }

  /// Sets the target for the acceptances rates.
  void SetAcceptanceRatesTarget(ScalarType const &t) { m_AcceptanceRatesTarget = t; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples. ("RER" Random Effects Realization).
  virtual void Sample(LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER,
                      VectorType &acceptanceRates) = 0;
  /// Variation with a temperature parameter.
  virtual void Sample(LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER,
                      VectorType &acceptanceRates,
                      const ScalarType &temperature) {
    std::cout << "Warning : the tempered version of AbstractSampler::Sample is not available yet." << std::endl;
    Sample(popRER, indRER, acceptanceRates);
  }

  /// Initializes the sampler.
  virtual void Initialize(LinearVariableMapType const &popRER,
                          LinearVariablesMapType const &indRER) = 0;

  /// Adapts the proposal distributions based on the detected acceptance rates and the iteration number.
  virtual void AdaptProposalDistributions(const VectorType &detectedAcceptanceRates,
                                          const unsigned int &iterationNumber,
                                          const bool verbose = false) = 0;

  /// Saves the current state of the sampler.
  virtual void SaveState(def::utils::DeformationState &state) const {
    std::cout << "Warning : the SaveState() method is not available yet for the used sampler. Ignoring."
              << std::endl;
  }
  /// Recovers a previously saved state.
  virtual void RecoverState(def::utils::DeformationState &state) {
    std::cout << "Warning : the RecoverState() method is not available yet for the used sampler. Ignoring."
              << std::endl;
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Type of the sampler.
  SamplerType m_Type;

  /// Target statistical model.
  std::shared_ptr<StatisticalModelType> m_StatisticalModel;

  /// Target data set.
  LongitudinalDataSetType *m_DataSet;
  /// Number of subjects in the dataset.
  unsigned long m_NumberOfSubjects;

  /// Targeted acceptance rate.
  ScalarType m_AcceptanceRatesTarget;

}; /* class AbstractSampler */
