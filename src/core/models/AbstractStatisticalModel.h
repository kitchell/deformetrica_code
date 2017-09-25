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
#include "src/core/observations/deformable_objects/DeformableMultiObject.h"
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"
#include "LongitudinalDataSet.h"

/// Support files.
#include <src/support/utilities/GeneralSettings.h>

/// Librairies files.
#include <lib/ThreadPool/ThreadPool.h>
#include <vector>
#include <map>
#include <iostream>
#include <cstring>
#include <sstream>
#include <memory>

using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      AbstractStatisticalModel object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A statistical model is a generative function, which tries to explain an observed stochastic process.
 */
template<class ScalarType, unsigned int Dimension>
class AbstractStatisticalModel {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Possible types of statistical model.
  typedef enum {
    null,                   /*!< Null value. */
    DeterministicAtlas,     /*!< Simple deterministic atlas model. (See DeterministicAtlas). */
    BayesianAtlas,          /*!< Bayesian atlas model. (See BayesianAtlas). */
    BayesianAtlasMixture,   /*!< Bayesian atlas mixture model. (See BayesianAtlasMixture). */
    Regression,             /*!< Regression model. (See Regression). */
    LongitudinalAtlas,      /*!< Longitudinal model. (See LongitudinalAtlas). */
    LongitudinalMatching,   /*!< Longitudinal matching. (See LongitudinalMatching). */
    LdaAtlas                /*!< Linear discriminant analysis atlas. (See LdaAtlas). */
    // EBM, ...
  } StatisticalModelType;

  /// Longitudinal dataset type.
  typedef LongitudinalDataSet<ScalarType, Dimension> LongitudinalDataSetType;

  /// Multi-deformable object type.
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  AbstractStatisticalModel();

  /// Copy constructor.
  AbstractStatisticalModel(const AbstractStatisticalModel &other);

  /// Makes a copy of the object.
  std::shared_ptr<AbstractStatisticalModel> Clone() const { return doClone(); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<AbstractStatisticalModel> doClone() const = 0;
 public:

  /// Destructor.
  virtual ~AbstractStatisticalModel() = 0;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the model name.
  std::string GetName() const { return m_Name; }
  /// Sets the model name.
  void SetName(const std::string name) { m_Name = name; }


  /// Returns true if the statistical model is a deterministic atlas one.
  bool IsDeterministicAtlas() const { return (m_Type == DeterministicAtlas); }
  /// Sets deterministic atlas type.
  void SetDeterministicAtlasType() { m_Type = DeterministicAtlas; }

  /// Returns true if the statistical model is a bayesian atlas one.
  bool IsBayesianAtlas() const { return (m_Type == BayesianAtlas); }
  /// Sets bayesian atlas type.
  void SetBayesianAtlasType() { m_Type = BayesianAtlas; }

  /// Returns true if the statistical model is a bayesian atlas mixture one.
  bool IsBayesianAtlasMixture() const { return (m_Type == BayesianAtlasMixture); }
  /// Sets bayesian atlas mixture type.
  void SetBayesianAtlasMixtureType() { m_Type = BayesianAtlasMixture; }

  /// Returns true if the statistical model is a regression one.
  bool IsRegression() const { return (m_Type == Regression); }
  /// Sets regression type.
  void SetRegressionType() { m_Type = Regression; }

  /// Returns true if the statistical model is a longitudinal one.
  bool IsLongitudinalAtlas() const { return (m_Type == LongitudinalAtlas); }
  /// Sets longitudinal atlas type.
  void SetLongitudinalAtlasType() { m_Type = LongitudinalAtlas; }

  /// Returns true if the statistical model is a longitudinal matching one.
  bool IsLongitudinalMatching() const { return (m_Type == LongitudinalMatching); }
  /// Sets longitudinal matching type.
  void SetLongitudinalMatchingType() { m_Type = LongitudinalMatching; }

  /// Returns true if the statistical model is a Lda one.
  bool IsLdaAtlas() const { return (m_Type == LdaAtlas); }
  /// Sets longitudinal atlas type.
  void SetLdaAtlasType() { m_Type = LdaAtlas; }

  /// Returns the fixed effects.
  void GetFixedEffects(LinearVariableMapType &map) const { map = m_FixedEffects; }
  /// Sets the fixed effects.
  virtual void SetFixedEffects(const LinearVariableMapType &map) {
    for (auto it = map.begin(); it != map.end(); ++it) { m_FixedEffects[it->first] = it->second; }
  }

  /// Returns the population random effects.
  void GetPriors(ProbabilityDistributionMapType &map) const { map = m_Priors; }
  /// Sets the population random effects.
  void SetPriors(const ProbabilityDistributionMapType &map) {
    for (auto it = map.begin(); it != map.end(); ++it) { m_Priors[it->first] = it->second; }
  }

  /// Returns the population random effects.
  void GetPopulationRandomEffects(ProbabilityDistributionMapType &map) const {
    map = m_PopulationRandomEffects;
  }
  /// Sets the population random effects.
  void SetPopulationRandomEffects(const ProbabilityDistributionMapType &map) {
    for (auto it = map.begin(); it != map.end(); ++it) { m_PopulationRandomEffects[it->first] = it->second; }
  }

  /// Returns the individual random effects.
  void GetIndividualRandomEffects(ProbabilityDistributionMapType &map) const {
    map = m_IndividualRandomEffects;
  }
  /// Sets the individual random effects.
  void SetIndividualRandomEffects(const ProbabilityDistributionMapType &map) {
    for (auto it = map.begin(); it != map.end(); ++it) { m_IndividualRandomEffects[it->first] = it->second; }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates and initializes parameters if needed.
  virtual void Update() = 0;

  /// Prints information about the model.
  virtual void Print() const {}
  /// Saves the model. ("RER" Random Effects Realization).
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const = 0;

  /// Samples from the model to generate independent realizations.
  virtual bool Sample(LongitudinalDataSetType *dataSet) const {
    LinearVariableMapType popRER;
    LinearVariablesMapType indRER;
    return Sample(dataSet, popRER, indRER);
  }
  /// Samples from the model, and stores the used parameters in \e popRER and \e indRER.
  virtual bool Sample(LongitudinalDataSetType *dataSet,
                      LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER) const {
    std::cerr << "Error : the Sample() method is unimplemented for the considered model" << std::endl;
    return false;
  }

  /// Computes the residuals. ("RER" Random Effects Realization).
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                const LinearVariableMapType &popRER,
                                const LinearVariablesMapType &indRER,
                                std::vector<std::vector<std::vector<ScalarType>>> &residuals) = 0;

  /// Computes the complete log-likelihood, given an input random effects realization ("RER").
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER) = 0;
  /// Updates the fixed effects and computes the complete likelihood at once.
  virtual bool UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                                 const LinearVariableMapType &popRER,
                                                                 const LinearVariablesMapType &indRER,
                                                                 VectorType &logLikelihoodTerms) {
    std::cerr << "Error : the UpdateFixedEffectsAndComputeCompleteLogLikelihood() method is "
        "unimplemented for the considered model" << std::endl;
    return false;
  }

  /// Computes the model log-likelihood.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER) {
    std::vector<ScalarType> contributions;
    return ComputeModelLogLikelihood(dataSet, popRER, indRER, contributions);
  }
  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions) {
    std::cerr << "Error : the ComputeModelLogLikelihood() method is unimplemented for the considered model"
              << std::endl;
    return 0.0;
  }
  /// Optional variation, where only one population random effect realization (RER) has been modified since last call.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               const std::string &modifiedVar,
                                               std::vector<ScalarType> &contributions) {
    return ComputeModelLogLikelihood(dataSet, popRER, indRER, contributions);
  };
  /// Optional variation, featuring a temperature parameter.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               const ScalarType &temperature,
                                               const std::string &modifiedVar,
                                               std::vector<ScalarType> &contributions) {
    return ComputeModelLogLikelihood(dataSet, popRER, indRER, modifiedVar, contributions);
  };

  /// Computes the model log-likelihood, ignoring the contribution of subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i) {
    std::cerr << "Error : the ComputeModelLogLikelihoodForSubject() method is unimplemented for the considered model"
              << std::endl;
    return 0.0;
  }
  /// Optional variation, where only one individual random effect realization (RER) has been modified since last call.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i,
                                                         const std::string &modifiedVar) {
    return ComputeModelLogLikelihoodForSubject(dataSet, popRER, indRER, i);
  };

  /// Computes the gradient of the log-likelihood mode. ("RER" Random Effects Realization).
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad) {
    std::cerr << "Error : the ComputeCompleteLogLikelihoodGradient() method is unimplemented for the considered model"
              << std::endl;
  }
  /// Computes the gradient of the log-likelihood mode and the corresponding squared norms.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad,
                                                    VectorType &gradSquaredNorms) {
    std::cerr << "Error : the ComputeCompleteLogLikelihoodGradient() method is unimplemented for the considered model"
              << std::endl;
  }
  /// Computes the marginal log-likelihood gradient with respect to the random effects of the model.
  virtual void ComputeLogLikelihoodGradientWrtRandomEffects(const LongitudinalDataSetType *const dataSet,
                                                            const LinearVariableMapType &popRER,
                                                            const LinearVariablesMapType &indRER,
                                                            LinearVariableMapType &popGrad,
                                                            LinearVariablesMapType &indGrad) {
    std::cerr << "Error : the ComputeLogLikelihoodGradientWrtRandomEffects() method is "
        "unimplemented for the considered model" << std::endl;
  }

  /// Computes the sufficient statistics of the statistical model. ("RER" Random Effects Realization).
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics) {
    std::cerr << "Error : the ComputeSufficientStatistics() method is unimplemented for the considered model"
              << std::endl;
  }

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics) {
    std::cerr << "Error : the UpdateFixedEffects() method is unimplemented for the considered model"
              << std::endl;
  }
  /// Optional variation, featuring a temperature parameter.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics,
                                  const ScalarType &temperature) {
    UpdateFixedEffects(dataSet, sufficientStatistics);
  }

  /// Recovers the memorized random effect realizations-based state.
  virtual void RecoverMemorizedState() {}

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Instances count, for automatic naming
  static int count;

  /// Name of the model.
  std::string m_Name;

  /// Type of model.
  StatisticalModelType m_Type;

  /// Fixed effects.
  LinearVariableMapType m_FixedEffects;
  /// Priors on fixed effects.
  ProbabilityDistributionMapType m_Priors;

  /// Population random effects.
  ProbabilityDistributionMapType m_PopulationRandomEffects;
  /// Individual random effects.
  ProbabilityDistributionMapType m_IndividualRandomEffects;

}; /* class AbstractStatisticalModel */


