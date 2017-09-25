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
#include "AbstractStatisticalModel.h"
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"

/// Core files.
#include "BayesianAtlas.h"

/// Support files.
#include "DirichletDistribution.h"

using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      BayesianAtlasMixture object class.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    An atlas is the combination of a template shape, control points, and parameters of variability
 *	            such as covariance matrices of the momentum vectors.
 */
template<class ScalarType, unsigned int Dimension>
class BayesianAtlasMixture : public AbstractStatisticalModel<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Statistical model type.
  typedef AbstractStatisticalModel<ScalarType, Dimension> Superclass;

  /// Longitudinal dataset type.
  typedef typename Superclass::LongitudinalDataSetType LongitudinalDataSetType;
  /// Cross-sectional data set type.
  typedef CrossSectionalDataSet<ScalarType, Dimension> CrossSectionalDataSetType;

  /// Abstract atlas type.
  typedef BayesianAtlas<ScalarType, Dimension> BayesianAtlasType;
  /// Multi-deformable object type.
  typedef typename BayesianAtlasType::DeformableMultiObjectType DeformableMultiObjectType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  BayesianAtlasMixture();

  /// Copy constructor.
  BayesianAtlasMixture(const BayesianAtlasMixture &other);

  /// Makes a copy of the object.
  std::shared_ptr<BayesianAtlasMixture> Clone() const {
    return std::static_pointer_cast<BayesianAtlasMixture>(doClone());
  }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<Superclass> doClone() const {
    return std::static_pointer_cast<Superclass>(std::make_shared<BayesianAtlasMixture>(*this));
  };
 public:

  /// Destructor.
  virtual ~BayesianAtlasMixture();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the fixed effects.
  virtual void SetFixedEffects(const LinearVariableMapType &map);

  /// Gets the member atlases.
  std::vector<std::shared_ptr<BayesianAtlasType>> GetAtlases() const { return m_Atlases; }
  /// Sets the member atlases.
  void SetAtlases(const std::vector<std::shared_ptr<BayesianAtlasType>> atlases) {
    m_Atlases = atlases;
    m_NumberOfAtlases = atlases.size();
  }

  /// Gets the atlas coefficients concentration paramters fixed effect.
  VectorType GetConcentrationParameters() const {
    return this->m_FixedEffects.at("ConcentrationParameters").vectorize();
  }
  /// Sets the atlas coefficients concentration paramters fixed effect.
  void SetConcentrationParameters(const VectorType &vec) {
    this->m_FixedEffects["ConcentrationParameters"] = vec;
    SetMixtureCoefficientsRandomEffectConcentrationParameters(vec);
  }

  /// Gets the concentration parameters of the atlas coefficients random effect.
  VectorType GetMixtureCoefficientsRandomEffectConcentrationParameters() const {
    return std::static_pointer_cast<DirichletDistributionType>(
        this->m_IndividualRandomEffects.at("MixtureCoefficients"))->GetConcentrationParameters();
  }
  /// Sets the concentration parameters of the atlas coefficients random effect.
  void SetMixtureCoefficientsRandomEffectConcentrationParameters(const VectorType &cp) {
    return std::static_pointer_cast<DirichletDistributionType>(
        this->m_IndividualRandomEffects.at("MixtureCoefficients"))->SetConcentrationParameters(cp);
  }
  /// Sets the concentration parameters of the atlas coefficients random effect (input as a linear variable).
  void SetMixtureCoefficientsRandomEffectConcentrationParameters(const LinearVariableType &cp) {
    return std::static_pointer_cast<DirichletDistributionType>(
        this->m_IndividualRandomEffects.at("MixtureCoefficients"))->SetConcentrationParameters(cp);
  }

  /// Gets the mean of the atlas coefficients prior.
  VectorType GetMixtureCoefficientsPriorMean() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("MixtureCoefficients"))->GetMean();
  }
  /// Sets the mean of the atlas coefficients prior.
  void SetMixtureCoefficientsPriorMean(const VectorType &mean) {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("MixtureCoefficients"))->SetMean(mean);
  }

  /// Gets the covariance matrix of the atlas coefficients prior.
  MatrixType GetMixtureCoefficientsPriorCovariance() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("MixtureCoefficients"))->GetCovariance();
  }
  /// Gets the sqrt of the covariance matrix of the atlas coefficients prior.
  MatrixType GetMixtureCoefficientsPriorCovarianceSqrt() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("MixtureCoefficients"))->GetCovarianceSqrt();
  }
  /// Sets the sqrt of the covariance matrix of the atlas coefficients prior.
  void SetMixtureCoefficientsPriorCovarianceSqrt(const MatrixType &covSqrt) {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("MixtureCoefficients"))->SetCovarianceSqrt(covSqrt);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates the bounding box and initializes the control points if needed.
  virtual void Update();

  /// Saves the model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const;

  /// Applies the atlas mixture model (at specified time points), given specific random \e effects.
  virtual bool Sample(LongitudinalDataSetType *realizations,
                      LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER) const {
    std::cout << "Warning : BayesianAtlasMixture::Sample method is not available yet." << std::endl;
    return true;
  }

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                LinearVariableMapType const &popRER,
                                LinearVariablesMapType const &indRER,
                                std::vector<std::vector<std::vector<ScalarType>>> &residuals);
  /// Variation with a residuals structure adapted to cross-sectional data.
  bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                        LinearVariableMapType const &popRER,
                        LinearVariablesMapType const &indRER,
                        std::vector<std::vector<ScalarType>> &residuals);

  /// Computes the complete log-likelihood, given an input random effects realization ("RER").
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER);

  /// Updates the fixed effects and computes the complete likelihood at once.
  virtual bool UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                                 const LinearVariableMapType &popRER,
                                                                 const LinearVariablesMapType &indRER,
                                                                 VectorType &logLikelihoodTerms);

  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions);

  /// Computes the model log-likelihood, ignoring the contribution of subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i);

  /// Computes the gradient of the log-likelihood mode. ("RER" Random Effects Realization).
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad);
  /// Computes the gradient of the log-likelihood mode and the corresponding squared norms.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad,
                                                    VectorType &gradSquaredNorms);
  /// Computes the marginal log-likelihood gradient with respect to the random effects of the model.
  virtual void ComputeLogLikelihoodGradientWrtRandomEffects(const LongitudinalDataSetType *const dataSet,
                                                            const LinearVariableMapType &popRER,
                                                            const LinearVariablesMapType &indRER,
                                                            LinearVariableMapType &popGrad,
                                                            LinearVariablesMapType &indGrad);

  /// Computes the sufficient statistics of the atlas model.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics);

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initializes the atlas coefficients random effect.
  void InitializeMixtureCoefficientsRandomEffect();
  /// Initializes the prior on the atlas coefficients.
  void InitializeMixtureCoefficientsPrior();

  /// Extracts the map components relevant to the atlases of the mixture.
  void ExtractAtlasInformation(LinearVariableMapType const &map,
                               std::vector<LinearVariableMapType> &out) const;
  void ExtractAtlasInformation(LinearVariablesMapType const &map,
                               std::vector<LinearVariablesMapType> &out) const;
  void ExtractAtlasInformation(ProbabilityDistributionMapType const &map,
                               std::vector<ProbabilityDistributionMapType> &out) const;

  /// Merges the map components relevant to the atlases of the mixture.
  void MergeAtlasInformation(std::vector<LinearVariableMapType> const &in,
                             LinearVariableMapType &map) const;
  void MergeAtlasInformation(std::vector<LinearVariablesMapType> const &in,
                             LinearVariablesMapType &map) const;
  void MergeAtlasInformation(std::vector<ProbabilityDistributionMapType> const &in,
                             ProbabilityDistributionMapType &map) const;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Atlas models.
  std::vector<std::shared_ptr<BayesianAtlasType>> m_Atlases;

  /// Number of atlases.
  unsigned int m_NumberOfAtlases;

}; /* class BayesianAtlasMixture */


