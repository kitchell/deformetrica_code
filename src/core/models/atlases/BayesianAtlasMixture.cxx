/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "BayesianAtlasMixture.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class ScalarType, unsigned int Dimension>
BayesianAtlasMixture<ScalarType, Dimension>
::BayesianAtlasMixture() : Superclass(), m_NumberOfAtlases(1) {
  this->SetBayesianAtlasMixtureType();

  /// atlas coefficients random effect.
  std::shared_ptr<DirichletDistributionType> ccRE = std::make_shared<DirichletDistributionType>();
  Superclass::m_IndividualRandomEffects["MixtureCoefficients"] = ccRE;

  /// Prior on the atlas coefficients.
  std::shared_ptr<NormalDistributionType> ccPrior = std::make_shared<NormalDistributionType>();
  Superclass::m_Priors["MixtureCoefficients"] = ccPrior;
}

template<class ScalarType, unsigned int Dimension>
BayesianAtlasMixture<ScalarType, Dimension>
::~BayesianAtlasMixture() {}

template<class ScalarType, unsigned int Dimension>
BayesianAtlasMixture<ScalarType, Dimension>
::BayesianAtlasMixture(const BayesianAtlasMixture &other) : Superclass(other) {
  m_NumberOfAtlases = other.m_NumberOfAtlases;
  m_Atlases.resize(m_NumberOfAtlases);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    m_Atlases[k] = other.m_Atlases[k]->Clone();
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::SetFixedEffects(const LinearVariableMapType &map) {
  /// Sets the mixture fixed effects.
  Superclass::SetFixedEffects(map);

  /// Updates the fixed effects of the member atlases, and gets the impacted random ones.
  std::vector<LinearVariableMapType> atlasFixedEffects(m_NumberOfAtlases);
  std::vector<ProbabilityDistributionMapType> atlasPopRE(m_NumberOfAtlases);
  std::vector<ProbabilityDistributionMapType> atlasIndRE(m_NumberOfAtlases);

  ExtractAtlasInformation(map, atlasFixedEffects);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    m_Atlases[k]->SetFixedEffects(atlasFixedEffects[k]);
    m_Atlases[k]->GetPopulationRandomEffects(atlasPopRE[k]);
    m_Atlases[k]->GetIndividualRandomEffects(atlasIndRE[k]);
  }

  /// Updates the mixture random effects.
  ProbabilityDistributionMapType popRE;
  ProbabilityDistributionMapType indRE;
  MergeAtlasInformation(atlasPopRE, popRE);
  MergeAtlasInformation(atlasIndRE, indRE);

  Superclass::SetPopulationRandomEffects(popRE);
  Superclass::SetIndividualRandomEffects(indRE);

  /// Mixture-specific fixed effects.
  SetMixtureCoefficientsRandomEffectConcentrationParameters(map.at("ConcentrationParameters"));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::Update() {
  /// Gathers information about the atlases.
  std::vector<LinearVariableMapType> atlasFixedEffects(m_NumberOfAtlases);
  std::vector<ProbabilityDistributionMapType> atlasPriors(m_NumberOfAtlases);
  std::vector<ProbabilityDistributionMapType> atlasPopRE(m_NumberOfAtlases);
  std::vector<ProbabilityDistributionMapType> atlasIndRE(m_NumberOfAtlases);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    std::cout << "-------- Initialization of atlas " << k << " --------" << std::endl;
    m_Atlases[k]->Update();
    m_Atlases[k]->GetFixedEffects(atlasFixedEffects[k]);
    m_Atlases[k]->GetPriors(atlasPriors[k]);
    m_Atlases[k]->GetPopulationRandomEffects(atlasPopRE[k]);
    m_Atlases[k]->GetIndividualRandomEffects(atlasIndRE[k]);
  }
  std::cout << "------ End of atlases initialization. ------" << std::endl;

  /// Merges this infomation.
  LinearVariableMapType fixedEffects;
  ProbabilityDistributionMapType priors;
  ProbabilityDistributionMapType popRE;
  ProbabilityDistributionMapType indRE;

  MergeAtlasInformation(atlasFixedEffects, fixedEffects);
  MergeAtlasInformation(atlasPriors, priors);
  MergeAtlasInformation(atlasPopRE, popRE);
  MergeAtlasInformation(atlasIndRE, indRE);

  Superclass::SetFixedEffects(fixedEffects);
  Superclass::SetPriors(priors);
  Superclass::SetPopulationRandomEffects(popRE);
  Superclass::SetIndividualRandomEffects(indRE);

  /// Mixture-specific initialization.
  InitializeMixtureCoefficientsRandomEffect();
  InitializeMixtureCoefficientsPrior();
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::Write(const LongitudinalDataSetType *const dataSet,
        LinearVariableMapType const &popRER,
        LinearVariablesMapType const &indRER) const {
  const std::string outputDir = def::utils::settings.output_dir;

  /// Write atlases information.
  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
    m_Atlases[k]->Write(dataSet, atlasPopRERs[k], atlasIndRERs[k]);

  /// Writes the atlas coefficients.
  const std::vector<VectorType> mixtureCoefficients = recast<VectorType>(indRER.at("MixtureCoefficients"));
  const unsigned int nbSubjects = indRER.at(indRER.begin()->first).size();

  MatrixType cc(nbSubjects, m_NumberOfAtlases);
  for (int i = 0; i < nbSubjects; i++) {
    for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
      cc(i, k) = mixtureCoefficients[i][k];
  }

  std::ostringstream oss;
  oss << outputDir << this->m_Name << "_MixtureCoefficients.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss.str().c_str(), cc);
}

template<class ScalarType, unsigned int Dimension>
bool
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<std::vector<ScalarType>>> &residuals) {
  std::vector<std::vector<ScalarType>> aux;

  bool oob = ComputeResiduals(dataSet, popRER, indRER, aux);
  for (unsigned int i = 0; i < aux.size(); ++i) {
    residuals[i].resize(1);
    residuals[i][0] = aux[i];
  }

  return oob;
}

template<class ScalarType, unsigned int Dimension>
bool
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<ScalarType>> &residuals) {
  const std::vector<VectorType> mixtureCoefficients = recast<VectorType>(indRER.at("MixtureCoefficients"));

  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  /// Computes the residuals for each atlas.
  std::vector<std::vector<std::vector<ScalarType>>> atlasResiduals(m_NumberOfAtlases);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
    if (m_Atlases[k]->ComputeResiduals(dataSet, atlasPopRERs[k], atlasIndRERs[k], atlasResiduals[k]))
      return true;

  /// Initializes the output container.
  const unsigned int nbSubjects = atlasResiduals[0].size();
  const unsigned int nbObjects = atlasResiduals[0][0].size();
  residuals.resize(nbSubjects);
  for (unsigned int i = 0; i < nbSubjects; ++i) {
    residuals[i].resize(nbObjects);
    std::fill(residuals[i].begin(), residuals[i].end(), 0.0);
  }

  /// Averages over all the atlas models.
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    for (unsigned int i = 0; i < nbSubjects; ++i) {
      for (unsigned int j = 0; j < nbObjects; ++j) {
        residuals[i][j] += mixtureCoefficients[i][k] * atlasResiduals[k][i][j];
      }
    }
  }

  return false;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                               const LinearVariableMapType &popRER,
                               const LinearVariablesMapType &indRER) {
  const std::vector<VectorType> mixtureCoefficients = recast<VectorType>(indRER.at("MixtureCoefficients"));
  const unsigned int nbSubjects = mixtureCoefficients.size();

  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  ScalarType out = 0.0;

//    /// Log-likelihood mixture.
//    for (unsigned int k = 0 ; k < m_NumberOfAtlases ; ++k)
//    {
//        ProbabilityDistributionMapType atlasFEs = m_Atlases[k]->GetFixedEffects();
//        ProbabilityDistributionMapType atlasPopREDs = m_Atlases[k]->GetPopulationRandomEffects();
//        ProbabilityDistributionMapType atlasIndREDs = m_Atlases[k]->GetIndividualRandomEffects();
//
//        ResidualsType residuals;
//        m_Atlases[k]->ComputeResiduals(dataSet, atlasPopRERs, atlasIndRERs, residuals);
//
//        const unsigned int nbObjects = residuals[0].size();
//        const VectorType noiseVariances = recast<VectorType>(atlasFEs.at("NoiseVariance"));
//        const std::vector<MatrixType> momentas = recast<MatrixType>(atlasIndRERs.at("Momenta"));
//
//        for (unsigned int i = 0 ; i < nbSubjects ; ++i)
//        {
//            /// Data (residuals) term.
//            out -= mixtureCoefficients[i][k] * 0.5 *
//                   * m_Atlases[k]->ComputeModelLogLikelihoodForSubject(dataSet, atlasPopRERs[k], atlasIndRERs[k]);
//
//            /// Regularity (momenta random effect) term.
//            out += mixtureCoefficients[i][k] * atlasIndREDs["Momenta"]->ComputeLogLikelihood(momentas[i]);
//        }
//
//        /// Prior on the noise variance.
//        const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
//        out += this->m_Priors["NoiseVariance"]->ComputeLogLikelihood(noiseVariances);
//
//        /// Prior on the momenta covariance matrix.
//        out += this->m_Priors["CovarianceMomenta"]
//            ->ComputeLogLikelihood(this->m_FixedEffects["CovarianceMomentaInverse"]);
//
//        /// Optional prior on the photometric weights.
//        if (m_UseParametricTemplate)
//            out += this->m_Priors["ParametricTemplate"]
//                ->ComputeLogLikelihood(this->GetTemplateData()[this->m_Template->GetImageIndex()]);
//
//        /// Optional prior on the control points.
//        if (m_UsePriorOnControlPoints)
//            out += this->m_Priors["ControlPoints"]->ComputeLogLikelihood(this->GetControlPoints());
//
//        /// Optional term on the control points random effect.
//        if (m_UseRandomControlPoints)
//            out += this->m_PopulationRandomEffects["ControlPoints"]->ComputeLogLikelihood(popRER.at("ControlPoints"));
//    }


  return 0.0;
}

template<class ScalarType, unsigned int Dimension>
bool
BayesianAtlasMixture<ScalarType, Dimension>
::UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    VectorType &logLikelihoodTerms) {
  /// Quick and dirty.
  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  bool oob = false;
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    oob = oob or m_Atlases[k]->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
        dataSet, atlasPopRERs[k], atlasIndRERs[k], logLikelihoodTerms);
  }

  return oob;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                            const LinearVariableMapType &popRER,
                            const LinearVariablesMapType &indRER,
                            std::vector<ScalarType> &contributions) {
  const std::vector<VectorType> mixtureCoefficients = recast<VectorType>(indRER.at("MixtureCoefficients"));
  const unsigned int nbSubjects = indRER.at(indRER.begin()->first).size();

  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  /// Computes the model log-likelihood for each atlas.
  std::vector<std::vector<ScalarType>> atlasModelLogLikelihoodForSubject(m_NumberOfAtlases);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    atlasModelLogLikelihoodForSubject[k].resize(nbSubjects);
    m_Atlases[k]->ComputeModelLogLikelihood(dataSet, atlasPopRERs[k], atlasIndRERs[k],
                                            atlasModelLogLikelihoodForSubject[k]);
  }

  /// Averages over all the atlas models.
  contributions.resize(nbSubjects);
  for (unsigned int i = 0; i < nbSubjects; ++i) {
    contributions[i] = 0.0;
    for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
      contributions[i] += mixtureCoefficients[i][k] * atlasModelLogLikelihoodForSubject[k][i];
  }

  ScalarType out = std::accumulate(contributions.begin(), contributions.end(), 0.0);
  return out;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                      const LinearVariableMapType &popRER,
                                      const LinearVariablesMapType &indRER,
                                      const unsigned int &i) {
  const std::vector<VectorType> mixtureCoefficients = recast<VectorType>(indRER.at("MixtureCoefficients"));

  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  /// Computes the model log-likelihood for each atlas.
  std::vector<ScalarType> atlasModelLogLikelihoodForSubject(m_NumberOfAtlases);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
    atlasModelLogLikelihoodForSubject[k]
        = m_Atlases[k]->ComputeModelLogLikelihoodForSubject(dataSet, atlasPopRERs[k], atlasIndRERs[k], i);

  /// Averages over all the atlas models.
  ScalarType out = 0.0;
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
    out += mixtureCoefficients[i][k] * atlasModelLogLikelihoodForSubject[k];

  return out;
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad) {

}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad,
                                       VectorType &gradSquaredNorms) {

}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeLogLikelihoodGradientWrtRandomEffects(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               LinearVariableMapType &popGrad,
                                               LinearVariablesMapType &indGrad) {

}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                              const LinearVariableMapType &popRER,
                              const LinearVariablesMapType &indRER,
                              LinearVariableMapType &sufficientStatistics) {
  /// Common initializations.
  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> targets = crossSectionalDataSet->GetDeformableMultiObjects();

  const std::vector<VectorType> mixtureCoefficients = recast<VectorType>(indRER.at("MixtureCoefficients"));
  const unsigned int nbSubjects = mixtureCoefficients.size();
  const unsigned int imgIndex = targets[0]->GetImageIndex();

  /// Computes the sufficient statistics for each atlas.
  std::vector<LinearVariableMapType> atlasPopRERs;
  std::vector<LinearVariablesMapType> atlasIndRERs;
  ExtractAtlasInformation(popRER, atlasPopRERs);
  ExtractAtlasInformation(indRER, atlasIndRERs);

  std::vector<LinearVariableMapType> atlasSufficientStatistics(m_NumberOfAtlases);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    std::vector<MatrixType> momentas = recast<MatrixType>(atlasIndRERs[k].at("Momenta"));

    // Initialization before the loop.
    m_Atlases[k]->ComputeSufficientStatisticsSubject(momentas[0], targets[0], atlasSufficientStatistics[k]);
    atlasSufficientStatistics[k]["S2"] *= mixtureCoefficients[0][k];
    atlasSufficientStatistics[k]["S3"] *= mixtureCoefficients[0][k];
    atlasSufficientStatistics[k]["S4"] *= mixtureCoefficients[0][k];

    LinearVariableMapType aux;
    ScalarType S0 = log(mixtureCoefficients[0][k]);
    for (unsigned int i = 1; i < nbSubjects; ++i) {
      m_Atlases[k]->ComputeSufficientStatisticsSubject(momentas[i], targets[i], aux);
      aux["S2"] *= mixtureCoefficients[i][k];
      aux["S3"] *= mixtureCoefficients[i][k];
      aux["S4"] *= mixtureCoefficients[i][k];
      atlasSufficientStatistics[k] += aux;
      S0 += log(mixtureCoefficients[i][k]);
    }

    //atlasSufficientStatistics[k] *= nbSubjects / S0;
    atlasSufficientStatistics[k]["S0"] = S0;
    atlasSufficientStatistics[k]["S1"] = recast<MatrixType>(atlasPopRERs[k]["ControlPoints"]);
  }

  /// Merges the information.
  MergeAtlasInformation(atlasSufficientStatistics, sufficientStatistics);
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                     const LinearVariableMapType &sufficientStatistics) {
  /// Updates the member atlases fixed effects, and gets the updated random ones.
  std::vector<LinearVariableMapType> atlasSufficientStatistics;
  std::vector<ProbabilityDistributionMapType> atlasPopRE(m_NumberOfAtlases);
  std::vector<ProbabilityDistributionMapType> atlasIndRE(m_NumberOfAtlases);

  ExtractAtlasInformation(sufficientStatistics, atlasSufficientStatistics);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    m_Atlases[k]->UpdateFixedEffects(dataSet, atlasSufficientStatistics[k]);
    m_Atlases[k]->GetPopulationRandomEffects(atlasPopRE[k]);
    m_Atlases[k]->GetIndividualRandomEffects(atlasIndRE[k]);
  }

  /// Updates the mixture random effects.
  ProbabilityDistributionMapType popRE;
  ProbabilityDistributionMapType indRE;
  MergeAtlasInformation(atlasPopRE, popRE);
  MergeAtlasInformation(atlasIndRE, indRE);

  Superclass::SetPopulationRandomEffects(popRE);
  Superclass::SetIndividualRandomEffects(indRE);

  /// Update of the mixture-specific fixed effects.
  VectorType S0(m_NumberOfAtlases, 0.0);
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k)
    S0(k) = recast<ScalarType>(atlasSufficientStatistics[k]["S0"]);
  VectorType a_bar = GetMixtureCoefficientsPriorMean();
  ScalarType sigmaSquared_a = GetMixtureCoefficientsPriorCovariance()(0, 0);

  VectorType aux = a_bar + sigmaSquared_a * S0;
  SetMixtureCoefficientsRandomEffectConcentrationParameters(aux);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

//template<class ScalarType, unsigned int Dimension>
//bool
//BayesianAtlasMixture<ScalarType, Dimension>
//::Sample(DataSetType* realizations,
//         LinearVariableMapType  const& popRER,
//         LinearVariablesMapType const& indRER) const
//{
//    const VectorType mixtureCoefficients = recast<VectorType>(popRER["MixtureCoefficients"]);
//    const unsigned int nbSubjects = indRER.at(indRER.begin()->first).size();
//
//    std::vector<LinearVariableMapType>  atlasPopRERs;
//    std::vector<LinearVariablesMapType> atlasIndRERs;
//    ExtractAtlasInformation(popRER, atlasPopRERs);
//    ExtractAtlasInformation(indRER, atlasIndRERs);
//
//    /// For now, it is assumed that we are dealing with single parametric images as templates.
//    std::vector<MatrixType> atlasImagePoints(nbSubjects);
//
//    /// Initialization with the first atlas.
//    DataSetType* realizations_0 = new DataSetType();
//    m_Atlases[0]->Sample(realizations_0, atlasPopRERs[0], atlasIndRERs[0], IndicesType());
//
//    std::vector<std::shared_ptr<DeformableMultiObjectType>> data = realizations_0->GetDeformableMultiObjects();
//    for (unsigned int s = 0 ; s < nbSubjects ; ++s)
//        atlasImagePoints[s]
//            = mixtureCoefficients[0]
//              * data[s]->GetImageIntensityAndLandmarkPointCoordinates()[data[s]->GetImageIndex()];
//
//    delete realizations_0;
//
//    /// Summing over the other atlases.
//    std::vector<DataSetType*> atlasRealizations;
//    for (unsigned int k = 0 ; k < m_NumberOfAtlases ; ++k)
//    {
//        DataSetType* realizations_k = new DataSetType();
//        m_Atlases[k]->Sample(realizations_k, atlasPopRERs[k], atlasIndRERs[k], IndicesType());
//
//        std::vector<std::shared_ptr<DeformableMultiObjectType>> data = realizations_k->GetDeformableMultiObjects();
//        for (unsigned int s = 0 ; s < nbSubjects ; ++s)
//            atlasImagePoints[s]
//                += mixtureCoefficients[k]
//                   * data[s]->GetImageIntensityAndLandmarkPointCoordinates()[data[s]->GetImageIndex()];
//
//        delete realizations_k;
//    }
//}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::InitializeMixtureCoefficientsRandomEffect() {
  /// Concentration parameters.
  SetConcentrationParameters(VectorType(m_NumberOfAtlases, 1.0));
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::InitializeMixtureCoefficientsPrior() {
  /// Mean.
  SetMixtureCoefficientsPriorMean(VectorType(m_NumberOfAtlases, 1.0));

  /// Sqrt of the covariance matrix.
  SetMixtureCoefficientsPriorCovarianceSqrt(diagonal_matrix<ScalarType>(m_NumberOfAtlases, 0.5));
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ExtractAtlasInformation(LinearVariableMapType const &map,
                          std::vector<LinearVariableMapType> &out) const {
  out.resize(m_NumberOfAtlases);
  for (auto it = map.begin(); it != map.end(); ++it) {
    std::size_t pos = it->first.find("_Atlas");
    if (pos != it->first.npos)
      out[std::stoi(it->first.substr(pos + 6, pos + 7))][it->first.substr(0, pos)] = it->second;
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ExtractAtlasInformation(LinearVariablesMapType const &map,
                          std::vector<LinearVariablesMapType> &out) const {
  out.resize(m_NumberOfAtlases);
  for (auto it = map.begin(); it != map.end(); ++it) {
    std::size_t pos = it->first.find("_Atlas");
    if (pos != it->first.npos)
      out[std::stoi(it->first.substr(pos + 6, pos + 7))][it->first.substr(0, pos)] = it->second;
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::ExtractAtlasInformation(ProbabilityDistributionMapType const &map,
                          std::vector<ProbabilityDistributionMapType> &out) const {
  out.resize(m_NumberOfAtlases);
  for (auto it = map.begin(); it != map.end(); ++it) {
    std::size_t pos = it->first.find("_Atlas");
    if (pos != it->first.npos)
      out[std::stoi(it->first.substr(pos + 6, pos + 7))][it->first.substr(0, pos)] = it->second;
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::MergeAtlasInformation(std::vector<LinearVariableMapType> const &in,
                        LinearVariableMapType &map) const {
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    for (auto it = in[k].begin(); it != in[k].end(); ++it)
      map[it->first + "_Atlas" + std::to_string(k)] = it->second;
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::MergeAtlasInformation(std::vector<LinearVariablesMapType> const &in,
                        LinearVariablesMapType &map) const {
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    for (auto it = in[k].begin(); it != in[k].end(); ++it)
      map[it->first + "_Atlas" + std::to_string(k)] = it->second;
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlasMixture<ScalarType, Dimension>
::MergeAtlasInformation(std::vector<ProbabilityDistributionMapType> const &in,
                        ProbabilityDistributionMapType &map) const {
  for (unsigned int k = 0; k < m_NumberOfAtlases; ++k) {
    for (auto it = in[k].begin(); it != in[k].end(); ++it)
      map[it->first + "_Atlas" + std::to_string(k)] = it->second;
  }
}

template
class BayesianAtlasMixture<double, 2>;
template
class BayesianAtlasMixture<double, 3>;
