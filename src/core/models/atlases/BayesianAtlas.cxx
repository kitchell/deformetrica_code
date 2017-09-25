/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "BayesianAtlas.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
BayesianAtlas<ScalarType, Dimension>
::BayesianAtlas() : Superclass(), m_UseParametricTemplate(false),
                    m_UsePriorOnControlPoints(false), m_UseRandomControlPoints(false) {
  this->SetBayesianAtlasType();

  /// Momenta individual random effect.
  std::shared_ptr<NormalDistributionType> momRE = std::make_shared<NormalDistributionType>();
  this->m_IndividualRandomEffects["Momenta"] = momRE;

  /// Prior on the momenta covariance.
  std::shared_ptr<InverseWishartDistributionType> cmPrior = std::make_shared<InverseWishartDistributionType>();
  this->m_Priors["CovarianceMomenta"] = cmPrior;

  /// Prior on the noise variance.
  std::shared_ptr<MultiScalarInverseWishartDistributionType>
      nvPrior = std::make_shared<MultiScalarInverseWishartDistributionType>();
  this->m_Priors["NoiseVariance"] = nvPrior;

  /// Control points population random effect. (Will be removed at Update() if not needed).
  std::shared_ptr<NormalDistributionType> cpRE = std::make_shared<NormalDistributionType>();
  this->m_PopulationRandomEffects["ControlPoints"] = cpRE;

  /// Prior on the control points fixed effect. (Will be removed at Update() if not needed).
  std::shared_ptr<NormalDistributionType> cpPrior = std::make_shared<NormalDistributionType>();
  this->m_Priors["ControlPoints"] = cpPrior;
}

template<class ScalarType, unsigned int Dimension>
BayesianAtlas<ScalarType, Dimension>
::~BayesianAtlas() {}

template<class ScalarType, unsigned int Dimension>
BayesianAtlas<ScalarType, Dimension>
::BayesianAtlas(const BayesianAtlas &other) : Superclass(other) {
  m_NoiseDimension = other.m_NoiseDimension;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::SetCovarianceMomenta_Prior_Inverse(std::string const &fn) {
  if (strlen(fn.c_str())) {
    const MatrixType covarianceMomenta_Prior_Inverse = readMatrixDLM<ScalarType>(fn.c_str());
    SetCovarianceMomenta_Prior_Inverse(covarianceMomenta_Prior_Inverse);
    std::cout << "Using a prior on the momenta covariance of size " << covarianceMomenta_Prior_Inverse.rows()
              << " x " << covarianceMomenta_Prior_Inverse.cols() << " from file: " << fn << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::Update() {
  Superclass::Update();

  /// Prior on the momenta covariance.
  if (GetCovarianceMomenta_Prior_Inverse().n_elem() == 1) {
    std::cout << "No prior for momenta covariance matrix given : set to kernel matrix with deformation kernel width = "
              << this->m_Def->GetKernelWidth() << std::endl;

    InitializeCovMomPrior();
  }

  /// Momenta individual random effect.
  InitializeMomentaRandomEffect();

  /// If a parametric template is used, a prior is defined on the photometric weights.
  if (this->m_Template->GetNumberOfImageKindObjects() &&
      this->m_Template->GetObjectList()[this->m_Template->GetImageIndex()]->IsParametricImage()) {
    m_UseParametricTemplate = true;

    std::shared_ptr<NormalDistributionType> ptPrior = std::make_shared<NormalDistributionType>();
    this->m_Priors["ParametricTemplate"] = ptPrior;

    InitializeParametricTemplatePrior();
  }

  /// If needed, define a prior on the control points fixed effect.
  if (m_UsePriorOnControlPoints)
    InitializeControlPointsPrior();
  else
    this->m_Priors.erase("ControlPoints");

  /// If needed, define the control points as a random effect.
  if (m_UseRandomControlPoints)
    InitializeControlPointsRandomEffect();
  else
    this->m_PopulationRandomEffects.erase("ControlPoints");
}

template<class ScalarType, unsigned int Dimension>
bool
BayesianAtlas<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<ScalarType>> &residuals) {
  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));

  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();

  if (m_UseRandomControlPoints) {
    const MatrixType cp = recast<MatrixType>(popRER.at("ControlPoints"));
    return Superclass::ComputeResiduals(cp, momentas, target, residuals);
  } else {
    return Superclass::ComputeResiduals(this->GetControlPoints(), momentas, target, residuals);
  }
}

template<class ScalarType, unsigned int Dimension>
ScalarType
BayesianAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                               const LinearVariableMapType &popRER,
                               const LinearVariablesMapType &indRER) {
  ScalarType out = 0.0;

  /// Data (residuals) term.
  std::vector<std::vector<ScalarType>> residuals;
  ComputeResiduals(dataSet, popRER, indRER, residuals);

  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const unsigned int nbSubjects = residuals.size();
  const unsigned int nbObjects = this->m_NumberOfObjects;
  for (unsigned int s = 0; s < nbSubjects; ++s) {
    for (unsigned int i = 0; i < nbObjects; ++i) {
      out -= 0.5 * (residuals[s][i] / noiseVariances(i) + m_NoiseDimension[i] * log(noiseVariances(i)));
    }
  }

  /// Prior on the noise variance.
  /// HERE PROBLEM IF SEVERAL OBJECTS IN TEMPLATE. PROBLEM ALSO WITH THE NEW CONSTRUCTOR OF MSINVERSEWISHART ??
  /// See LongitudinalAtlas::InitializeNoiseVariables().
  out += this->m_Priors["NoiseVariance"]->ComputeLogLikelihood(noiseVariances);

  /// Regularity (momenta random effect) term.
  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));
  assert(nbSubjects == momentas.size());

  for (unsigned int s = 0; s < nbSubjects; ++s)
    out += this->m_IndividualRandomEffects["Momenta"]->ComputeLogLikelihood(momentas[s]);

  /// Prior on the momenta covariance matrix.
  out += this->m_Priors["CovarianceMomenta"]
      ->ComputeLogLikelihood(this->m_FixedEffects["CovarianceMomentaInverse"]);

  /// Optional prior on the photometric weights.
  if (m_UseParametricTemplate)
    out += this->m_Priors["ParametricTemplate"]
        ->ComputeLogLikelihood(this->GetTemplateData()[this->m_Template->GetImageIndex()]);

  /// Optional prior on the control points.
  if (m_UsePriorOnControlPoints)
    out += this->m_Priors["ControlPoints"]->ComputeLogLikelihood(this->GetControlPoints());

  /// Optional term on the control points random effect.
  if (m_UseRandomControlPoints)
    out += this->m_PopulationRandomEffects["ControlPoints"]->ComputeLogLikelihood(popRER.at("ControlPoints"));

  return out;
}

template<class ScalarType, unsigned int Dimension>
bool
BayesianAtlas<ScalarType, Dimension>
::UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    VectorType &logLikelihoodTerms) {
  std::vector<std::vector<ScalarType>> residuals;
  bool oob = this->ComputeResiduals(dataSet, popRER, indRER, residuals);
  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));

  if (momentas.size() != residuals.size()) {
    throw std::runtime_error(
        "Number of Momentas and target multi-object mismatch in BayesianAtlas::ComputeLogLikelihood");
  }

  ScalarType covMomTerm, covMomPriorTerm, dataTerm, dataPriorTerm;
  UpdateCovMomInverseAndComputeTerms(momentas, covMomTerm, covMomPriorTerm);
  UpdateDataSigmaSquaredAndComputeTerms(residuals, dataTerm, dataPriorTerm);

  logLikelihoodTerms.set_size(2);
  logLikelihoodTerms.fill(0.0);

  logLikelihoodTerms[1] += covMomTerm;
  logLikelihoodTerms[1] += covMomPriorTerm;
  logLikelihoodTerms[0] += dataTerm;
  logLikelihoodTerms[1] += dataPriorTerm;

  if (m_UseParametricTemplate)
    logLikelihoodTerms[1] += this->m_Priors["ParametricTemplate"]
        ->ComputeLogLikelihood(this->GetTemplateData()[this->m_Template->GetImageIndex()]);

  if (m_UsePriorOnControlPoints)
    logLikelihoodTerms[1] += this->m_Priors["ControlPoints"]->ComputeLogLikelihood(this->GetControlPoints());

  if (m_UseRandomControlPoints)
    logLikelihoodTerms[1] +=
        this->m_PopulationRandomEffects["ControlPoints"]->ComputeLogLikelihood(popRER.at("ControlPoints"));

  return oob; // Out of box flag.
}

template<class ScalarType, unsigned int Dimension>
ScalarType
BayesianAtlas<ScalarType, Dimension>
::ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                            const LinearVariableMapType &popRER,
                            const LinearVariablesMapType &indRER,
                            std::vector<ScalarType> &contributions) {
  /// Data (residuals) term.
  std::vector<std::vector<ScalarType>> residuals;
  ComputeResiduals(dataSet, popRER, indRER, residuals);

  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const unsigned int nbSubjects = residuals.size();
  const unsigned int nbObjects = this->m_NumberOfObjects;

  contributions.resize(nbSubjects);
  for (unsigned int s = 0; s < nbSubjects; ++s) {
    contributions[s] = 0.0;
    for (unsigned int i = 0; i < nbObjects; ++i)
      contributions[s] -= 0.5 * residuals[s][i] / noiseVariances(i);
  }

  ScalarType out = std::accumulate(contributions.begin(), contributions.end(), 0.0);
  return out;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
BayesianAtlas<ScalarType, Dimension>
::ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                      const LinearVariableMapType &popRER,
                                      const LinearVariablesMapType &indRER,
                                      const unsigned int &i) {
  ScalarType out = 0.0;

  /// Data (residuals) term.
  MatrixType cp;
  if (m_UseRandomControlPoints) { cp = recast<MatrixType>(popRER["ControlPoints"]); }
  else { cp = this->GetControlPoints(); }

  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> targets = crossSectionalDataSet->GetDeformableMultiObjects();

  const MatrixType momenta = recast<MatrixType>(indRER.at("Momenta"))[i];
  std::vector<ScalarType> residuals;
  Superclass::ComputeResidualsSubject(cp, momenta, targets[i], residuals);

  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const unsigned int nbObjects = this->m_NumberOfObjects;
  for (unsigned int k = 0; k < nbObjects; ++k)
    out -= 0.5 * residuals[k] / noiseVariances(k);

  return out;
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad) {
  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();

  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));
  const unsigned long nbSubjects = target.size();

  MatrixType cp;
  if (m_UseRandomControlPoints)
    cp = recast<MatrixType>(popRER.at("ControlPoints"));
  else
    cp = this->GetControlPoints();

  if (momentas.size() != nbSubjects)
    throw std::runtime_error("Number of Momentas and target multi-object mismatch in BayesianAtlas::ComputeLikelihood");

  const MatrixType covarianceMomentaInverse = GetCovarianceMomentaInverse();
  const VectorType dataSigmaSquared = GetDataSigmaSquared();

  MatrixType gradPos;
  MatrixListType gradTempL_L2, gradTempL_Sob;
  std::vector<MatrixType> gradMom(nbSubjects);

  this->UpdateDeformationAndKernelDataDomain(target);

  if (!this->m_FreezeTemplateFlag) {
    this->ComputeDataTermGradient(cp, momentas, dataSigmaSquared, target, gradPos, gradMom, gradTempL_L2);
    if (m_UseParametricTemplate) {
      const VectorType photoWeights = (this->GetTemplateData()[this->m_Template->GetImageIndex()]).get_column(0);
      gradTempL_L2[this->m_Template->GetImageIndex()].get_column(0)
          -= GetParametricTemplatePriorCovarianceInverse() * (photoWeights - GetParametricTemplatePriorMean());
      gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
    }
  } else { this->ComputeDataTermGradient(cp, momentas, dataSigmaSquared, target, gradPos, gradMom); }

  for (unsigned int s = 0; s < momentas.size(); s++) {
    VectorType aux = covarianceMomentaInverse * this->Vectorize(momentas[s]);
    gradMom[s] -= this->VectorToMatrix(aux);
  }

  if (m_UseRandomControlPoints) {
    const VectorType cp = this->Vectorize(recast<MatrixType>(popRER.at("ControlPoints")));
    gradPos -= this->VectorToMatrix(GetControlPointsRandomEffectCovarianceInverse()
                                        * (cp - GetControlPointsRandomEffectMean()));
  } else if (m_UsePriorOnControlPoints) {
    const VectorType cp = this->Vectorize(this->GetControlPoints());
    gradPos -= this->VectorToMatrix(GetControlPointsPriorCovarianceInverse()
                                        * (cp - GetControlPointsPriorMean()));
  }

  if (!this->m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  if (!this->m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }
  indGrad["Momenta"] = gradMom;
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad,
                                       VectorType &gradSquaredNorms) {
  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();

  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));
  const unsigned long nbSubjects = target.size();

  MatrixType cp;
  if (m_UseRandomControlPoints)
    cp = recast<MatrixType>(popRER.at("ControlPoints"));
  else
    cp = this->GetControlPoints();

  if (momentas.size() != nbSubjects)
    throw std::runtime_error("Number of Momentas and target multi-object mismatch in BayesianAtlas::ComputeLikelihood");

  const MatrixType covarianceMomentaInverse = GetCovarianceMomentaInverse();
  const VectorType dataSigmaSquared = GetDataSigmaSquared();

  MatrixType gradPos;
  MatrixListType gradTempL_L2, gradTempL_Sob;
  std::vector<MatrixType> gradMom(nbSubjects);

  this->UpdateDeformationAndKernelDataDomain(target);

  if (!this->m_FreezeTemplateFlag) {
    this->ComputeDataTermGradient(cp, momentas, dataSigmaSquared, target, gradPos, gradMom, gradTempL_L2);
    if (m_UseParametricTemplate) {
      const VectorType photoWeights = (this->GetTemplateData()[this->m_Template->GetImageIndex()]).get_column(0);
      gradTempL_L2[this->m_Template->GetImageIndex()].get_column(0)
          -= GetParametricTemplatePriorCovarianceInverse() * (photoWeights - GetParametricTemplatePriorMean());
    }
    gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
  } else { this->ComputeDataTermGradient(cp, momentas, dataSigmaSquared, target, gradPos, gradMom); }

  for (unsigned int s = 0; s < momentas.size(); s++) {
    VectorType aux = covarianceMomentaInverse * this->Vectorize(momentas[s]);
    gradMom[s] -= this->VectorToMatrix(aux);
  }

  if (m_UseRandomControlPoints) {
    const VectorType cp = this->Vectorize(recast<MatrixType>(popRER.at("ControlPoints")));
    gradPos -= this->VectorToMatrix(GetControlPointsRandomEffectCovarianceInverse()
                                        * (cp - GetControlPointsRandomEffectMean()));
  } else if (m_UsePriorOnControlPoints) {
    const VectorType cp = this->Vectorize(this->GetControlPoints());
    gradPos -= this->VectorToMatrix(GetControlPointsPriorCovarianceInverse()
                                        * (cp - GetControlPointsPriorMean()));
  }

  if (!this->m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  if (!this->m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }
  indGrad["Momenta"] = gradMom;

  // Must be the population effects in alphabetical order first, the individual effects in alphabetical order then.
  gradSquaredNorms.set_size(0);
  if (!this->m_FreezeControlPointsFlag) { gradSquaredNorms.push_back(gradPos.sum_of_squares()); }
  if (!this->m_FreezeTemplateFlag) { gradSquaredNorms.push_back(dot_product(gradTempL_L2, gradTempL_Sob)); }
  gradSquaredNorms.push_back(MatrixListType(gradMom).sum_of_squares());
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::ComputeLogLikelihoodGradientWrtRandomEffects(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               LinearVariableMapType &popGrad,
                                               LinearVariablesMapType &indGrad) {
  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();

  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));
  const unsigned long nbSubjects = target.size();

  MatrixType cp;
  if (m_UseRandomControlPoints)
    cp = recast<MatrixType>(popRER.at("ControlPoints"));
  else
    cp = this->GetControlPoints();

  assert(momentas.size() == nbSubjects);

  const MatrixType covarianceMomentaInverse = GetCovarianceMomentaInverse();
  const VectorType dataSigmaSquared = GetDataSigmaSquared();

  MatrixType gradPos;
  MatrixListType gradTempL_L2;
  std::vector<MatrixType> gradMom(nbSubjects);

  this->UpdateDeformationAndKernelDataDomain(target);
  this->ComputeDataTermGradient(cp, momentas, dataSigmaSquared, target, gradPos, gradMom);

  for (unsigned int s = 0; s < momentas.size(); s++) {
    VectorType aux = covarianceMomentaInverse * this->Vectorize(momentas[s]);
    gradMom[s] -= this->VectorToMatrix(aux);
  }
  indGrad["Momenta"] = gradMom;

  if (m_UseRandomControlPoints) {
    const VectorType cp = this->Vectorize(recast<MatrixType>(popRER.at("ControlPoints")));
    gradPos -= this->VectorToMatrix(GetControlPointsRandomEffectCovarianceInverse()
                                        * (cp - GetControlPointsRandomEffectMean()));
    popGrad["ControlPoints"] = gradPos;
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                              const LinearVariableMapType &popRER,
                              const LinearVariablesMapType &indRER,
                              LinearVariableMapType &sufficientStatistics) {
  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);

  sufficientStatistics.clear();
  this->UpdateDeformationAndKernelDataDomain(crossSectionalDataSet->GetDeformableMultiObjects());

  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));
  ComputeSufficientStatisticsSubject(momentas[0], crossSectionalDataSet->GetDataForSubject(0), sufficientStatistics);

  LinearVariableMapType aux;
  for (unsigned int s = 1; s < momentas.size(); ++s) {
    ComputeSufficientStatisticsSubject(momentas[s], crossSectionalDataSet->GetDataForSubject(s), aux);
    sufficientStatistics += aux;
  }

  //sufficientStatistics["S0"] = ScalarType(crossSectionalDataSet->GetNumberOfSubjects());
  sufficientStatistics["S1"] = recast<MatrixType>(popRER["ControlPoints"]);
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                     const LinearVariableMapType &sufficientStatistics) {
  //const ScalarType S0 = recast<ScalarType>(sufficientStatistics["S0"]);
  const MatrixType S1 = recast<MatrixType>(sufficientStatistics["S1"]);
  const ScalarType S2 = recast<ScalarType>(sufficientStatistics["S2"]);
  const MatrixType S3 = recast<MatrixType>(sufficientStatistics["S3"]);
  const MatrixType S4 = recast<MatrixType>(sufficientStatistics["S4"]);
  const MatrixType S5 = recast<MatrixType>(sufficientStatistics["S5"]);

  /// Control points random effect mean close-form update.
  SetControlPoints(sufficientStatistics["S1"]);

  /// Covariance of the momenta closed-form update.
  const unsigned int nbSubjects = dataSet->GetNumberOfSubjects();
  const ScalarType covMomHP = GetCovarianceMomenta_HyperParameter();
  const MatrixType covMomPrior = GetCovarianceMomenta_Prior();

  SetCovarianceMomentaInverse((nbSubjects + covMomHP)
                                  * inverse_sympd(S5 + covMomHP * covMomPrior));

  /// Intricated update of the photometric weights w and data noise variance ss.
  const unsigned int imgIndex = this->m_Template->GetImageIndex();
  const unsigned long nbVoxels = m_NoiseDimension[imgIndex];
  const MatrixType paramTempMean = MatrixType(GetParametricTemplatePriorMean());
  const MatrixType paramTempCovInv = GetParametricTemplatePriorCovarianceInverse();
  const VectorType noiseVarHP = GetDataSigmaSquared_HyperParameter();
  const VectorType noiseVarPrior = GetDataSigmaSquared_Prior();

  const MatrixType w0 = this->GetTemplateData()[imgIndex];
  const ScalarType ss0 = GetDataSigmaSquared(imgIndex);

  MatrixType w_old(w0), w(w0);
  ScalarType ss_old(ss0), ss(ss0);
  const unsigned int maxNbIterations = 20;
  const ScalarType tolerance = 1e-10;
  unsigned int k = 0;
  for (; k < maxNbIterations; ++k) {
    ss = (1 / (nbSubjects * nbVoxels + noiseVarHP[imgIndex]))
        * (S2 + dot_product(w_old, S4 * w_old) - 2 * dot_product(w_old, S3)
            + noiseVarHP[imgIndex] * noiseVarPrior[imgIndex]);
    w = inverse_sympd(S4 + ss * paramTempCovInv) * (S3 + ss * paramTempCovInv * paramTempMean);

    if (pow((ss - ss_old), 2) < tolerance && (w - w_old).sum_of_squares() < tolerance)
      break;
    else
      ss_old = ss;
    w_old = w;
  }
  assert(k < maxNbIterations - 1);

  SetDataSigmaSquared(VectorType(1, ss));
  this->SetTemplateData(MatrixListType(std::vector<MatrixType>(1, w)));

  /// Quick and dirty code for frozen template case.
  if (Superclass::m_FreezeTemplateFlag) {
    ss = (1 / (nbSubjects * nbVoxels + noiseVarHP[imgIndex]))
        * (S2 + dot_product(w0, S4 * w0) - 2 * dot_product(w0, S3)
            + noiseVarHP[imgIndex] * noiseVarPrior[imgIndex]);
    SetDataSigmaSquared(VectorType(1, ss));
    this->SetTemplateData(MatrixListType(std::vector<MatrixType>(1, w0)));
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s)
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::UpdateCovMomInverseAndComputeTerms(const std::vector<MatrixType> &momentas,
                                     ScalarType &covMomTerm,
                                     ScalarType &covMomPriorTerm) {
  const unsigned long nbSubjects = momentas.size();
  const unsigned int nbControlPoints = this->GetControlPoints().rows();

  /// Compute \sum_i \alpha_i\alpha_i^T.
  MatrixType AlphaPart(nbControlPoints * Dimension, nbControlPoints * Dimension, 0.0);
  for (unsigned int s = 0; s < nbSubjects; s++) {
    if (momentas[s].rows() != nbControlPoints)
      throw std::runtime_error(
          "Number of Control points and size of momenta mismatch in BayesianAtlas::ComputeOptimalCovMomInverse");

    MatrixType alpha_temp(nbControlPoints * Dimension, 1, 0.0);
    alpha_temp.set_column(0, this->Vectorize(momentas[s]));

    AlphaPart += (alpha_temp * alpha_temp.transpose());
  }

  /// Compute the inverse of the optimal momenta covariance matrix.
  const ScalarType covMomHP = GetCovarianceMomenta_HyperParameter();
  const MatrixType covMomPriorInv = GetCovarianceMomenta_Prior_Inverse();

  MatrixType aux = inverse(diagonal_matrix(nbControlPoints * Dimension, covMomHP) + covMomPriorInv * AlphaPart);
  aux *= (covMomHP + nbSubjects);

  MatrixType optimalCovMomInverse = aux * covMomPriorInv;
  SetCovarianceMomentaInverse(optimalCovMomInverse);

  /// Computes the log-likelihood terms.
  ScalarType logDetCovMom = -log_det(optimalCovMomInverse);

  covMomPriorTerm = -0.5 * covMomHP * (trace(aux) + logDetCovMom);
  covMomTerm = 0.0;
  for (unsigned int s = 0; s < nbSubjects; s++) {
    VectorType Moms = this->Vectorize(momentas[s]);
    covMomTerm -= dot_product(Moms, optimalCovMomInverse * Moms) + logDetCovMom;
  }
  covMomTerm *= 0.5;

}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::UpdateDataSigmaSquaredAndComputeTerms(const std::vector<std::vector<ScalarType>> &residuals,
                                        ScalarType &dataTerm,
                                        ScalarType &dataPriorTerm) {
  const unsigned int nbObjects = this->m_NumberOfObjects;
  const unsigned long nbSubjects = residuals.size();

  /// Compute optimal data sigma squared.
  const VectorType noiseVarHP = GetDataSigmaSquared_HyperParameter();
  const VectorType noiseVarPrior = GetDataSigmaSquared_Prior();

  VectorType optimalDataSigmaSquared(nbObjects, 0.0);
  for (int s = 0; s < nbSubjects; s++) {
    if (residuals[s].size() != nbObjects)
      throw std::runtime_error(
          "Number of objects in Residual Norms and Template mismatch in BayesianAtlas::ComptueOptimalDataSigmaSquared");

    for (int i = 0; i < nbObjects; i++)
      optimalDataSigmaSquared(i) += residuals[s][i];
  }

  for (int i = 0; i < nbObjects; i++)
    optimalDataSigmaSquared(i)
        = (optimalDataSigmaSquared[i] + noiseVarHP[i] * noiseVarPrior[i])
        / (noiseVarHP[i] + m_NoiseDimension[i] * nbSubjects);

  SetDataSigmaSquared(optimalDataSigmaSquared);

  /// Compute dataTerm and dataPriorTerm.
  dataPriorTerm = 0.0;
  dataTerm = 0.0;
  for (int i = 0; i < nbObjects; i++) {
    ScalarType OptimalDSS_i = optimalDataSigmaSquared(i);
    for (unsigned int s = 0; s < nbSubjects; s++)
      dataTerm -= 0.5 * (residuals[s][i] / OptimalDSS_i + m_NoiseDimension[i] * log(OptimalDSS_i));

    dataPriorTerm -= 0.5 * noiseVarHP[i] * (log(OptimalDSS_i) + noiseVarPrior[i] / OptimalDSS_i);
  }
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::ComputeSufficientStatisticsSubject(const MatrixType &momenta,
                                     const std::shared_ptr<const DeformableMultiObjectType> target,
                                     LinearVariableMapType &ss_i) const {
  /// Target.
  const unsigned int imgIndex = target->GetImageIndex();
  const MatrixType Ii = target->GetImageIntensityAndLandmarkPointCoordinates()[imgIndex];

  /// Deformed template.
  std::shared_ptr<DiffeosType> subjectDef = this->m_Def->Clone();
  subjectDef->SetDeformableMultiObject(this->m_Template);
  subjectDef->SetStartPositions(this->GetControlPoints());
  subjectDef->SetStartMomentas(momenta);
  subjectDef->Update();
  const std::shared_ptr<const ParametricImageType> deformedTemplate
      = std::static_pointer_cast<const ParametricImageType>(subjectDef->GetDeformedObject()->GetObjectList()[imgIndex]);
  const MatrixType deformedKernelField = deformedTemplate->GetPhotometricKernelField();

  /// Regularity.
  MatrixType alpha(momenta.rows() * Dimension, 1, 0.0);
  alpha.set_column(0, this->Vectorize(momenta));

  /// Compute the sufficient statistics.
  ss_i["S2"] = Ii.sum_of_squares();
  ss_i["S3"] = deformedKernelField.transpose() * Ii;
  ss_i["S4"] = deformedKernelField.transpose() * deformedKernelField;
  ss_i["S5"] = alpha * alpha.transpose();
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::InitializeCovMomPrior() {
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> ker = kfac->CreateKernelObject(this->m_Def->GetKernelType());
  ker->SetKernelWidth(this->m_Def->GetKernelWidth());

  const MatrixType controlPoints = this->GetControlPoints();
  const unsigned int nbCPs = controlPoints.rows();

  // Be careful: this should be consistent with the way this->Vectorize works!
  MatrixType covarianceMomenta_Prior_Inverse(nbCPs * Dimension, nbCPs * Dimension, 0.0);

  for (unsigned int i = 0; i < nbCPs; i++) {
    for (unsigned int j = 0; j < nbCPs; j++) {
      VectorType CPi = controlPoints.get_row(i);
      VectorType CPj = controlPoints.get_row(j);

      ScalarType k = ker->EvaluateKernel(CPi, CPj);
      MatrixType Aux = diagonal_matrix(Dimension, k);

      covarianceMomenta_Prior_Inverse.update(Aux, Dimension * i, Dimension * j);
      covarianceMomenta_Prior_Inverse.update(Aux, Dimension * j, Dimension * i);
    }
  }
  SetCovarianceMomenta_Prior_Inverse(covarianceMomenta_Prior_Inverse);
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::InitializeParametricTemplatePrior() {
  std::shared_ptr<ParametricImageType> parametricTemplate
      =
      std::static_pointer_cast<ParametricImageType>(this->m_Template->GetObjectList()[this->m_Template->GetImageIndex()]);
  const unsigned long nbPhotoCPs = parametricTemplate->GetNumberOfPhotometricControlPoints();

  /// Mean of the prior normal distribution.
  VectorType mean(nbPhotoCPs, 0.0);
  SetParametricTemplatePriorMean(mean);

  /// Inverse of the covariance of the prior normal distribution.
  MatrixType photometricCPs = parametricTemplate->GetPhotometricControlPoints();

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> ker = kfac->CreateKernelObject(this->m_Def->GetKernelType());
  ker->SetKernelWidth(parametricTemplate->GetPhotometricKernelWidth());

  MatrixType covInv(nbPhotoCPs, nbPhotoCPs);
  for (unsigned int i = 0; i < nbPhotoCPs; i++) {
    for (unsigned int j = 0; j < i; j++) {
      VectorType CPi = photometricCPs.get_row(i);
      VectorType CPj = photometricCPs.get_row(j);
      ScalarType k = ker->EvaluateKernel(CPi, CPj);

      covInv(i, j) = k;
      covInv(j, i) = k;
    }
    covInv(i, i) = 1;
  }

  SetParametricTemplatePriorCovarianceInverse(covInv);
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::InitializeControlPointsPrior() {
  /// Mean of the prior normal distribution.
  SetControlPointsPriorMean(this->GetControlPoints().vectorize());

  /// Sqrt of the covariance of the prior normal distribution.
  const std::shared_ptr<const DeformableMultiObjectType> temp = this->GetTemplate();
  const VectorType Xmin = temp->GetBoundingBox().get_column(0);
  const VectorType Xmax = temp->GetBoundingBox().get_column(1);
  const unsigned long nbCPs = this->GetControlPoints().rows();

  VectorType aux(Dimension);
  for (unsigned int dim = 0; dim < Dimension; ++dim)
    aux(dim) = 0.5 * (Xmax(dim) - Xmin(dim));

  MatrixType covSqrt(nbCPs * Dimension, nbCPs * Dimension, 0.0);
  for (unsigned int k = 0; k < nbCPs; ++k) {
    for (unsigned int dim = 0; dim < Dimension; ++dim) {
      covSqrt(Dimension * k + dim, Dimension * k + dim) = aux(dim);
    }
  }
  SetControlPointsPriorCovarianceSqrt(covSqrt);
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::InitializeMomentaRandomEffect() {
  /// Mean of the momenta random effect normal distribution.
  const unsigned long nbCPs = this->GetControlPoints().rows();
  SetMomentaRandomEffectMean(VectorType(nbCPs * Dimension, 0.0));

  /// Covariance of the momenta random effect normal distribution.
  SetCovarianceMomentaInverse(GetCovarianceMomenta_Prior_Inverse());
}

template<class ScalarType, unsigned int Dimension>
void
BayesianAtlas<ScalarType, Dimension>
::InitializeControlPointsRandomEffect() {
  /// Mean of the control points random effect normal distribution.
  SetControlPointsRandomEffectMean(this->GetControlPoints().vectorize());

  /// Covariance of the control points random effect normal distribution.
  if (GetControlPointsRandomEffectCovarianceInverse().n_elem() == 1) {
    std::cout << "No standard deviation given for the control points random effect : "
        "defaulting to 1/100th of the bounding box size on each coordinate." << std::endl;

    const std::shared_ptr<const DeformableMultiObjectType> temp = this->GetTemplate();
    const VectorType Xmin = temp->GetBoundingBox().get_column(0);
    const VectorType Xmax = temp->GetBoundingBox().get_column(1);
    const unsigned long nbCPs = this->GetControlPoints().rows();

    VectorType aux(Dimension);
    for (unsigned int dim = 0; dim < Dimension; ++dim)
      aux(dim) = 0.01 * (Xmax(dim) - Xmin(dim)); // Here 1/100th.

    MatrixType covSqrt(nbCPs * Dimension, nbCPs * Dimension, 0.0);
    for (unsigned int k = 0; k < nbCPs; ++k) {
      for (unsigned int dim = 0; dim < Dimension; ++dim) {
        covSqrt(Dimension * k + dim, Dimension * k + dim) = aux(dim);
      }
    }

    SetControlPointsRandomEffectCovarianceSqrt(covSqrt);
  }
}

template
class BayesianAtlas<double, 2>;
template
class BayesianAtlas<double, 3>;
