/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "LdaAtlas.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
LdaAtlas<ScalarType, Dimension>
::LdaAtlas() : Superclass(){
  this->SetLdaAtlasType();

  MatrixType G;
  MatrixType F;
  MatrixListType alpha;
  this->m_FixedEffects["F"] = F;
  this->m_FixedEffects["alpha"] = alpha;
  this->m_FixedEffects["G"] = G;

  ///F pop effect
  std::shared_ptr<NormalDistributionType> FDistribution = std::make_shared<NormalDistributionType>();
  this->m_Priors["F"] = FDistribution;

  ///G pop effect
  std::shared_ptr<AutomaticRelevanceDeterminationDistributionType> GDistribution = std::make_shared<AutomaticRelevanceDeterminationDistributionType>();
  this->m_Priors["G"] = GDistribution;

  ///alpha pop effect
  std::shared_ptr<NormalDistributionType> alphaDistribution = std::make_shared<NormalDistributionType>();
  this->m_Priors["alpha"] = alphaDistribution;

  ///Beta individual random effect
  std::shared_ptr<NormalDistributionType> betaDistribution = std::make_shared<NormalDistributionType>();
  this->m_IndividualRandomEffects["beta"] = betaDistribution;

  m_FreezeG = false;

}

template<class ScalarType, unsigned int Dimension>
LdaAtlas<ScalarType, Dimension>
::~LdaAtlas() {}

template<class ScalarType, unsigned int Dimension>
LdaAtlas<ScalarType, Dimension>
::LdaAtlas(const LdaAtlas &other) : Superclass(other) {
  m_NumberOfClasses = other.m_NumberOfClasses;
  m_IntraClassPCADimension = other.m_IntraClassPCADimension;
}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::SetFixedEffects(const LinearVariableMapType &map){
  const MatrixType G = GetG();
  const MatrixType F = GetF();

  Superclass::SetFixedEffects(map);

  if (m_FreezeG) SetG(G);
  if (m_FreezeF) SetF(F);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::Update() {
  Superclass::Update();
  if (m_NumberOfClasses > 1){
    this->InitializeF();
    this->InitializeAlpha();
  }
  else
    std::cout << "Lda model with one class : Bayesian PGA model" << std::endl;
  this->InitializeG();
  this->InitializeBeta();

  ///Initializing the diffeo bounding box (mostly for simulate-lda-atlas)
  this->InitializeBoundingBox();
}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::InitializeF(){

  MatrixType F = GetF();
  ///Initialization of F:

  if (F.size()==0) {
    F = MatrixType(this->GetControlPoints().rows() * Dimension, m_NumberOfClasses - 1, 1.);

    std::shared_ptr<NormalDistributionType> FDistribution = std::make_shared<NormalDistributionType>();
    VectorType mean(this->GetControlPoints().rows() * Dimension, 0.);
    FDistribution->SetMean(mean);
    MatrixType cov = diagonal_matrix(this->GetControlPoints().rows() * Dimension, 0.2);
    FDistribution->SetCovariance(cov);

    for (unsigned int i = 0; i < m_NumberOfClasses - 1; ++i)
      F.set_column(i, FDistribution->Sample());

    this->m_FixedEffects["F"] = F;
  }

  ///F pop effect
  std::shared_ptr<NormalDistributionType> FDistribution = std::make_shared<NormalDistributionType>();
  VectorType mean(this->GetControlPoints().rows() * Dimension, 0.);
  FDistribution->SetMean(mean);
  MatrixType cov = diagonal_matrix(this->GetControlPoints().rows() * Dimension, 1.);
  FDistribution->SetCovariance(cov);
  this->m_Priors["F"] = FDistribution;
}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::InitializeG(){

  MatrixType G = GetG();


  if (G.rows()==0) {
    ///Initialization of G:
    G = MatrixType(this->GetControlPoints().rows() * Dimension, m_IntraClassPCADimension, 1.);

    std::shared_ptr<NormalDistributionType> alphaDistribution = std::make_shared<NormalDistributionType>();
    MatrixType mean(this->GetControlPoints().rows() * Dimension, 1, 0.);
    alphaDistribution->SetMean(mean);
    MatrixType cov = diagonal_matrix(this->GetControlPoints().rows() * Dimension, 0.002);
    alphaDistribution->SetCovariance(cov);

    for (unsigned int i = 0; i < m_IntraClassPCADimension; ++i)
      G.set_column(i, alphaDistribution->Sample());

    this->m_FixedEffects["G"] = G;
  }

  ///G pop effect
  std::shared_ptr<AutomaticRelevanceDeterminationDistributionType> GDistribution = std::make_shared<AutomaticRelevanceDeterminationDistributionType>();
  ///Initializing gamma with unit variance for each column.
  VectorType gamma(m_IntraClassPCADimension, 1.);
  GDistribution->SetGamma(gamma);
  this->m_Priors["G"] = GDistribution;
}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::InitializeAlpha(){

  MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);


  if (alpha.size()==0) {
    ///alpha pop effect
    std::shared_ptr<NormalDistributionType> alphaDistribution = std::make_shared<NormalDistributionType>();
    MatrixType mean(m_NumberOfClasses - 1, 1, 0.);
    alphaDistribution->SetMean(mean);
    MatrixType cov = diagonal_matrix(m_NumberOfClasses - 1, 1.);
    alphaDistribution->SetCovariance(cov);
    this->m_Priors["alpha"] = alphaDistribution;

    ///Initialization of alpha:
    alpha.resize(m_NumberOfClasses);
    for (unsigned int i = 0; i < m_NumberOfClasses; ++i) {
      VectorType aux = alphaDistribution->Sample();
      alpha[i] = MatrixType(aux);
    }
    this->m_FixedEffects["alpha"] = alpha;
  }
}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::InitializeBeta(){
  ///This distribution on beta is best seen as a regularization.
  std::shared_ptr<NormalDistributionType> betaDistribution = std::make_shared<NormalDistributionType>();
  MatrixType mean(m_IntraClassPCADimension, 1, 0.);
  betaDistribution->SetMean(mean);
  MatrixType cov = diagonal_matrix(m_IntraClassPCADimension, 1.);
  betaDistribution->SetCovariance(cov);
  this->m_IndividualRandomEffects["beta"] = betaDistribution;
}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::InitializeBoundingBox(){
  /// Initialization.
  MatrixType boundingBox = this->m_Template->GetBoundingBox();
  MatrixType controlPoints = this->GetControlPoints();

  /// Compute the model bounding box.
  for (unsigned int i = 0; i < controlPoints.rows(); i++) {
    for (unsigned int dim = 0; dim < Dimension; dim++) {
      boundingBox(dim, 0) =
          (this->m_BoundingBox(dim, 0) < controlPoints(i, dim) ? this->m_BoundingBox(dim, 0) : controlPoints(i, dim));
      boundingBox(dim, 1) =
          (this->m_BoundingBox(dim, 1) > controlPoints(i, dim) ? this->m_BoundingBox(dim, 1) : controlPoints(i, dim));
    }
  }

  /// Set the related bounding boxes.
 this->m_Def->SetDataDomain(boundingBox);
}

template<class ScalarType, unsigned int Dimension>
bool
LdaAtlas<ScalarType, Dimension>
::Sample(LongitudinalDataSetType *dataSet,
         LinearVariableMapType &popRER,
         LinearVariablesMapType &indRER) const {
   ///Variables declaration
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType >> > samples(m_NumberOfClasses * m_NumberOfSubjectsPerClass);
  std::vector<MatrixType> betas(m_NumberOfClasses * m_NumberOfSubjectsPerClass);
  std::shared_ptr<NormalDistributionType> betaDistribution = std::make_shared<NormalDistributionType>();
  MatrixType mean(m_IntraClassPCADimension,1,0.);
  ///Distribution for noise component (pga part)
  betaDistribution->SetMean(mean);
  MatrixType cov = diagonal_matrix(m_IntraClassPCADimension, 0.1);
  betaDistribution->SetCovariance(cov);

  std::vector<MatrixType> classesMomenta(m_NumberOfClasses);
  MatrixType F = GetF();
  MatrixType G = GetG();
  MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
  std::vector<unsigned int> sizeParams(2);
  sizeParams[0] = this->GetControlPoints().rows();
  sizeParams[1] = Dimension;
  this->m_Def->SetDeformableMultiObject(this->GetTemplate());


  for (unsigned int c=0;c<m_NumberOfClasses;++c){
    if (m_NumberOfClasses>1) {
      VectorType aux = F * vectorize(LinearVariableType(alpha[c]));
      classesMomenta[c] = recast<MatrixType>(unvectorize(aux, sizeParams));
    } else
      classesMomenta[c] = MatrixType(this->GetControlPoints().rows(), Dimension, 0.);
    for (unsigned int i=0;i<m_NumberOfSubjectsPerClass;++i){
      ///Sampling a noise component beta
      betas[c*m_NumberOfSubjectsPerClass + i] = MatrixType(m_IntraClassPCADimension, 1., 0.);
      betas[c*m_NumberOfSubjectsPerClass + i].set_column(0, betaDistribution->Sample());

      ///Generating the corresponding momenta for this class
      VectorType vectorizedMomentaForSubject = classesMomenta[c].vectorize() + G * betas[c*m_NumberOfSubjectsPerClass + i].vectorize();
      MatrixType momentaForSubject = recast<MatrixType>(unvectorize(vectorizedMomentaForSubject, sizeParams));

      ///Feeding the corresponding momenta to the diffeo
      this->m_Def->SetStartMomentas(momentaForSubject);
      this->m_Def->SetStartPositions(this->GetControlPoints());
      this->m_Def->Update();
      
      ///Storing output
      samples[c*m_NumberOfSubjectsPerClass + i].push_back(this->m_Def->GetDeformedObject());
    }
  }

  ///Output
  dataSet->SetDeformableMultiObjects(samples);
  popRER["F"] = F;
  popRER["G"] = G;
  popRER["alpha"] = alpha;
  indRER["beta"] = betas;

  return true;
}

template<class ScalarType, unsigned int Dimension>
bool
LdaAtlas<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<ScalarType>> &residuals) {

  const std::vector<MatrixType> betas = recast<MatrixType>(indRER.at("beta"));
  const CrossSectionalDataSetType *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  const std::vector<unsigned int> classes = crossSectionalDataSet->GetClasses();
  std::vector<MatrixType> momentas = this->GetMomentasFromLdaDescription(betas, classes);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();
  return Superclass::ComputeResiduals(this->GetControlPoints(), momentas, target, residuals);
}

template<class ScalarType, unsigned int Dimension>
ScalarType
LdaAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                               const LinearVariableMapType &popRER,
                               const LinearVariablesMapType &indRER) {
  ScalarType out = 0.0;

  const CrossSectionalDataSetType *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  const std::vector<unsigned int> classes = crossSectionalDataSet->GetClasses();

  ///List of terms in the likelihood : noise variance, attachment, F, G, alpha, beta.

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

  if (m_NumberOfClasses > 1){
    ///Regularity on classes positions:
    const MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
    assert(m_NumberOfClasses == alpha.size());
    for (unsigned int c=0; c<m_NumberOfClasses;++c)
      out += this->m_Priors["alpha"]->ComputeLogLikelihood(alpha[c]);

    ///Regularity on F
    const MatrixType F = GetF();
    for (unsigned int c=0; c<F.cols(); ++c)
      out += this->m_Priors["F"]->ComputeLogLikelihood(F.get_column(c));
  }

  ///Regularity on the intra-class positions:
  std::vector<MatrixType> betas = recast<MatrixType>(indRER.at("beta"));
  assert(nbSubjects == betas.size());

  for (unsigned int s=0; s<nbSubjects;++s)
    out += this->m_IndividualRandomEffects["beta"]->ComputeLogLikelihood(betas[s]);


  ///Regularity on G
  const MatrixType G = GetG();
  out += this->m_Priors["G"]->ComputeLogLikelihood(G);

  return out;
}

template<class ScalarType, unsigned int Dimension>
bool
LdaAtlas<ScalarType, Dimension>
::UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    VectorType &logLikelihoodTerms) {

  const CrossSectionalDataSetType *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  const std::vector<unsigned int> classes = crossSectionalDataSet->GetClasses();

  std::vector<std::vector<ScalarType>> residuals;
  unsigned int nbSubjects = residuals.size();
  bool oob = this->ComputeResiduals(dataSet, popRER, indRER, residuals);

  ScalarType gammaTerm, dataTerm, FTerm, GTerm, alphaTerm, betaTerm, noiseVarianceTerm;

  UpdateGamma(gammaTerm, nbSubjects);
  UpdateNoiseVarianceAndComputeTerms(residuals, dataTerm, GTerm, FTerm, alphaTerm, betaTerm, noiseVarianceTerm, indRER, classes);

  logLikelihoodTerms.set_size(2);
  logLikelihoodTerms.fill(0.0);

  logLikelihoodTerms[1] += gammaTerm;
  logLikelihoodTerms[1] += GTerm;
  logLikelihoodTerms[1] += FTerm;
  logLikelihoodTerms[1] += alphaTerm;
  logLikelihoodTerms[1] += betaTerm;
  logLikelihoodTerms[1] += noiseVarianceTerm;

  logLikelihoodTerms[0] += dataTerm;

  return oob; // Out of box flag.
}





template<class ScalarType, unsigned int Dimension>
ScalarType
LdaAtlas<ScalarType, Dimension>
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
void
LdaAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad) {

  const CrossSectionalDataSetType
      *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();
  const std::vector<unsigned int> classes = crossSectionalDataSet->GetClasses();

  const std::vector<MatrixType> betas = recast<MatrixType>(indRER.at("beta"));
  const unsigned long nbSubjects = target.size();
  const MatrixType G = GetG();
  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const VectorType gamma = GetGamma();
  const VectorType noiseVariance = GetNoiseVariance();

  printMatrix<ScalarType>("G", G);
  std::cout << "Noise variance : "<< noiseVariance[0] << std::endl;

  const unsigned int nbObjects = this->m_NumberOfObjects;
  assert(betas.size() == nbSubjects);

  MatrixType cp = this->GetControlPoints();

  if (betas.size() != nbSubjects)
    throw std::runtime_error("Number of Momentas and target multi-object mismatch in LdaAtlas::ComputeLikelihood");

  if (betas.size() == 0)
    throw std::runtime_error("Cannot run LDA with no subject");

  std::vector<MatrixType> momentas = this->GetMomentasFromLdaDescription(betas, classes);

  MatrixType gradPos;
  MatrixListType gradTempL_L2, gradTempL_Sob;
  std::vector<MatrixType> gradBetas(betas.size());
  MatrixType gradG(G.rows(), G.cols(), 0.);
  std::vector<MatrixType> gradMom(nbSubjects);
  MatrixListType gradAlpha;
  MatrixType gradF;

  this->UpdateDeformationAndKernelDataDomain(target);

  if (!this->m_FreezeTemplateFlag) {
    this->ComputeDataTermGradient(cp, momentas, noiseVariance, target, gradPos, gradMom, gradTempL_L2);
    gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
  } else { this->ComputeDataTermGradient(cp, momentas, noiseVariance, target, gradPos, gradMom); }

  if (m_NumberOfClasses > 1) {
    const MatrixType F = GetF();
    const MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
    gradAlpha.resize(m_NumberOfClasses);
    gradF = MatrixType(F.rows(), F.cols(), 0.);
    if (not(m_FreezeF)) {
      ///Gradient for F, data term part:
      for (unsigned int s = 0; s < nbSubjects; ++s)
        gradF += this->VectorizeAndMultiply(gradMom[s], alpha[classes[s]]);

      ///Gradient for F : regularity part
      gradF -= F;
    }

    ///Gradient for classPositions : regularity
    for (unsigned int c = 0; c < m_NumberOfClasses; ++c) {
      gradAlpha[c] = -1. * alpha[c];
    }

    for (unsigned int s = 0; s < nbSubjects; ++s) {
      ///This is the data term
      MatrixType aux = this->MultiplyMatrixTransposeByVector(F, gradMom[s]);
      gradAlpha[classes[s]] += aux;
    }

  }

  if (not(m_FreezeG)) {
    ///Gradient for G, data term part.
    for (unsigned int s = 0; s < nbSubjects; ++s) {
      MatrixType aux = this->VectorizeAndMultiply(gradMom[s], betas[s]);
      gradG = aux;
    }

    ///Gradient for G : regularity part.
    MatrixType gammaMatrix = diagonal_matrix(gamma);
    printMatrix<ScalarType>("gamma matrix", gammaMatrix);
    gradG -= (G * gammaMatrix);
  }

  std::vector<unsigned int> sizeParams;
  vectorize(LinearVariableType(betas[0]), sizeParams);

  ///Gradient for intraclassPositions : data term
  for (unsigned int s = 0; s < nbSubjects; ++s) {
    MatrixType aux = this->MultiplyMatrixTransposeByVector(G, gradMom[s]);
    gradBetas[s] = aux;
    gradBetas[s] -= betas[s];
  }


  if (!this->m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  if (!this->m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }

  indGrad["beta"] = gradBetas;
  popGrad["G"] = gradG;

  if (m_NumberOfClasses >1) {
    popGrad["alpha"] = gradAlpha;
    popGrad["F"] = gradF;
  }

  std::cout << "Squared norm gradient momenta :" << MatrixListType(gradMom).sum_of_squares() << std::endl;
  std::cout << "Squared norm gradient beta : " << MatrixListType(gradBetas).sum_of_squares() << std::endl;
  if (m_NumberOfClasses >1) {
    std::cout << "Squared norm grad alpha : " << gradAlpha.sum_of_squares() << std::endl;
    std::cout << "Squared norm grad F : " << gradF.sum_of_squares() << std::endl;
  }
  std::cout << "Squared norm grad G : " << gradG.sum_of_squares() << std::endl;



}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad,
                                       VectorType &gradSquaredNorms) {


  const CrossSectionalDataSetType
      *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();
  const std::vector<unsigned int> classes = crossSectionalDataSet->GetClasses();

  const std::vector<MatrixType> betas = recast<MatrixType>(indRER.at("beta"));
  const unsigned long nbSubjects = target.size();
  const MatrixType G = GetG();
  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const VectorType gamma = GetGamma();
  const VectorType noiseVariance = GetNoiseVariance();

  printMatrix<ScalarType>("G", G);
  std::cout << "Noise variance : "<< noiseVariance[0] << std::endl;

  const unsigned int nbObjects = this->m_NumberOfObjects;
  assert(betas.size() == nbSubjects);

  MatrixType cp = this->GetControlPoints();

  if (betas.size() != nbSubjects)
    throw std::runtime_error("Number of Momentas and target multi-object mismatch in LdaAtlas::ComputeLikelihood");

  if (betas.size() == 0)
    throw std::runtime_error("Cannot run LDA with no subject");

  std::vector<MatrixType> momentas = this->GetMomentasFromLdaDescription(betas, classes);

  MatrixType gradPos;
  MatrixListType gradTempL_L2, gradTempL_Sob;
  std::vector<MatrixType> gradBetas(betas.size());
  MatrixType gradG(G.rows(), G.cols(), 0.);
  std::vector<MatrixType> gradMom(nbSubjects);
  MatrixListType gradAlpha;
  MatrixType gradF;

  this->UpdateDeformationAndKernelDataDomain(target);

  if (!this->m_FreezeTemplateFlag) {
    this->ComputeDataTermGradient(cp, momentas, noiseVariance, target, gradPos, gradMom, gradTempL_L2);
    gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
  } else { this->ComputeDataTermGradient(cp, momentas, noiseVariance, target, gradPos, gradMom); }

  if (m_NumberOfClasses > 1) {
    const MatrixType F = GetF();
    const MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
    gradAlpha.resize(m_NumberOfClasses);
    gradF = MatrixType(F.rows(), F.cols(), 0.);
    if (not(m_FreezeF)) {
      ///Gradient for F, data term part:
      for (unsigned int s = 0; s < nbSubjects; ++s)
        gradF += this->VectorizeAndMultiply(gradMom[s], alpha[classes[s]]);

      ///Gradient for F : regularity part
      gradF -= F;
    }

    ///Gradient for classPositions : regularity
    for (unsigned int c = 0; c < m_NumberOfClasses; ++c) {
      gradAlpha[c] = -1. * alpha[c];
    }

    for (unsigned int s = 0; s < nbSubjects; ++s) {
      ///This is the data term
      MatrixType aux = this->MultiplyMatrixTransposeByVector(F, gradMom[s]);
      gradAlpha[classes[s]] += aux;
    }

  }

  if (not(m_FreezeG)) {
    ///Gradient for G, data term part.
    for (unsigned int s = 0; s < nbSubjects; ++s) {
      MatrixType aux = this->VectorizeAndMultiply(gradMom[s], betas[s]);
      gradG = aux;
    }

    ///Gradient for G : regularity part.
    MatrixType gammaMatrix = diagonal_matrix(gamma);
    printMatrix<ScalarType>("gamma matrix", gammaMatrix);
    gradG -= (G * gammaMatrix);
  }

  std::vector<unsigned int> sizeParams;
  vectorize(LinearVariableType(betas[0]), sizeParams);

  ///Gradient for intraclassPositions : data term
  for (unsigned int s = 0; s < nbSubjects; ++s) {
    MatrixType aux = this->MultiplyMatrixTransposeByVector(G, gradMom[s]);
    gradBetas[s] = aux;
    gradBetas[s] -= betas[s];
  }


  if (!this->m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  if (!this->m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }

  indGrad["beta"] = gradBetas;
  popGrad["G"] = gradG;

  if (m_NumberOfClasses >1) {
    popGrad["alpha"] = gradAlpha;
    popGrad["F"] = gradF;
  }

  std::cout << "Squared norm gradient momenta :" << MatrixListType(gradMom).sum_of_squares() << std::endl;

  gradSquaredNorms.set_size(0);
  gradSquaredNorms.push_back(MatrixListType(gradBetas).sum_of_squares());
  std::cout << "Squared norm gradient beta : " << MatrixListType(gradBetas).sum_of_squares() << std::endl;
  if (m_NumberOfClasses >1) {
    gradSquaredNorms.push_back(gradAlpha.sum_of_squares());
    std::cout << "Squared norm grad alpha : " << gradAlpha.sum_of_squares() << std::endl;
  }
  if (!this->m_FreezeControlPointsFlag) { gradSquaredNorms.push_back(gradPos.sum_of_squares()); }
  if (m_NumberOfClasses >1){
    gradSquaredNorms.push_back(gradF.sum_of_squares());
    std::cout << "Squared norm grad F : " << gradF.sum_of_squares() << std::endl;
  }
  gradSquaredNorms.push_back(gradG.sum_of_squares());
  std::cout << "Squared norm grad G : " << gradG.sum_of_squares() << std::endl;
  if (!this->m_FreezeTemplateFlag) { gradSquaredNorms.push_back(dot_product(gradTempL_L2, gradTempL_Sob)); }

}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::UpdateGamma(ScalarType& gammaTerm, unsigned int& nbSubjects){

  const MatrixType G = GetG();
  const unsigned int nbCPs = this->GetControlPoints().rows();


  VectorType newGamma(m_IntraClassPCADimension, 0.);
  for (unsigned int i=0; i<m_IntraClassPCADimension; ++i)
    if (G.get_column(i).sum_of_squares() != 0.) {
      assert(G.get_column(i).size() == nbCPs * Dimension);
      newGamma[i] = nbCPs * Dimension / G.get_column(i).sum_of_squares();
    }
    else
      newGamma[i] = 1.;

  this->SetGamma(newGamma);

  gammaTerm = 0;
  for (unsigned int i=0; i<newGamma.size(); ++i)
    gammaTerm += log(newGamma.get(i)/ 2*M_PI);

  gammaTerm *= 0.5 * nbSubjects * Dimension * this->GetControlPoints().rows();


}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::UpdateNoiseVarianceAndComputeTerms(std::vector<std::vector<ScalarType>> const& residuals, ScalarType& dataTerm,
                                     ScalarType& GTerm, ScalarType& FTerm, ScalarType& alphaTerm,
                                     ScalarType& betaTerm, ScalarType& noiseVarianceTerm, const LinearVariablesMapType& indRER, std::vector<unsigned int> const& classes){


  const unsigned int nbObjects = this->m_NumberOfObjects;
  const unsigned int nbSubjects = residuals.size();

  alphaTerm = 0.;
  betaTerm = 0.;
  GTerm = 0.;
  FTerm = 0.;
  noiseVarianceTerm = 0.;
  dataTerm = 0.;

  MatrixType G = GetG();
  MatrixType F;
  if (m_NumberOfClasses>1)
    F = GetF();

  VectorType newNoiseVariance(nbObjects, 0.);

  for (unsigned int o=0;o<nbObjects;++o) {
    for (unsigned int s = 0; s < nbSubjects; ++s)
      newNoiseVariance[o] += residuals[s][o];
    newNoiseVariance[o] /= (m_NoiseDimension[o] * nbSubjects);
  }

  this->SetNoiseVariance(newNoiseVariance);

  for (unsigned int s = 0;s<nbSubjects; ++s)
    for (unsigned int o=0;o<nbObjects; ++o)
      dataTerm -= 0.5 * residuals[s][o] / newNoiseVariance(o);

  std::cout << "Sum of residuals : " << (dataTerm *2. * newNoiseVariance(0)) << std::endl;

  for (unsigned int o=0; o<nbObjects; ++o)
    noiseVarianceTerm -= nbSubjects * m_NoiseDimension[o] * log(newNoiseVariance[o]);

  GTerm = this->m_Priors["G"]->ComputeLogLikelihood(G);

  if (m_NumberOfClasses > 1) {
    for (unsigned int c=0;c<F.cols();++c)
      FTerm += this->m_Priors["F"]->ComputeLogLikelihood(F.get_column(c));

    ///Regularity on classes positions:
    const MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
    assert(m_NumberOfClasses = alpha.size());

    ///Go through the map. For each label, classesPositions[label] is the vector indicating the position of the class.
    for (unsigned int c = 0; c <m_NumberOfClasses; ++c)
      alphaTerm += this->m_Priors["alpha"]->ComputeLogLikelihood(alpha[c]);
  }

  ///Regularity on the intra-class positions:
  const std::vector<MatrixType> betas = recast<MatrixType>(indRER.at("beta"));
  assert(nbSubjects == betas.size());

//  printMatrix<ScalarType>("first beta", betas[0]);

  for (unsigned int s=0; s<nbSubjects;++s)
    betaTerm += this->m_IndividualRandomEffects["beta"]->ComputeLogLikelihood(betas[s]);

//  std::cout << "Noise variance term : " << noiseVarianceTerm <<  "GTerm : " << GTerm << " Beta term : " << betaTerm  << " alpha term : " << alphaTerm << " FTerm = " << FTerm << std::endl;

}


template<class ScalarType, unsigned int Dimension>
std::vector<MatrixType>
LdaAtlas<ScalarType, Dimension>
::GetMomentasFromLdaDescription(std::vector<MatrixType> const& betas, std::vector<unsigned int> classes) const{

  if (m_NumberOfClasses>1) {
    MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
    MatrixType F = GetF();
    MatrixType G = GetG();

    std::vector<MatrixType> momentas(betas.size());
    if (betas.size() == 0)
      return momentas;

    std::vector<unsigned int> sizeParams(2);
    sizeParams[0] = this->GetControlPoints().rows();
    sizeParams[1] = Dimension;

    for (unsigned int subj = 0; subj < momentas.size(); ++subj) {
      VectorType aux = F * vectorize(LinearVariableType(alpha[classes[subj]])) + G * betas[subj].vectorize();
      momentas[subj] = recast<MatrixType>(unvectorize(aux, sizeParams));
    }
    return momentas;
  }

  else {
    MatrixType G = recast<MatrixType>(this->m_FixedEffects["G"]);
    std::vector<MatrixType> momentas(betas.size());

    if (betas.size() == 0)
      return momentas;

    std::vector<unsigned int> sizeParams(2);
    sizeParams[0] = this->GetControlPoints().rows();
    sizeParams[1] = Dimension;

    for (unsigned int subj = 0; subj < momentas.size(); ++subj) {
      VectorType aux = G * betas[subj].vectorize();
      momentas[subj] = recast<MatrixType>(unvectorize(aux, sizeParams));
    }

    return momentas;
  }
}

template<class ScalarType, unsigned int Dimension>
MatrixType
LdaAtlas<ScalarType, Dimension>
::VectorizeAndMultiply(MatrixType const& gradMom, MatrixType const& alpha) {

  VectorType gradMomVectorized = gradMom.vectorize();
  VectorType alphaVectorized = alpha.vectorize();

  MatrixType out(gradMom.cols()*gradMom.rows(), alpha.cols()*alpha.rows(), 0.);
  for (unsigned int i=0; i<gradMom.cols()*gradMom.rows(); ++i)
    for (unsigned j=0; j<alpha.cols()*alpha.rows(); ++j)
      out(i,j) = gradMomVectorized.get(i) * alphaVectorized.get(j);

  return out;
}

template<class ScalarType, unsigned int Dimension>
MatrixType
LdaAtlas<ScalarType, Dimension>
::MultiplyMatrixTransposeByVector(MatrixType const& F, MatrixType const& gradMom) {

  VectorType vectorizedMom = vectorize(LinearVariableType(gradMom));
  MatrixType out = F.transpose() * vectorizedMom;
  return out;

}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::Write(const LongitudinalDataSetType *const dataSet,
        LinearVariableMapType const &popRER,
        LinearVariablesMapType const &indRER) const {

  const std::string outputDir = def::utils::settings.output_dir;

  const CrossSectionalDataSetType *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();
  const std::vector<unsigned int> classes = crossSectionalDataSet->GetClasses();
  const std::vector<MatrixType> betas = recast<MatrixType>(indRER.at("beta"));
  const MatrixType G = GetG();

  const unsigned int nbObjects = this->m_NumberOfObjects;
  unsigned int nbSubjects = betas.size();

  MatrixType cp = this->GetControlPoints();

  if (betas.size() != nbSubjects)
    throw std::runtime_error("Number of Momentas and target multi-object mismatch in LdaAtlas::ComputeLikelihood");

  if (betas.size() == 0)
    throw std::runtime_error("Cannot run LDA with no subject");

  std::vector<MatrixType> momentas = this->GetMomentasFromLdaDescription(betas, classes);

  MatrixListType beta(betas.size());
  for (unsigned int i=0;i<betas.size();++i){
    beta[i] = betas[i];
  }


  ///Les matrices F et G
  if (m_NumberOfClasses > 1) {
    const MatrixType F = GetF();
    const MatrixListType alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
    writeMatrixDLM<ScalarType>(outputDir+"F.txt", F);
    writeMultipleMatrixDLM<ScalarType>(outputDir+"alpha.txt", alpha);
  }

  writeMatrixDLM<ScalarType>(outputDir+"G.txt", G);
  writeMultipleMatrixDLM<ScalarType>(outputDir+"betas.txt", beta);

  ///Writing the gammas
  MatrixType gamma = MatrixType(GetGamma());
  writeMatrixDLM<ScalarType>(outputDir+"gamma.txt", gamma);


  ///Need to write : les objets reconstitués (toute la séquence ?)
  Superclass::WriteAtlasToSubjectDeformations(momentas, dataSet);
  ///le template
  Superclass::WriteTemplateData();


  ///les positions des classes, et les objets moyens de chaque classe
  WriteClassesInfo(G);

    /// Write control points.
  std::ostringstream oss1;
  oss1 << outputDir << "Atlas_ControlPoints.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss1.str().c_str(), this->GetControlPoints());

  /// Write initial momentas.
  std::ostringstream oss2;
  oss2 << outputDir << "Atlas_Momentas.txt" << std::ends;
  writeMultipleMatrixDLM<ScalarType>(oss2.str().c_str(), momentas);

}

template<class ScalarType, unsigned int Dimension>
void
LdaAtlas<ScalarType, Dimension>
::WriteClassesInfo(MatrixType const& G) const {

  std::vector<unsigned int> sizeParams(2);
  sizeParams[0] = this->GetControlPoints().rows();
  sizeParams[1] = Dimension;

  MatrixType F;
  MatrixListType alpha;
  if (m_NumberOfClasses > 1) {
    F = GetF();
    alpha = recast<MatrixListType>(this->m_FixedEffects["alpha"]);
  }

  const std::string outputDir = def::utils::settings.output_dir;

  std::vector<DeformableMultiObjectType> classesTemplates = std::vector<DeformableMultiObjectType>(m_NumberOfClasses);
  std::vector<MatrixType> classesMomenta = std::vector<MatrixType>(m_NumberOfClasses);

  this->m_Def->SetDeformableMultiObject(this->GetTemplate());
  this->m_Def->SetStartPositions(this->GetControlPoints());



  for (unsigned int i=0;i<m_NumberOfClasses;++i){
    if (m_NumberOfClasses > 1)
      classesMomenta[i] = recast<MatrixType>(unvectorize(F * vectorize(LinearVariableType(alpha[i])), sizeParams));
    else
      classesMomenta[i] = MatrixType(this->GetControlPoints().rows(), Dimension, 0.);
    this->m_Def->SetStartMomentas(classesMomenta[i]);
    this->m_Def->Update();
    std::shared_ptr<DeformableMultiObjectType> classTemplate = this->m_Def->GetDeformedObject();
    std::vector<std::string> names(this->m_TemplateObjectsName.size());
    for (int o = 0; o < this->m_TemplateObjectsName.size(); o++) {
      std::ostringstream oss;
      oss << outputDir << "classTemplate_" << i << this->m_TemplateObjectsNameExtension[o] << std::ends;
      names[o] = oss.str();
    }
    classTemplate->WriteMultiObject(names);
  }

  unsigned int numberOfTimePoints = this->m_Def->GetNumberOfTimePoints();


  std::vector<std::vector<std::shared_ptr<typename LdaAtlas<ScalarType,Dimension>::DeformableMultiObjectType>>> components = GetComponents();


  for (unsigned int i=0;i<m_IntraClassPCADimension;++i) {
    for (unsigned int t = 0; t < components[i].size(); ++t) {
      std::vector<std::string> names(1);
      std::ostringstream oss;
      oss << outputDir << "principal_component_" << i << "_time_" << t << ".vtk" << std::ends;
      names[0] = oss.str();
      components[i][t]->WriteMultiObject(names);
    }
  }
}

///Method used in the sampling.
template<class ScalarType, unsigned int Dimension>
std::vector<std::vector<std::shared_ptr<typename LdaAtlas<ScalarType,Dimension>::DeformableMultiObjectType>>>
LdaAtlas<ScalarType, Dimension>
::GetComponents() const {

  std::vector<unsigned int> sizeParams(2);
  sizeParams[0] = this->GetControlPoints().rows();
  sizeParams[1] = Dimension;
  unsigned int numberOfTimePoints = this->m_Def->GetNumberOfTimePoints();

  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> components(m_IntraClassPCADimension);
  MatrixType G = GetG();


  for (unsigned int i=0;i<m_IntraClassPCADimension;++i)
  {
    std::cout << i << std::endl;
    ///We need to write the trajectories corresponding to the different columns of G.
    MatrixType principalComponent = recast<MatrixType>(unvectorize(G.get_column(i), sizeParams));
    ///Shooting backward first
    this->m_Def->SetStartMomentas(-1.*principalComponent);
    this->m_Def->SetStartPositions(this->GetControlPoints());
    this->m_Def->Update();
    for (unsigned int t=numberOfTimePoints-1;t>0;--t) {
      components[i].push_back(this->m_Def->GetDeformedObjectAt(t));
    }
    ///Shooting forward now.
    this->m_Def->SetStartMomentas(principalComponent);
    this->m_Def->Update();
    for (unsigned int t=0;t<numberOfTimePoints;++t) {
      components[i].push_back(this->m_Def->GetDeformedObjectAt(t));
    }
  }

  return components;
}


template class LdaAtlas<double, 2>;
template class LdaAtlas<double, 3>;
