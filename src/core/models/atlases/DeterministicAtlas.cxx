/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeterministicAtlas.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
DeterministicAtlas<ScalarType, Dimension>
::DeterministicAtlas() : Superclass() {
  this->SetDeterministicAtlasType();
  m_UseRKHSNormForRegularization = true;
}

template<class ScalarType, unsigned int Dimension>
DeterministicAtlas<ScalarType, Dimension>
::~DeterministicAtlas() {}

template<class ScalarType, unsigned int Dimension>
DeterministicAtlas<ScalarType, Dimension>
::DeterministicAtlas(const DeterministicAtlas &other) : Superclass(other) {
  m_CovarianceMomentaInverse = other.m_CovarianceMomentaInverse;
  m_DataSigmaSquared = other.m_DataSigmaSquared;
  m_UseRKHSNormForRegularization = other.m_UseRKHSNormForRegularization;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class ScalarType, unsigned int Dimension>
void
DeterministicAtlas<ScalarType, Dimension>
::SetCovarianceMomentaInverse(std::string const &fn) {
  if (strlen(fn.c_str())) {
    m_CovarianceMomentaInverse = readMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using the inverse matrix of the momenta covariance of size " << m_CovarianceMomentaInverse.rows()
              << " x " << m_CovarianceMomentaInverse.cols() << " from file: " << fn << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
DeterministicAtlas<ScalarType, Dimension>
::Update() {
  Superclass::Update();

  const unsigned int nbControlPoints = this->GetControlPoints().rows();

  if (!m_UseRKHSNormForRegularization) {
    if (m_CovarianceMomentaInverse.rows() == 0) {
      std::cout << "Warning: no covariance momenta matrix set for computing regularization term. "
          "We will use the RKHS norm for computing the regularization term" << std::endl;
      m_UseRKHSNormForRegularization = true;
    }

    if (m_CovarianceMomentaInverse.rows() != Dimension * nbControlPoints) {
      std::cout << " Warning: the size of the covariance momenta matrix that you set does not match "
          "the number of control points. Matrix set to identity" << std::endl;
      m_CovarianceMomentaInverse.set_size(nbControlPoints * Dimension, nbControlPoints * Dimension);
      m_CovarianceMomentaInverse.set_identity();
    }
  }
}

template<class ScalarType, unsigned int Dimension>
bool
DeterministicAtlas<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<ScalarType>> &residuals) {
  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));

  const CrossSectionalDataSetType *const crossSectionalDataSet
      = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();

  return Superclass::ComputeResiduals(this->GetControlPoints(), momentas, target, residuals);
}

template<class ScalarType, unsigned int Dimension>
bool
DeterministicAtlas<ScalarType, Dimension>
::UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    VectorType &logLikelihoodTerms) {
  /// Initialization.
  logLikelihoodTerms.set_size(2);
  logLikelihoodTerms.fill(0.0);

  const MatrixType controlPoints = this->GetControlPoints();
  const std::vector<MatrixType> momentas = recast<MatrixType>(indRER.at("Momenta"));

  /// Data (residuals) term.
  std::vector<std::vector<ScalarType>> residuals;
  bool oob = this->ComputeResiduals(dataSet, popRER, indRER, residuals);
  const unsigned long nbSubjects = residuals.size();

  if (momentas.size() != nbSubjects)
    throw std::runtime_error("Number of momenta matrices and residuals mismatch in"
                                 " DeterministicAtlas::ComputeLogLikelihood");

  for (unsigned int s = 0; s < nbSubjects; s++) {
    for (unsigned int i = 0; i < this->m_NumberOfObjects; i++)
      logLikelihoodTerms[0] -= residuals[s][i] / m_DataSigmaSquared(i);
  }
  logLikelihoodTerms[0] *= 0.5;

  /// Regularity term
  if (m_UseRKHSNormForRegularization) // use the RKHS norm
  {
    KernelFactoryType *kfac = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> momKernelObj = kfac->CreateKernelObject(this->m_Def->GetKernelType());
    momKernelObj->SetKernelWidth(this->m_Def->GetKernelWidth());
    momKernelObj->SetSources(controlPoints);

    for (unsigned int s = 0; s < nbSubjects; s++) {
      momKernelObj->SetWeights(momentas[s]);
      MatrixType kMom = momKernelObj->Convolve(controlPoints);

      for (unsigned int i = 0; i < controlPoints.rows(); i++)
        logLikelihoodTerms[1] -= dot_product(kMom.get_row(i), momentas[s].get_row(i));
    }
  } else // covariance matrix given
  {
    for (unsigned int s = 0; s < nbSubjects; s++) {
      VectorType Moms = this->Vectorize(momentas[s]);
      logLikelihoodTerms[1] -= dot_product(Moms, this->GetCovarianceMomentaInverse() * Moms);
    }
  }
  logLikelihoodTerms[1] *= 0.5;

  return oob; // Out of box flag.
}

template<class ScalarType, unsigned int Dimension>
void
DeterministicAtlas<ScalarType, Dimension>
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

  if (momentas.size() != nbSubjects)
    throw std::runtime_error(
        "Number of Momentas and target multi-object mismatch in DeterministicAtlas::ComputeFunctionalGradient");

  MatrixType gradPos;
  MatrixListType gradTempL_L2, gradTempL_Sob;
  std::vector<MatrixType> gradMom(nbSubjects);

  this->UpdateDeformationAndKernelDataDomain(target);
  if (!this->m_FreezeTemplateFlag) {
    this->ComputeDataTermGradient(this->GetControlPoints(), momentas, m_DataSigmaSquared, target,
                                  gradPos, gradMom, gradTempL_L2);
    gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
  } else {
    this->ComputeDataTermGradient(this->GetControlPoints(), momentas, m_DataSigmaSquared, target,
                                  gradPos, gradMom);
  }
  AddGradientRegularityTerm(momentas, gradPos, gradMom);

  if (!this->m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  if (!this->m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }
  indGrad["Momenta"] = gradMom;
}

template<class ScalarType, unsigned int Dimension>
void
DeterministicAtlas<ScalarType, Dimension>
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

  if (momentas.size() != nbSubjects) {
    throw std::runtime_error(
        "Number of Momentas and target multi-object mismatch in DeterministicAtlas::ComputeFunctionalGradient");
  }

  MatrixType gradPos;
  MatrixListType gradTempL_L2, gradTempL_Sob;
  std::vector<MatrixType> gradMom(nbSubjects);

  this->UpdateDeformationAndKernelDataDomain(target);
  if (!this->m_FreezeTemplateFlag) {
    this->ComputeDataTermGradient(this->GetControlPoints(), momentas, m_DataSigmaSquared, target,
                                  gradPos, gradMom, gradTempL_L2);
    gradTempL_Sob = this->ConvolveGradTemplate(gradTempL_L2);
  } else {
    this->ComputeDataTermGradient(this->GetControlPoints(), momentas, m_DataSigmaSquared, target,
                                  gradPos, gradMom);
  }
  AddGradientRegularityTerm(momentas, gradPos, gradMom);

  if (!this->m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  if (!this->m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }
  indGrad["Momenta"] = gradMom;

  // Must be the population effects in alphabetical order first, the individual effects in alphabetical order then.
  gradSquaredNorms.set_size(0);
  if (!this->m_FreezeControlPointsFlag) { gradSquaredNorms.push_back(gradPos.sum_of_squares()); }
  if (!this->m_FreezeTemplateFlag) { gradSquaredNorms.push_back(dot_product(gradTempL_L2, gradTempL_Sob)); }
  gradSquaredNorms.push_back(MatrixListType(gradMom).sum_of_squares());
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected methods(s)
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
DeterministicAtlas<ScalarType, Dimension>
::AddGradientRegularityTerm(std::vector<MatrixType> const &momentas,
                            MatrixType &GradPos, std::vector<MatrixType> &GradMom) {
  const unsigned long numSubjects = momentas.size();

  // gradient of the regularity term
  if (m_UseRKHSNormForRegularization) // use the RKHS norm
  {
    const MatrixType controlPoints = this->GetControlPoints();

    typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
    typedef typename KernelFactoryType::KernelBaseType KernelType;
    KernelFactoryType *kfac = KernelFactoryType::Instantiate();

    std::shared_ptr<KernelType> momKernelObj = kfac->CreateKernelObject(this->m_Def->GetKernelType());
    momKernelObj->SetKernelWidth(this->m_Def->GetKernelWidth());
    momKernelObj->SetSources(controlPoints);

    for (unsigned int s = 0; s < numSubjects; s++) {
      momKernelObj->SetWeights(momentas[s]);
      GradMom[s] -= momKernelObj->Convolve(controlPoints);

      GradPos -= momKernelObj->ConvolveGradient(controlPoints, momentas[s]);
    }
  } else // covariance matrix given
  {
    const MatrixType covarianceMomentaInverse = this->GetCovarianceMomentaInverse();

    for (unsigned int s = 0; s < numSubjects; s++) {
      VectorType aux = covarianceMomentaInverse * this->Vectorize(momentas[s]);
      GradMom[s] -= this->VectorToMatrix(aux);
    }
  }
}

template class DeterministicAtlas<double,2>;
template class DeterministicAtlas<double,3>;
