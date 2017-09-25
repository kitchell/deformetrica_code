/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

/// Class file.
#include "AbstractAtlas.h"
#include "LinearAlgebra.h"

using namespace def::algebra;

/**
 *  \brief      Atlas object class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    A child class of the class Atlas for deterministic atlas estimation.
 */

template<class ScalarType, unsigned int Dimension>
class DeterministicAtlas : public AbstractAtlas<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Atlas Type.
  typedef AbstractAtlas<ScalarType, Dimension> Superclass;

  /// Longitudinal dataset type.
  typedef typename Superclass::LongitudinalDataSetType LongitudinalDataSetType;
  /// Cross-sectional data set type.
  typedef CrossSectionalDataSet<ScalarType, Dimension> CrossSectionalDataSetType;

  /// Multi-Deformable object type.
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;

  /// Kernel factory type.
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  /// Kernel type.
  typedef typename KernelFactoryType::KernelBaseType KernelType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  DeterministicAtlas();

  /// Copy constructor.
  DeterministicAtlas(const DeterministicAtlas &other);

  /// Makes a copy of the object.
  std::shared_ptr<DeterministicAtlas> Clone() const { return std::static_pointer_cast<DeterministicAtlas>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<typename Superclass::Superclass> doClone() const {
    return std::static_pointer_cast<typename Superclass::Superclass>(std::make_shared<DeterministicAtlas>(*this));
  };
 public:

  /// Destructor.
  virtual ~DeterministicAtlas();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

//	/// Sets control points stored in the file \e fn.
//	void SetControlPoints(const std::string& fn);

  /// Returns the inverse of the covariance matrix of momentas.
  MatrixType GetCovarianceMomentaInverse() const { return m_CovarianceMomentaInverse; }
  /// Sets inverse of covariance matrix of momentas.
  void SetCovarianceMomentaInverse(MatrixType const &CM) { m_CovarianceMomentaInverse = CM; }
  /// Sets inverse of covariance matrix of momentas from file fn.
  void SetCovarianceMomentaInverse(std::string const &fn);

  /// Returns data sigma squared of each object.
  VectorType GetDataSigmaSquared() const { return m_DataSigmaSquared; }
  /// Returns data sigma squared of object \e i.
  ScalarType GetDataSigmaSquared(unsigned int const &i) const { return m_DataSigmaSquared[i]; }
  /// Sets data sigma squared of each object to \e d2.
  void SetDataSigmaSquared(const VectorType d2) { m_DataSigmaSquared = d2; };

  /// Use the norm of the RKHS for computing the regularity term. Otherwise, use the pre-defined matrix \e m_CovMom_Inverse
  void SetRKHSNormForRegularization() { m_UseRKHSNormForRegularization = true; }
  /// Use the  pre-defined matrix \e m_CovMom_Inverse for computing the regularity term
  void UnSetRKHSNormForRegularization() { m_UseRKHSNormForRegularization = false; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates.
  virtual void Update();

  /// Saves the model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const {
    Superclass::Write(recast<MatrixType>(indRER.at("Momenta")), m_CovarianceMomentaInverse, m_DataSigmaSquared, dataSet);
  }

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                LinearVariableMapType const &popRER,
                                LinearVariablesMapType const &indRER,
                                std::vector<std::vector<ScalarType>> &residuals);

  /// Computes the log-likelihood first mode, given an input random effects realization.
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER) {
    VectorType decomposition;
    return UpdateFixedEffectsAndComputeCompleteLogLikelihood(dataSet, popRER, indRER, decomposition);
  }
  /// Updates the fixed effects and computes the likelihood at once. Allows to avoid repeted computations.
  virtual bool UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                                 const LinearVariableMapType &popRER,
                                                                 const LinearVariablesMapType &indRER,
                                                                 VectorType &logLikelihoodTerms);

  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions) {
    std::cout << "Warning : DeterministicAtlas::ComputeModelLogLikelihood method does not detail the contributions "
        "of each subject yet." << std::endl;
    return ComputeCompleteLogLikelihood(dataSet, popRER, indRER);
  }
  /// Computes the log-likelihood mode, ignoring contributions of priors and subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i) {
    std::cerr << "Exception : DeterministicAtlas::ComputeLogLikelihoodWithoutPriorsForSubject method "
        "is not available yet." << std::endl;
    return 0;
  }

  /// Computes the gradient of the log-likelihood first mode.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad);
  /// Computes the gradient of the log-likelihood first mode and the corresponding norms.
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
                                                            LinearVariablesMapType &indGrad) {
    std::cerr << "Exception : DeterministicAtlas::ComputeRandomEffectsMarginalLogLikelihoodGradient"
        " method is not available yet." << std::endl;
  }

  /// Computes the sufficient statistics of the deterministic atlas model. TODO.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics) {
    std::cerr << "Exception : the sufficient statistics are not available yet for the DeterministicAtlas model."
              << std::endl;
  }

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood. TODO.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics) {
    std::cerr << "Exception : the method DeterministicAtlas::UpdateFixedEffects is not available yet."
              << std::endl;
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected methods(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// ComputeFunctionalGradient auxiliary function.
  void AddGradientRegularityTerm(std::vector<MatrixType> const &momentas,
                                 MatrixType &GradPos, std::vector<MatrixType> &GradMom);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Optional covariance of the momenta inverse matrix, for regularization.
  MatrixType m_CovarianceMomentaInverse;

  /// Data sigma squared for each object.
  VectorType m_DataSigmaSquared;

  /// Flag to select the method to compute the regularity term.
  bool m_UseRKHSNormForRegularization;

}; /* class DeterministicAtlas */

