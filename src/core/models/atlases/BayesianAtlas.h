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

/// Core files.
#include "ParametricImage.h"
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"

/// Support files.
#include "ExactKernel.h"

using namespace def::algebra;
using namespace def::proba;

/**
 *  \brief      Atlas object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    A child class of the class Atlas for Bayesian atlas estimation.
 */
template<class ScalarType, unsigned int Dimension>
class BayesianAtlas : public AbstractAtlas<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Atlas Type.
  typedef AbstractAtlas<ScalarType, Dimension> Superclass;

  /// Longitudinal dataset type.
  typedef typename Superclass::LongitudinalDataSetType LongitudinalDataSetType;
  /// Cross-sectional data set type.
  typedef typename Superclass::CrossSectionalDataSetType CrossSectionalDataSetType;

  /// Parametric image type.
  typedef ParametricImage<ScalarType, Dimension> ParametricImageType;

  /// Multi-Deformable object type.
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;

  /// Kernel factory type.
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  /// Kernel type.
  typedef typename KernelFactoryType::KernelBaseType KernelType;

  /// Diffeos type.
  typedef Diffeos<ScalarType, Dimension> DiffeosType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  BayesianAtlas();

  /// Copy constructor.
  BayesianAtlas(const BayesianAtlas &other);

  /// Makes a copy of the object.
  std::shared_ptr<BayesianAtlas> Clone() const { return std::static_pointer_cast<BayesianAtlas>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<typename Superclass::Superclass> doClone() const {
    return std::static_pointer_cast<typename Superclass::Superclass>(std::make_shared<BayesianAtlas>(*this));
  };
 public:

  virtual ~BayesianAtlas();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// A prior on the control points will be used.
  void SetUsePriorOnControlPoints() { m_UsePriorOnControlPoints = true; }

  /// Returns a bool indicating if random control points are used.
  bool UseRandomControlPoints() { return m_UseRandomControlPoints; }
  /// Random controls points will be used.
  void SetUseRandomControlPoints() { m_UseRandomControlPoints = true; }

  /// Sets the control points fixed effect.
  virtual void SetControlPoints(MatrixType const &cp) {
    Superclass::SetControlPoints(cp);
    if (m_UseRandomControlPoints) { SetControlPointsRandomEffectMean(cp); }
  }
  /// Sets the control points fixed effect.
  void SetControlPoints(LinearVariableType const &cp) { SetControlPoints(recast<MatrixType>(cp)); }
//	/// Sets control points stored in the file \e fn.
//	void SetControlPoints(const std::string& fn);

  /// Gets the inverse of the covariance matrix of momentas.
  MatrixType GetCovarianceMomentaInverse() const {
    return recast<MatrixType>(this->m_FixedEffects.at("CovarianceMomentaInverse"));
  }
  /// Sets inverse of covariance matrix of momentas.
  void SetCovarianceMomentaInverse(MatrixType const &cmi) {
    this->m_FixedEffects["CovarianceMomentaInverse"] = cmi;
    SetMomentaRandomEffectCovarianceInverse(cmi);
  }

  /// Gets the data sigma squared of each object.
  VectorType GetDataSigmaSquared() const {
    return recast<VectorType>(this->m_FixedEffects.at("NoiseVariance"));
  }
  /// Returns data sigma squared of object \e i.
  ScalarType GetDataSigmaSquared(unsigned int const &i) const {
    return recast<VectorType>(this->m_FixedEffects["NoiseVariance"])[i];
  }
  /// Sets data sigma squared of each object to \e d2.
  void SetDataSigmaSquared(VectorType const &d2) {
    this->m_FixedEffects["NoiseVariance"] = d2;
  }

  /// Gets the mean of the momenta individual random effect.
  VectorType GetMomentaRandomEffectMean() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_IndividualRandomEffects.at("Momenta"))->GetMean();
  }
  /// Sets the mean of the momenta individual random effect.
  void SetMomentaRandomEffectMean(VectorType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_IndividualRandomEffects.at("Momenta"))->SetMean(m);
  }

  /// Gets the covariance of the momenta individual random effect.
  MatrixType GetMomentaRandomEffectCovariance() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_IndividualRandomEffects.at("Momenta"))->GetCovariance();
  }
  /// Gets the inverse of the covariance of the momenta individual random effect.
  MatrixType GetMomentaRandomEffectCovarianceInverse() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_IndividualRandomEffects.at("Momenta"))->GetCovarianceInverse();
  }
  /// Sets the covariance of the momenta individual random effect.
  void SetMomentaRandomEffectCovariance(MatrixType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_IndividualRandomEffects.at("Momenta"))->SetCovariance(m);
  }
  /// Sets the inverse of the covariance of the momenta individual random effect.
  void SetMomentaRandomEffectCovarianceInverse(MatrixType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_IndividualRandomEffects.at("Momenta"))->SetCovarianceInverse(m);
  }

  /// Gets the mean of the control points population random effect.
  VectorType GetControlPointsRandomEffectMean() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->GetMean();
  }
  /// Sets the mean of the control points population random effect.
  void SetControlPointsRandomEffectMean(VectorType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->SetMean(m);
  }
  /// Sets the mean of the control points population random effect, from a matrix.
  void SetControlPointsRandomEffectMean(MatrixType const &m) {
    SetControlPointsRandomEffectMean(m.vectorize());
  }
  /// Sets the mean of the control points population random effect, from a linear variable.
  void SetControlPointsRandomEffectMean(LinearVariableType const &m) {
    SetControlPointsRandomEffectMean(m.vectorize());
  }

  /// Gets the sqrt of the covariance of the control points population random effect.
  MatrixType GetControlPointsRandomEffectCovarianceSqrt() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->GetCovarianceSqrt();
  }
  /// Gets the inverse of the covariance of the control points population random effect.
  MatrixType GetControlPointsRandomEffectCovarianceInverse() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->GetCovarianceInverse();
  }
  /// Sets the sqrt of the covariance of the control points population random effect.
  void SetControlPointsRandomEffectCovarianceSqrt(MatrixType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->SetCovarianceSqrt(m);
  }

  /// Gets the momenta covariance matrix hyperparameter.
  ScalarType GetCovarianceMomenta_HyperParameter() const {
    return std::static_pointer_cast<InverseWishartDistributionType>(
        this->m_Priors.at("CovarianceMomenta"))->GetDegreesOfFreedom();
  }
  /// Sets the momenta covariance matrix hyperparameter to \e d.
  void SetCovarianceMomenta_HyperParameter(ScalarType const &d) {
    std::static_pointer_cast<InverseWishartDistributionType>(
        this->m_Priors.at("CovarianceMomenta"))->SetDegreesOfFreedom(d);
  }

  /// Gets the prior on the momenta covariance.
  MatrixType GetCovarianceMomenta_Prior() const {
    return std::static_pointer_cast<InverseWishartDistributionType>(
        this->m_Priors.at("CovarianceMomenta"))->GetScaleMatrix();
  }
  /// Gets the inverse of the prior on the momenta covariance.
  MatrixType GetCovarianceMomenta_Prior_Inverse() const {
    return std::static_pointer_cast<InverseWishartDistributionType>(
        this->m_Priors.at("CovarianceMomenta"))->GetScaleMatrixInverse();
  }
  /// Sets the prior on the momenta covariance.
  void SetCovarianceMomenta_Prior_Inverse(MatrixType const &M) {
    std::static_pointer_cast<InverseWishartDistributionType>(
        this->m_Priors.at("CovarianceMomenta"))->SetScaleMatrixInverse(M);
  }
  /// Loads the prior on the momenta covariance from a file \e fn.
  void SetCovarianceMomenta_Prior_Inverse(std::string const &fn);

  /// Gets the hyperparameter on the data variance.
  VectorType GetDataSigmaSquared_HyperParameter() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->GetDegreesOfFreedom();
  }
  /// Sets the hyperparameter on the data variance.
  void SetDataSigmaSquared_HyperParameter(VectorType const &d) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->SetDegreesOfFreedom(d);
  }

  /// Gets the prior on the data variance.
  VectorType GetDataSigmaSquared_Prior() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->GetScaleVector();
  }
  /// Sets the prior on the data variance.
  void SetDataSigmaSquared_Prior(VectorType const &V) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->SetScaleVector(V);
  }

  /// Gets the mean of the prior distribution on the parametric template.
  VectorType GetParametricTemplatePriorMean() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ParametricTemplate"))->GetMean();
  }
  /// Sets the mean of the prior distribution on the parametric template.
  void SetParametricTemplatePriorMean(VectorType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ParametricTemplate"))->SetMean(m);
  }

  /// Gets the covariance of the prior distribution on the parametric template.
  MatrixType GetParametricTemplatePriorCovarianceInverse() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ParametricTemplate"))->GetCovarianceInverse();
  }
  /// Sets the covariance of the prior distribution on the parametric template.
  void SetParametricTemplatePriorCovarianceInverse(MatrixType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ParametricTemplate"))->SetCovarianceInverse(m);
  }

  /// Gets the mean of the prior distribution on the random control points.
  VectorType GetControlPointsPriorMean() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->GetMean();
  }
  /// Sets the mean of the prior distribution on the random control points.
  void SetControlPointsPriorMean(VectorType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->SetMean(m);
  }

  /// Gets the sqrt of the covariance of the prior distribution on the random control points.
  MatrixType GetControlPointsPriorCovarianceSqrt() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->GetCovarianceSqrt();
  }
  /// Gets the inverse of the covariance of the prior distribution on the random control points.
  MatrixType GetControlPointsPriorCovarianceInverse() const {
    return std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->GetCovarianceInverse();
  }
  /// Sets the covariance of the prior distribution on the random control points.
  void SetControlPointsPriorCovarianceSqrt(MatrixType const &m) {
    std::static_pointer_cast<NormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->SetCovarianceSqrt(m);
  }

  /// Sets the dimension of the noise.
  void SetNoiseDimension(std::vector<unsigned long> V) { m_NoiseDimension = V; }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initializes.
  virtual void Update();

  /// Saves the model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const {
    Superclass::Write(recast<MatrixType>(indRER.at("Momenta")), GetCovarianceMomentaInverse(), GetDataSigmaSquared(), dataSet);
  }

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
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

  /// Computes the gradient of the log-likelihood mode.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad);
  /// Computes the gradient of the log-likelihood mode and the corresponding norms.
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

  /// Computes the sufficient statistics of the bayesian atlas model.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics);
  /// Computes the sufficient statistics for a given subject.
  void ComputeSufficientStatisticsSubject(const MatrixType &momenta,
                                          const std::shared_ptr<const DeformableMultiObjectType> target,
                                          LinearVariableMapType &ss_i) const;

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Compute the optimal covariance matrix of momentas, and associated log-likelihood terms.
  void UpdateCovMomInverseAndComputeTerms(const std::vector<MatrixType> &momentas,
                                          ScalarType &covMomTerm,
                                          ScalarType &covMomPriorTerm);
  /// Compute the optimal data sigma squared for each object, and associated log-likelihood terms.
  void UpdateDataSigmaSquaredAndComputeTerms(const std::vector<std::vector<ScalarType>> &residuals,
                                             ScalarType &dataTerm,
                                             ScalarType &dataPriorTerm);

  /// Initializes the prior on the momenta covariance.
  void InitializeCovMomPrior();
  /// Initializes the prior on the parametric template.
  void InitializeParametricTemplatePrior();
  /// Initializes the prior on the random control points.
  void InitializeControlPointsPrior();
  /// Initializes the momenta individual random effect.
  void InitializeMomentaRandomEffect();
  /// Initializes the control points population random effect.
  void InitializeControlPointsRandomEffect();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Dimension of the noise.
  std::vector<unsigned long> m_NoiseDimension;
  /// Flag that indicates if the template is parametric. A prior is then defined on the photometric weights.
  bool m_UseParametricTemplate;
  /// Flag that indicates if a prior should be used on the control points fixed effect.
  bool m_UsePriorOnControlPoints;
  /// Flag that indicates if the control points are considered as a random effect.
  bool m_UseRandomControlPoints;

}; /* class BayesianAtlas */

