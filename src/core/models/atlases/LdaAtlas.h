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
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"

/// Support files.
#include "ExactKernel.h"
#include <boost/filesystem.hpp>


using namespace def::algebra;
using namespace def::proba;

/**
 *  \brief      Atlas object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    A child class of the class Atlas for Bayesian Non Linear Discriminant Analysis
 */
template<class ScalarType, unsigned int Dimension>
class LdaAtlas : public AbstractAtlas<ScalarType, Dimension> {
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
  LdaAtlas();

  /// Copy constructor.
  LdaAtlas(const LdaAtlas &other);

  /// Makes a copy of the object.
  std::shared_ptr<LdaAtlas> Clone() const { return std::static_pointer_cast<LdaAtlas>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<typename Superclass::Superclass> doClone() const {
    return std::static_pointer_cast<typename Superclass::Superclass>(std::make_shared<LdaAtlas>(*this));
  };
 public:

  virtual ~LdaAtlas();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the control points fixed effect.
  virtual void SetControlPoints(MatrixType const &cp) {
    Superclass::SetControlPoints(cp);
  }

  /// Sets the control points fixed effect.
  void SetControlPoints(LinearVariableType const &cp) { SetControlPoints(recast<MatrixType>(cp)); }


  /// Gets the data sigma squared of each object.
  VectorType GetNoiseVariance() const {
    return recast<VectorType>(this->m_FixedEffects.at("NoiseVariance"));
  }

  /// Returns data sigma squared of object \e i.
  ScalarType GetNoiseVariance(unsigned int const &i) const {
    return recast<VectorType>(this->m_FixedEffects["NoiseVariance"])[i];
  }
  /// Sets data sigma squared of each object to \e d2.
  void SetNoiseVariance(VectorType const &d2) {
    this->m_FixedEffects["NoiseVariance"] = d2;
  }

  ///Sets the number of classes in the LDA
  void SetNumberOfClasses(unsigned int nbClasses){m_NumberOfClasses = nbClasses;}
  ///Gets the number of classes in the LDA
  unsigned int GetNumberOfClasses(){return m_NumberOfClasses;}


  ///Sets the number of classes in the LDA
  void SetNoiseDimension(std::vector<unsigned long> noiseDimension){m_NoiseDimension = noiseDimension;}

  void SetIntraClassPCADimension(unsigned int intraClassPCADimension){m_IntraClassPCADimension = intraClassPCADimension;}

  void SetNumberOfSubjectsPerClass(unsigned int numberOfSubjectsPerClass){m_NumberOfSubjectsPerClass = numberOfSubjectsPerClass;}

  void SetFreezeF(const bool b){
    m_FreezeF = b;
    if (m_FreezeF)
      std::cout << "F will be frozen during estimation" << std::endl;
  }
  void SetFreezeG(const bool b){
    m_FreezeG = b;
    if (m_FreezeG)
      std::cout << "G will be frozen during estimation" << std::endl;
  }

  ///Set the parameter gamma in the ARD distribution of G, both in the fixed effects and in the prior distribution used for regularization
  void SetGamma(VectorType const& newGamma){
    this->m_FixedEffects["gamma"] = newGamma;
    std::static_pointer_cast<AutomaticRelevanceDeterminationDistributionType>(this->m_Priors.at("G"))->SetGamma(newGamma);
  }

  VectorType GetGamma() const  {return std::static_pointer_cast<AutomaticRelevanceDeterminationDistributionType>(this->m_Priors.at("G"))->GetGamma(); }

  void SetG(MatrixType const& G){
    this->m_FixedEffects["G"] = G;
  }

  MatrixType GetG() const {
    return recast<MatrixType>(this->m_FixedEffects["G"]);
  }


  void SetF(MatrixType const& F){
    this->m_FixedEffects["F"] = F;
  }

  MatrixType GetF() const{
    return recast<MatrixType>(this->m_FixedEffects["F"]);
  }

  void SetAlpha(MatrixListType const& Alpha){
    this->m_FixedEffects["alpha"] = Alpha;
  }

  virtual void SetFixedEffects(const LinearVariableMapType &map);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initializes.
  virtual void Update();

  /// Saves the model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const;

  /// Saves the classes info of the moel.
  virtual void WriteClassesInfo(MatrixType const& G) const;

  /// Samples from the model, and stores the used parameters in \e popRER and \e indRER.
  virtual bool Sample(LongitudinalDataSetType *dataSet,
                      LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER) const;

  ///Method specific to this atlas, writing the pga components
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> GetComponents() const;

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

  /// Computes the model log-likelihood.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER) {
    return Superclass::Superclass::ComputeModelLogLikelihood(dataSet, popRER, indRER);
  }
  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions);
  /// Optional variation, where only one population random effect realization (RER) has been modified since last call.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               const std::string &modifiedVar,
                                               std::vector<ScalarType> &contributions) {
    return Superclass::Superclass::ComputeModelLogLikelihood(dataSet, popRER, indRER, modifiedVar, contributions);
  };

  /// Optional variation, featuring a temperature parameter.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               const ScalarType &temperature,
                                               const std::string &modifiedVar,
                                               std::vector<ScalarType> &contributions) {
    return Superclass::Superclass::ComputeModelLogLikelihood(
        dataSet, popRER, indRER, temperature, modifiedVar, contributions);
  };

  /// Computes the model log-likelihood, ignoring the contribution of subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i){  throw(std::runtime_error("ComputeModelLogLikelihoodForSubject Not implemented for LDA atlas"));};

  /// Optional variation, where only one individual random effect realization (RER) has been modified since last call.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i,
                                                         const std::string &modifiedVar) {
    return Superclass::Superclass::ComputeModelLogLikelihoodForSubject(dataSet, popRER, indRER, i, modifiedVar);
  };

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
                                                            LinearVariablesMapType &indGrad) {throw(std::runtime_error("ComputeLogLikelihoodGradientWrtRandomEffects not implemented for LdaAtlas")); }

  /// Computes the sufficient statistics of the bayesian atlas model.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics){throw(std::runtime_error("ComputeSufficientStatistics not implemented for LdaAtlas"));};

  /// Computes the sufficient statistics for a given subject.
  void ComputeSufficientStatisticsSubject(const MatrixType &momenta,
                                          const std::shared_ptr<const DeformableMultiObjectType> target,
                                          LinearVariableMapType &ss_i) const {throw(std::runtime_error("ComputeSufficientStatisticsSubject not implemented for LdaAtlas"));};

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics){throw(std::runtime_error("ComputeSufficientStatisticsSubject not implemented for LdaAtlas"));};

  /// Optional variation, featuring a temperature parameter.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics,
                                  const ScalarType &temperature) {
    throw(std::runtime_error("UpdateFixedEffects Not implemented for LDA atlas"));
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

///Compute the optimal gammas in the automatic relevance determination for the intraclasspca.
  void UpdateGamma(ScalarType& gammaTerm, unsigned int& nbSubjects);

  void UpdateNoiseVarianceAndComputeTerms(std::vector<std::vector<ScalarType>> const& residuals, ScalarType& dataTerm,
                                          ScalarType& GTerm, ScalarType& FTerm, ScalarType& alphaTerm, ScalarType& betaTerm,
                                          ScalarType& noiseVarianceTerm, const LinearVariablesMapType& indRER, std::vector<unsigned int> const& classes);

  ///Initialization of the prior distribution and of the initial value of F
  void InitializeF();
  ///Initialization of the prior distribution and of the initial value of G
  void InitializeG();
  ///Initialization of the prior distribution and of the initial value of alpha
  void InitializeAlpha();
  ///Initialization of the distribution beta (it is a regularization on beta)
  void InitializeBeta();
  ///Initialization of the diffeo bounding box
  void InitializeBoundingBox();



  ///Generate the momenta for each subject m=F alpha + G beta
  ///Note that classes is std::vector<unsigned int> which contains the number of the class for each subject
  ///Note that alpha is a matrix list type of length the number of classes. Each matrix is of length number of classes -1 * dimension.
  std::vector<MatrixType> GetMomentasFromLdaDescription(std::vector<MatrixType> const& betas, std::vector<unsigned int> classes) const;

  ///Returns the matrix gradMom * transpose(alpha)
  MatrixType VectorizeAndMultiply(MatrixType const& gradMom, MatrixType const& alpha);

  ///Return transpose(F) * gradMom, where gradMom is a matrix representing the coordinates of the momenta (needs to be flattened)
  MatrixType MultiplyMatrixTransposeByVector(MatrixType const& F, MatrixType const& gradMom);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Number of principal components in the intraclass reduction.
  unsigned int m_IntraClassPCADimension;
  /// Number of classes in the linear discriminant analysis
  unsigned int m_NumberOfClasses;
  ///Number of subject per class (used to sample from the model)
  unsigned int m_NumberOfSubjectsPerClass;
  /// Noise dimension for objects
  std::vector<unsigned long> m_NoiseDimension;
  ///Whether learn G or just to fit the observations to the given G.
  bool m_FreezeG;
  ///Whether learn F or just to fit the observations to the given F.
  bool m_FreezeF;

}; /* class LdaAtlas */

