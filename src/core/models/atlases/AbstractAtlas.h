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
#include "src/core/observations/data_sets/CrossSectionalDataSet.h"
#include "Diffeos.h"
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"
#include "AbstractStatisticalModel.h"

/// Support files.
#include "KernelFactory.h"

/// Input-output files.
#include "DeformationFieldIO.h"
#include "MatrixDLM.h"

/// Librairies files.
#include "itkImage.h"
#include "itkSimpleFastMutexLock.h"

using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      AbstractAtlas object class.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    An atlas is the combination of a template shape, control points, and parameters of variability
 *	            such as covariance matrices of the momentum vectors.
 */
template<class ScalarType, unsigned int Dimension>
class AbstractAtlas : public AbstractStatisticalModel<ScalarType, Dimension> {
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

  /// Multi-deformable object type.
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;

  /// ITK image type.
  typedef itk::Image<ScalarType, Dimension> ImageType;
  /// ITK image pointer type.
  typedef typename ImageType::Pointer ImageTypePointer;
  /// ITK image index type.
  typedef typename ImageType::IndexType ImageIndexType;
  /// ITK image point type.
  typedef typename ImageType::PointType ImagePointType;
  /// ITK image spacing type.
  typedef typename ImageType::SpacingType ImageSpacingType;

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
  AbstractAtlas();

  /// Copy constructor.
  AbstractAtlas(const AbstractAtlas &other);

  /// Makes a copy of the object.
  std::shared_ptr<AbstractAtlas> Clone() const { return std::static_pointer_cast<AbstractAtlas>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<Superclass> doClone() const = 0;
 public:

  /// Destructor.
  virtual ~AbstractAtlas();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the fixed effects.
  virtual void SetFixedEffects(const LinearVariableMapType &map);

  /// Returns the template deformable objects.
  std::shared_ptr<DeformableMultiObjectType> GetTemplate() const { return m_Template; }
  /// Sets the template deformable objects to \e objects.
  void SetTemplate(std::shared_ptr<DeformableMultiObjectType> const temp) {
    m_Template = temp;
    Superclass::m_FixedEffects["TemplateData"] = temp->GetImageIntensityAndLandmarkPointCoordinates();
    ///TODO : why not do a m_def->SetDeformableMultiObject here ?
  }

  /// Returns image intensity and landmark point coordinates of the template.
  MatrixListType GetTemplateData() const {
    return recast<MatrixListType>(Superclass::m_FixedEffects.at("TemplateData"));
  }
  /// Sets image intensity and landmark point coordinates of the template.
  void SetTemplateData(const MatrixListType &tempData) {
    Superclass::m_FixedEffects["TemplateData"] = tempData;
    m_Template->UpdateImageIntensityAndLandmarkPointCoordinates(tempData);
    m_Template->Update();
  }

  /// Returns the name of template objects.
  std::vector<std::string> GetTemplateObjectsName() const { return m_TemplateObjectsName; }
  /// Sets template objects name for saving.
  void SetTemplateObjectsName(const std::vector<std::string> &ObjectsName) { m_TemplateObjectsName = ObjectsName; }

  /// Returns the extension of template objects name.
  std::vector<std::string> GetTemplateObjectsNameExtension() const { return m_TemplateObjectsNameExtension; }
  /// Sets template objects name extension for saving.
  void SetTemplateObjectsNameExtension(const std::vector<std::string> &ObjectsNameExt) {
    m_TemplateObjectsNameExtension = ObjectsNameExt;
  }

  /// Returns a tight bounding box including template objects and control points.
  MatrixType GetBoundingBox() const { return m_BoundingBox; }

  /// Sets the freeze template flag.
  void SetFreezeTemplateFlag(bool const &flag) { m_FreezeTemplateFlag = flag; }
  /// Sets the freeze control points flag.
  void SetFreezeControlPointsFlag(bool const &flag) { m_FreezeControlPointsFlag = flag; }

  /// Returns control points positions.
  MatrixType GetControlPoints() const {
    return recast<MatrixType>(Superclass::m_FixedEffects.at("ControlPoints"));
  }
  /// Sets the control points.
  virtual void SetControlPoints(MatrixType const &CP) { Superclass::m_FixedEffects["ControlPoints"] = CP; }
  /// Sets control points stored in the file \e fn.
  void SetControlPoints(const std::string &fn);

  /// Returns the number of objects in the template.
  int GetNumberOfObjects() const { return m_NumberOfObjects; }

  /// Gets control point spacing.
  ScalarType GetCPSpacing() const { return m_CPSpacing; }
  /// Sets control point spacing. Used in case no control points have been set to define a regular lattice of control points.
  void SetCPSpacing(const ScalarType s) { m_CPSpacing = s; }

  /// Sets Diffeos to deform the template.
  void SetDiffeos(std::shared_ptr<DiffeosType> const def) { m_Def = def; }

  /// Sets the size of the smoothing kernel to \e d. See Atlas::ConvolveGradTemplate().
  void SetSmoothingKernelWidth(const ScalarType d) { m_SmoothingKernelWidth = d; }

  /// Returns the number of threads.
  unsigned int GetNumberOfThreads() const { return m_NumberOfThreads; }
  /// Sets the number of threads to \e n.
  void SetNumberOfThreads(const unsigned int n) { m_NumberOfThreads = n; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates the bounding box and initializes the control points if needed.
  virtual void Update();

  /// Saves the model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const = 0;

  /// Samples from the model, and stores the used parameters in \e popRER and \e indRER.
  virtual bool Sample(LongitudinalDataSetType *realizations,
                      LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER) const {
    std::cerr << "Error : AbstractAtlas::Sample method is not available yet." << std::endl;
    return true;
  }

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                LinearVariableMapType const &popRER,
                                LinearVariablesMapType const &indRER,
                                std::vector<std::vector<std::vector<ScalarType>>> &residuals);
  /// Variation with a residuals structure adapted to cross-sectional data.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                LinearVariableMapType const &popRER,
                                LinearVariablesMapType const &indRER,
                                std::vector<std::vector<ScalarType>> &residuals) = 0;
  /// Variation of the ComputeResiduals method with natural parameters.
  virtual bool ComputeResiduals(const MatrixType &controlPoints,
                                const std::vector<MatrixType> &momentas,
                                const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                                std::vector<std::vector<ScalarType>> &residuals);
  /// Computes the residuals for a given subject.
  bool ComputeResidualsSubject(const MatrixType &controlPoints,
                               const MatrixType &momenta,
                               const std::shared_ptr<DeformableMultiObjectType> target,
                               std::vector<ScalarType> &residuals) const;

  /// Computes the complete log-likelihood, given an input random effects realization ("RER").
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER) = 0;
  /// Updates the fixed effects which have a closed-form update and computes the complete likelihood at once.
  virtual bool UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                                 const LinearVariableMapType &popRER,
                                                                 const LinearVariablesMapType &indRER,
                                                                 VectorType &logLikelihoodTerms) = 0;
  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions) = 0;
  /// Computes the model log-likelihood, ignoring the contribution of subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i) = 0;

  /// Computes the gradient of the log-likelihood mode. ("RER" Random Effects Realization).
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad) = 0;
  /// Computes the gradient of the log-likelihood mode and the corresponding squared norms.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad,
                                                    VectorType &gradSquaredNorms) = 0;
  /// Computes the marginal log-likelihood gradient with respect to the random effects of the model.
  virtual void ComputeLogLikelihoodGradientWrtRandomEffects(const LongitudinalDataSetType *const dataSet,
                                                            const LinearVariableMapType &popRER,
                                                            const LinearVariablesMapType &indRER,
                                                            LinearVariableMapType &popGrad,
                                                            LinearVariablesMapType &indGrad) = 0;

  /// Computes the sufficient statistics of the atlas model.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics) = 0;

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics) = 0;

  /// If no control points have been set, this method generates a regular lattice of control points.
  void InitializeControlPoints(bool optimize = false);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Saves the model.
  void Write(std::vector<MatrixType> const &momentas,
             MatrixType const &covarianceMomentaInverse,
             VectorType const &dataSigmaSquared,
             const LongitudinalDataSetType *const dataSet) const;
  /// Saves the deformation of the template to every subject, and the residuals.
  void WriteAtlasToSubjectDeformations(std::vector<MatrixType> const &momentas, const LongitudinalDataSetType *const dataSet) const;
  /// Saves atlas parameters (such as control points, initial momentas, covariance matrices, etc).
  void WriteAtlasParameters(std::vector<MatrixType> const &momentas, MatrixType const &covarianceMomentaInverse,
                            VectorType const &dataSigmaSquared) const;
  /// Saves the template.
  void WriteTemplateData() const;

  /// Computes the gradient of the sum of subject's data terms with respect to momenta, control points and template objects.
  // The data term is the sum of the residuals normalized by 2 times the data sigma squared.
  void ComputeDataTermGradient(const MatrixType &controlPoints,
                               const std::vector<MatrixType> &momentas,
                               const VectorType &dataSigmaSquared,
                               const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                               MatrixType &gradPos, std::vector<MatrixType> &gradMom, MatrixListType &gradTempL_L2);
  /// Computes the gradient of the data term with respect to momentas, control points and template objects.
  void ComputeDataTermGradientSubject(const MatrixType &controlPoints,
                                      const MatrixType &momenta,
                                      const VectorType &dataSigmaSquared,
                                      const std::shared_ptr<DeformableMultiObjectType> target,
                                      MatrixType &dPos, MatrixType &dMom, MatrixListType &dTempL);

  /// Computes the gradient of the sum of subject's data terms with respect to the control points and the momentas.
  void ComputeDataTermGradient(const MatrixType &controlPoints,
                               const std::vector<MatrixType> &momentas,
                               const VectorType &dataSigmaSquared,
                               const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                               MatrixType &gradPos, std::vector<MatrixType> &gradMom);
  /// Computes the gradient of the data term with respect to the control points and the momentas.
  void ComputeDataTermGradientSubject(const MatrixType &controlPoints,
                                      const MatrixType &momenta,
                                      const VectorType &dataSigmaSquared,
                                      const std::shared_ptr<DeformableMultiObjectType> target,
                                      MatrixType &dPos, MatrixType &dMom);

  /// Computes the Sobolev gradient of template objects of landmark type (or children types) from the \f$L^2\f$ gradient
  MatrixListType ConvolveGradTemplate(MatrixListType &GradTemplate_L2);

  /// Updates the data domain of the deformation and the kernel.
  void UpdateDeformationAndKernelDataDomain(const std::vector<std::shared_ptr<DeformableMultiObjectType>> target);

  /// Converts a matrix of size N x Dimension to a vector \e V of length Dimension x N.
  VectorType Vectorize(const MatrixType &M) const { return M.vectorise_row_wise(); }
  /// Converts a vector \e V of length Dimension x N to a matrix of size N x Dimension. Inverse operation of Vectorize.
  MatrixType VectorToMatrix(const VectorType &V) const {
    int nrow = V.size() / Dimension;
    return V.convert_to_matrix_row_wise(nrow, Dimension);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Template object.
  std::shared_ptr<DeformableMultiObjectType> m_Template;
  /// Template objects name.
  std::vector<std::string> m_TemplateObjectsName;
  /// Extension of template objects name.
  std::vector<std::string> m_TemplateObjectsNameExtension;
  /// Template bounding box: the union of the template bounding box and the bounding box around control points.
  MatrixType m_BoundingBox;

  /// Number of objects in the template.
  unsigned int m_NumberOfObjects;
  /// CP spacing used if no set of control points is given.
  ScalarType m_CPSpacing;

  /// Flag which freezes the template when true.
  bool m_FreezeTemplateFlag;
  /// Flag which freezes the control points when true.
  bool m_FreezeControlPointsFlag;

  /// Working deformation to compute atlas deformation.
  std::shared_ptr<DiffeosType> m_Def;

  /// Kernel width to compute the Sobolev gradient of the log-likelihood w.r.t. template variable.
  ScalarType m_SmoothingKernelWidth;

  /// Number of threads.
  unsigned int m_NumberOfThreads;

 protected :

  /// \cond HIDE_FOR_DOXYGEN

  static ITK_THREAD_RETURN_TYPE _residualsThread(void *arg);

  static ITK_THREAD_RETURN_TYPE _gradientResidualsThread(void *arg);

  static ITK_THREAD_RETURN_TYPE _gradientWithoutTemplateResidualsThread(void *arg);

  itk::SimpleFastMutexLock m_Mutex;

  unsigned int m_MT_SubjectCounter;
  unsigned int m_MT_NumberOfSubjects;

  bool m_MT_OutOfBox;

  MatrixType m_MT_ControlPoints;
  std::vector<MatrixType> m_MT_Momentas;
  VectorType m_MT_DataSigmaSquared;
  std::vector<std::shared_ptr<DeformableMultiObjectType>> m_MT_Target;
  std::vector<std::vector<ScalarType>> m_MT_Residuals;

  MatrixType m_MT_GradPos;
  std::vector<MatrixType> m_MT_GradMom;
  MatrixListType m_MT_GradTempL_L2;

  /// \endcond


}; /* class AbstractAtlas */


