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

/// Core files.
#include "TimeSeriesDataSet.h"
#include "Diffeos.h"

/// Support files.
#include "NormalDistribution.h"
#include "LinearAlgebra.h"

/// Input-output files.
#include "DeformationFieldIO.h"
#include "MatrixDLM.h"
#include "MatrixDLM.h"

/// Librairies files.
#include "itkImage.h"
#include "itkSimpleFastMutexLock.h"

using namespace def::algebra;

/**
 *	\brief      Regression object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A Regression is the combination of a template shape along with control points/initial momenta
 *		    	that deform the template shape to match closely time-indexed shape observations.
 */
template<class ScalarType, unsigned int Dimension>
class Regression : public AbstractStatisticalModel<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Statistical model type.
  typedef AbstractStatisticalModel<ScalarType, Dimension> Superclass;

  /// Longitudinal dataset type.
  typedef typename Superclass::LongitudinalDataSetType LongitudinalDataSetType;
  /// Time series data set type.
  typedef TimeSeriesDataSet<ScalarType, Dimension> TimeSeriesDataSetType;

  /// Multi-deformable object type.
  typedef typename Superclass::DeformableMultiObjectType DeformableMultiObjectType;
  /// Deformable object type.
  typedef typename DeformableMultiObjectType::AbstractGeometryType AbstractGeometryType;
  /// List of deformable objects type.
  typedef typename DeformableMultiObjectType::AbstractGeometryList AbstractGeometryList;

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

  /// Diffeos type.
  typedef Diffeos<ScalarType, Dimension> DiffeosType;

  /// Kernel factory type.
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  /// Kernel type.
  typedef typename KernelFactoryType::KernelBaseType KernelType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  Regression();

  /// Copy constructor.
  Regression(const Regression &other);

  /// Makes a copy of the object.
  std::shared_ptr<Regression> Clone() const { return std::static_pointer_cast<Regression>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<Superclass> doClone() const {
    return std::static_pointer_cast<Superclass>(std::make_shared<Regression>(*this));
  };
 public:

  /// Destructor
  ~Regression();


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
  }

  /// Returns image intensity and landmark point coordinates of the template.
  MatrixListType GetTemplateData() const {
    return recast<MatrixListType>(Superclass::m_FixedEffects["TemplateData"]);
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

  /// Returns a tight bounding box including baseline shape objects and control points.
  MatrixType GetBoundingBox() { return m_BoundingBox; }

  /// Sets the freeze template flag.
  void SetFreezeTemplateFlag(bool const &flag) { m_FreezeTemplateFlag = flag; }
  /// Sets the freeze control points flag.
  void SetFreezeControlPointsFlag(bool const &flag) { m_FreezeControlPointsFlag = flag; }

  /// Returns control points positions.
  MatrixType GetControlPoints() const {
    return recast<MatrixType>(Superclass::m_FixedEffects["ControlPoints"]);
  }
  /// Updates positions of control points
  void SetControlPoints(MatrixType const &CP) { Superclass::m_FixedEffects["ControlPoints"] = CP; }
  /// Sets control points stored in the file \e fn
  void SetControlPoints(const std::string &fn);

  /// Returns the number of objects in the template.
  int GetNumberOfObjects() const { return m_NumberOfObjects; }

  /// Gets control point spacing.
  ScalarType GetCPSpacing() const { return m_CPSpacing; }
  /// Sets control point spacing. Used in case no control points have been set to define a regular lattice of control points.
  void SetCPSpacing(const ScalarType s) { m_CPSpacing = s; }

  /// Returns the initial momenta.
  MatrixType GetInitialMomenta() const {
    return recast<MatrixType>(Superclass::m_FixedEffects["InitialMomenta"]);
  }
  /// Sets the initial momenta to \e M.
  void SetInitialMomenta(MatrixType const &M) { Superclass::m_FixedEffects["InitialMomenta"] = M; }
  /// Sets the initial momenta by reading the file \e fn.
  void SetInitialMomenta(std::string const &fn);

  /// Returns data sigma squared of each object.
  VectorType GetDataSigmaSquared() const { return m_DataSigmaSquared; }
  /// Returns data sigma squared of object \e i.
  ScalarType GetDataSigmaSquared(int i) const { return m_DataSigmaSquared[i]; }
  /// Sets data sigma squared of each object to \e d2.
  void SetDataSigmaSquared(const VectorType d2) { m_DataSigmaSquared = d2; };
  /// Sets data sigma squared of object \e i to \e d2.
  void SetDataSigmaSquared(const ScalarType d2, const int i) { m_DataSigmaSquared[i] = d2; };

  /// Sets diffeos to deform the template.
  void SetDiffeos(std::shared_ptr<DiffeosType> def) { m_Def = def; }

  /// Sets the size of the smoothing kernel to \e d. See Regression::ConvolveGradTemplate().
  void SetSmoothingKernelWidth(ScalarType d) { m_SmoothingKernelWidth = d; }

  /// Sets the write-full-trajectories flag.
  void SetWriteFullTrajectoriesFlag(bool const &flag) { m_WriteFullTrajectoriesFlag = flag; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates the bounding box and initializes the control points if needed.
  virtual void Update();

  /// Saves the regression model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const;

  /// Samples from the model, and stores the used parameters in \e popRER and \e indRER.
  virtual bool Sample(LongitudinalDataSetType *realizations,
                      LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER) const {
    std::cerr << "Error : Regression::Sample method is not available yet." << std::endl;
    return true;
  }

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                const LinearVariableMapType &popRER,
                                const LinearVariablesMapType &indRER,
                                std::vector<std::vector<std::vector<ScalarType>>> &residuals);
  /// Variation with a residuals structure adapted to time series data.
  bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                        LinearVariableMapType const &popRER,
                        LinearVariablesMapType const &indRER,
                        std::vector<std::vector<ScalarType>> &residuals);

  /// Evaluates how well the regression model explains a specific time series dataset.
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER) {
    VectorType decomposition;
    return UpdateFixedEffectsAndComputeCompleteLogLikelihood(dataSet, popRER, indRER, decomposition);
  }
  /// Computes the functional first mode, given an input random effects realization.
  virtual bool UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                                 const LinearVariableMapType &popRER,
                                                                 const LinearVariablesMapType &indRER,
                                                                 VectorType &logLikelihoodTerms);

  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions);
  /// Computes the log-likelihood mode, ignoring contributions of priors and subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i) {
    std::cerr << "Exception : Regression::ComputeLogLikelihoodWithoutPriorsForSubject method "
        "makes no sense and should not be called." << std::endl;
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
    std::cerr << "Exception : Regression::ComputeRandomEffectsMarginalLogLikelihoodGradient"
        " method is not available yet." << std::endl;
  }

  /// Computes the sufficient statistics of the regression model.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics);

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics);

  /// If no control points have been set, this method generates a regular lattice of control points.
  void InitializeControlPoints(bool optimize = false);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Saves the current template shape(s).
  void WriteTemplateFlow() const;
  /// Saves regression parameters.
  void WriteRegressionParameters() const;
  /// Saves the current template.
  void WriteTemplateData() const;

  /// Computes the gradient of the data term.
  void ComputeDataTermGradient(const MatrixType &momentas,
                               const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                               std::vector<unsigned int> timeIndices,
                               MatrixType &gradPos,
                               MatrixType &gradMom,
                               MatrixListType &gradTempL_L2,
                               MatrixListType &gradTempL_Sob);
  /// Computes the gradient of regularity portion of the functional
  void AddGradientRegularityTerm(MatrixType gradPos, MatrixType gradMom, const MatrixType &initialMomenta);

  /// Computes the Sobolev gradient of template objects of landmark type (or children types) from the \f$L^2\f$ gradient
  MatrixListType ConvolveGradTemplate(MatrixListType &gradTemplate_L2);

  /// Updates the data domain of the deformation and the kernel.
  void UpdateDeformationAndKernelDataDomain(const std::vector<std::shared_ptr<DeformableMultiObjectType> > target);

  /// Converts a matrix of size N x Dimension to a vector \e V of length Dimension x N.
  VectorType Vectorize(const MatrixType &M) const { return M.vectorise_row_wise(); }
  /// Converts a vector \e V of length Dimension x N to a matrix of size N x Dimension. The inverse operation of Vectorize.
  MatrixType Unvectorize(const VectorType &V) const {
    unsigned int nrow = V.size() / Dimension;
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

  /// Flag which freezes the template when true.
  bool m_FreezeTemplateFlag;
  /// Flag wich freezes the control points when true.
  bool m_FreezeControlPointsFlag;

  /// Number of objects in the template.
  unsigned int m_NumberOfObjects;
  /// CP spacing used if no set of control points is given.
  ScalarType m_CPSpacing;

  /// Data sigma squared for each object.
  VectorType m_DataSigmaSquared;

  /// Working deformation to compute regression deformation.
  std::shared_ptr<DiffeosType> m_Def;

  /// Kernel width to compute the Sobolev gradient of the functional/likelihood w.r.t. template variable.
  ScalarType m_SmoothingKernelWidth;

  /// Flag that indicates wheter cp and momenta trajectories should be written.
  bool m_WriteFullTrajectoriesFlag;

}; /* class Regression */

