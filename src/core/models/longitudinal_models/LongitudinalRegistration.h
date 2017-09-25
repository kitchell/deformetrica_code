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
#include "LongitudinalDataSet.h"
#include "TimeSeriesDataSet.h"
#include "Diffeos.h"

/// Support files.
#include "MultiScalarNormalDistribution.h"
#include "MultiScalarInverseWishartDistribution.h"
#include "LinearAlgebra.h"

/// Input-output files.
#include "DeformationFieldIO.h"
#include "MatrixDLM.h"
#include "MatrixDLM.h"

/// Librairies files.
#include <boost/filesystem.hpp>
#include "itkImage.h"
#include "itkSimpleFastMutexLock.h"

using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      LongitudinalMatching object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    See [Schiratti et al. 2015].
 */
template<class ScalarType, unsigned int Dimension>
class LongitudinalRegistration : public AbstractStatisticalModel<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Statistical model type.
  typedef AbstractStatisticalModel<ScalarType, Dimension> Superclass;

  /// Longitudinal data set type. 
  typedef LongitudinalDataSet<ScalarType, Dimension> LongitudinalDataSetType;
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
  LongitudinalRegistration();

  /// Copy constructor.
  LongitudinalRegistration(const LongitudinalRegistration &other);

  /// Makes a copy of the object.
  std::shared_ptr<LongitudinalRegistration> Clone() const { return std::static_pointer_cast<LongitudinalRegistration>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<Superclass> doClone() const {
    return std::static_pointer_cast<Superclass>(std::make_shared<LongitudinalRegistration>(*this));
  };
 public:

  /// Destructor
  ~LongitudinalRegistration();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

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

  /// Gets the bounding box.
  MatrixType GetBoundingBox() const { return m_BoundingBox; }
  /// Sets the bounding box.
  void SetBoundingBox(const MatrixType &bb) {
    m_BoundingBox = bb;
    m_PerpendicularDeformation->SetDataDomain(bb);
    m_ForwardReferenceGeodesic->SetDataDomain(bb);
    m_BackwardReferenceGeodesic->SetDataDomain(bb);
    KernelFactoryType *kfac = KernelFactoryType::Instantiate();
    kfac->SetDataDomain(bb);
  }

  /// Sets diffeos to deform the template.
  void SetDiffeos(std::shared_ptr<DiffeosType> def) {
    m_PerpendicularDeformation = def;
    m_ForwardReferenceGeodesic = def->Clone();
    m_BackwardReferenceGeodesic = def->Clone();
  }

  /// Sets the concentration of time points for the reference geodesic.
  void SetConcentrationOfTimePointsForReferenceGeodesic(const unsigned int c) {
    m_ConcentrationOfTimePointsForReferenceGeodesic = c;
  }
  /// Sets the number of time points for the exponentiation operation.
  void SetNumberOfTimePointsForExponentiation(const unsigned int n) { m_NumberOfTimePointsForExponentiation = n; }
  /// Sets the margin on the geodesic length.
  void SetMarginOnGeodesicLength(const ScalarType m) { m_MarginOnGeodesicLength = m; }

  /// Returns the number of objects in the template.
  unsigned int GetNumberOfObjects() const { return m_NumberOfObjects; }

  /// Gets the dimension of the discretized objects.
  std::vector<unsigned long> GetDimensionOfDiscretizedObjects() const { return m_DimensionOfDiscretizedObjects; }
  /// Sets the dimension of the discretized objects.
  void SetDimensionOfDiscretizedObjects(const std::vector<unsigned long> ddo) { m_DimensionOfDiscretizedObjects = ddo; }

  /// Gets the number of sources.
  unsigned int GetNumberOfSources() const { return m_NumberOfSources; }
  /// Sets the number of sources.
  void SetNumberOfSources(const unsigned int ns) { m_NumberOfSources = ns; }

  /// Returns the template deformable objects.
  std::shared_ptr<DeformableMultiObjectType> GetTemplate() const { return m_Template; }
  /// Sets the template deformable objects.
  void SetTemplate(std::shared_ptr<DeformableMultiObjectType> const temp) { m_Template = temp; }

  /// Returns image intensity and landmark point coordinates of the template.
  MatrixListType GetTemplateData() const { return m_Template->GetImageIntensityAndLandmarkPointCoordinates(); }
  /// Sets image intensity and landmark point coordinates of the template.
  void SetTemplateData(const MatrixListType &tempData) {
    m_Template->UpdateImageIntensityAndLandmarkPointCoordinates(tempData);
    m_Template->Update();
  }

  /// Returns control points positions.
  MatrixType GetControlPoints() const { return m_ControlPoints; }
  /// Updates positions of control points.
  void SetControlPoints(MatrixType const &cp) { m_ControlPoints = cp; }
  /// Sets control points stored in the file \e fn
  void SetControlPoints(const std::string &fn);

  /// Returns the momenta.
  MatrixType GetMomenta() const { return m_Momenta; }
  /// Sets the momenta to \e M.
  void SetMomenta(MatrixType const &mom) { m_Momenta = mom; }
  /// Sets the momenta by reading the file \e fn.
  void SetMomenta(std::string const &fn);

  /// Gets the modulation matrix.
  MatrixType GetModulationMatrix() const { return m_ModulationMatrix; }
  /// Sets the modulation matrix.
  void SetModulationMatrix(MatrixType const &modMat) { m_ModulationMatrix = modMat; }
  /// Sets the modulation matrix stored in the file \e fn
  void SetModulationMatrix(const std::string &fn);

  /// Gets the reference time.
  ScalarType GetReferenceTime() const { return m_ReferenceTime; }
  /// Sets the reference time.
  void SetReferenceTime(ScalarType const &refTime) { m_ReferenceTime = refTime; }

  /// Gets the time shift variance.
  ScalarType GetTimeShiftVariance() const { return m_TimeShiftVariance; }
  /// Sets the time shift variance.
  void SetTimeShiftVariance(ScalarType const &var) { m_TimeShiftVariance = var; }

  /// Gets the log acceleration variance.
  ScalarType GetLogAccelerationVariance() const { return m_LogAccelerationVariance; }
  /// Sets the log acceleration variance.
  void SetLogAccelerationVariance(ScalarType const &var) { m_LogAccelerationVariance = var; }

  /// Gets the noise variance.
  VectorType GetNoiseVariance() const { return m_NoiseVariance; }
  /// Sets the noise variance.
  void SetNoiseVariance(VectorType const &var) { m_NoiseVariance = var; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates the bounding box and initializes the control points if needed.
  virtual void Update();

  /// Prints information about the LongitudinalMatching model.
  virtual void Print() const;
  /// Saves the LongitudinalMatching model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const;

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                const LinearVariableMapType &popRER,
                                const LinearVariablesMapType &indRER,
                                std::vector<std::vector<std::vector<ScalarType>>> &residuals) {
    residuals.resize(1);
    return ComputeResiduals(dataSet, residuals[0]);
  }
  /// Computes the residuals. Specific overloading.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                std::vector<std::vector<ScalarType>> &residuals);

  /// Evaluates how well the LongitudinalMatching model explains a specific time series dataset.
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Saves the current template shape(s).
  void WriteModelPredictions(const LongitudinalDataSetType *const dataSet) const;
  /// Saves LongitudinalMatching parameters.
  void WriteModelParameters() const;

  /// Updates the memorized absolute times.
  void UpdateAbsoluteTimeIncrements(const LongitudinalDataSetType *const dataSet);
  /// Updates the memorized reference geodesic.
  void UpdateReferenceGeodesic();

  /// Returns the deformed object and associated control points for a given time.
  void GetDeformedObjectAndControlPointsAt(ScalarType const &time,
                                           std::shared_ptr<DeformableMultiObjectType> &shape,
                                           MatrixType &cp) const;

  /// Initializes the bounding box, based on the template data and the control points.
  void InitializeBoundingBox();
  /// Projects the modulation matrix and save the result in the model attributes.
  void ProjectModulationMatrix();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s). They need to be specified.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Template object.
  std::shared_ptr<DeformableMultiObjectType> m_Template;
  /// Template objects name.
  std::vector<std::string> m_TemplateObjectsName;
  /// Extension of template objects name.
  std::vector<std::string> m_TemplateObjectsNameExtension;
  /// Template bounding box: the union of the template bounding box and the bounding box around control points.
  MatrixType m_BoundingBox;

  /// Control points.
  MatrixType m_ControlPoints;
  /// Momenta.
  MatrixType m_Momenta;
  /// Modulation matrix.
  MatrixType m_ModulationMatrix;
  /// Projected modulation matrix.
  MatrixType m_ProjectedModulationMatrix;
  /// Reference time.
  ScalarType m_ReferenceTime;

  /// Time-shift prior variance.
  ScalarType m_TimeShiftVariance;
  /// Log-acceleration prior variance.
  ScalarType m_LogAccelerationVariance;
  /// Noise prior variance.
  VectorType m_NoiseVariance;

  /// Number of time points per unit of time to be used for the reference geodesic computation.
  unsigned int m_ConcentrationOfTimePointsForReferenceGeodesic;
  /// Number of time points to be used for the Riemannian exponential computations.
  unsigned int m_NumberOfTimePointsForExponentiation;
  /// Unidirectional margin for the reference geodesic.
  ScalarType m_MarginOnGeodesicLength;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Working attributes. They could be avoided.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Number of objects in the template.
  unsigned int m_NumberOfObjects;
  /// Dimension of each discretized object.
  std::vector<unsigned long> m_DimensionOfDiscretizedObjects;
  /// Number of control points.
  unsigned int m_NumberOfControlPoints;
  /// Dimension of the tangent spaces to the manifold of diffeomorphisms = nbOfControlPoints x Dimension.
  unsigned int m_DimensionOfTangentSpaces;
  /// Number of sources in the ICA-like decomposition.
  unsigned int m_NumberOfSources;

  /// Working deformations.
  std::shared_ptr<DiffeosType> m_ForwardReferenceGeodesic;
  std::shared_ptr<DiffeosType> m_BackwardReferenceGeodesic;
  std::shared_ptr<DiffeosType> m_PerpendicularDeformation;

  /// Memorized absolute time increments.
  std::vector<ScalarType> m_AbsoluteTimeIncrements;
  /// Memorized minimum absolute time increment.
  ScalarType m_MinimumAbsoluteTimeIncrement;
  /// Memorized maximum absolute time increment.
  ScalarType m_MaximumAbsoluteTimeIncrement;

}; /* class LongitudinalRegistration */


