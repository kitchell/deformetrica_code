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
 *	\brief      LongitudinalAtlas object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    See [Schiratti et al. 2015].
 */
template<class ScalarType, unsigned int Dimension>
class LongitudinalAtlas : public AbstractStatisticalModel<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Statistical model type.
  typedef AbstractStatisticalModel<ScalarType, Dimension> Superclass;

  /// Longitudinal data set type. 
  typedef LongitudinalDataSet<ScalarType, Dimension> LongitudinalDataSetType;

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
  LongitudinalAtlas();

  /// Copy constructor.
  LongitudinalAtlas(const LongitudinalAtlas &other);

  /// Makes a copy of the object.
  std::shared_ptr<LongitudinalAtlas> Clone() const { return std::static_pointer_cast<LongitudinalAtlas>(doClone()); }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<Superclass> doClone() const {
    return std::static_pointer_cast<Superclass>(std::make_shared<LongitudinalAtlas>(*this));
  };
 public:

  /// Destructor
  ~LongitudinalAtlas();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the fixed effects.
  virtual void SetFixedEffects(const LinearVariableMapType &map);

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
    m_Def->SetDataDomain(bb);
    m_ForwardReferenceGeodesic->SetDataDomain(bb);
    m_BackwardReferenceGeodesic->SetDataDomain(bb);
    KernelFactoryType *kfac = KernelFactoryType::Instantiate();
    kfac->SetDataDomain(bb);
  }

  /// Sets diffeos to deform the template.
  void SetDiffeos(std::shared_ptr<DiffeosType> def) {
    m_Def = def;
    m_ForwardReferenceGeodesic = def->Clone();
    m_BackwardReferenceGeodesic = def->Clone();
    m_ForwardReferenceGeodesic_Memory = def->Clone();
    m_BackwardReferenceGeodesic_Memory = def->Clone();
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

  /// Gets control point spacing.
  ScalarType GetCPSpacing() const { return m_CPSpacing; }
  /// Sets control point spacing. Used in case no control points have been set to define a regular lattice of control points.
  void SetCPSpacing(const ScalarType s) { m_CPSpacing = s; }

  /// Sets the freeze control points flag.
  void SetFreezeControlPointsFlag(bool const &flag) { m_FreezeControlPointsFlag = flag; }

  /// Returns the template deformable objects.
  std::shared_ptr<DeformableMultiObjectType> GetTemplate() const { return m_Template; }
  /// Sets the template deformable objects.
  void SetTemplate(std::shared_ptr<DeformableMultiObjectType> const temp) {
    m_Template = temp;
    MatrixListType templateData = temp->GetImageIntensityAndLandmarkPointCoordinates();
    this->m_FixedEffects["TemplateData"] = templateData;
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("TemplateData"))->SetMean(templateData);
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
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("TemplateData"))->SetMean(tempData);
  }

  /// Gets the template data random effect variance.
  ScalarType GetTemplateDataRandomEffectVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("TemplateData"))->GetVariance();
  }
  /// Gets the template data random effect standard deviation.
  ScalarType GetTemplateDataRandomEffectVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("TemplateData"))->GetVarianceSqrt();
  }
  /// Sets the template data random effect standard deviation.
  void SetTemplateDataRandomEffectVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("TemplateData"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the template data prior mean.
  MatrixListType GetTemplateDataPriorMean() const {
    return recast<MatrixListType>(unvectorize(std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("TemplateData"))->GetMean(), m_TemplateDataSizeParameters));
  }
  /// Sets the template data prior mean.
  void SetTemplateDataPriorMean(MatrixListType const &tempData) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("TemplateData"))->SetMean(tempData);
  }
  /// Gets the template data prior variance.
  ScalarType GetTemplateDataPriorVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("TemplateData"))->GetVariance();
  }
  /// Gets the template data prior standard deviation.
  ScalarType GetTemplateDataPriorVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("TemplateData"))->GetVarianceSqrt();
  }
  /// Sets the template data random effect standard deviation.
  void SetTemplateDataPriorVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("TemplateData"))->SetVarianceSqrt(varSqrt);
  }

  /// Returns control points positions.
  MatrixType GetControlPoints() const {
    if (!m_FreezeControlPointsFlag) { return recast<MatrixType>(Superclass::m_FixedEffects.at("ControlPoints")); }
    else { return m_FrozenControlPoints; }
  }
  /// Updates positions of control points.
  void SetControlPoints(MatrixType const &cp) {
    if (!m_FreezeControlPointsFlag) {
      Superclass::m_FixedEffects["ControlPoints"] = cp;
      std::static_pointer_cast<MultiScalarNormalDistributionType>(
          this->m_PopulationRandomEffects.at("ControlPoints"))->SetMean(cp);
    } else { m_FrozenControlPoints = cp; }
  }
  /// Sets control points stored in the file \e fn
  void SetControlPoints(const std::string &fn);

  /// Gets the control points random effect variance.
  ScalarType GetControlPointsRandomEffectVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->GetVariance();
  }
  /// Gets the control points random effect standard deviation.
  ScalarType GetControlPointsRandomEffectVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->GetVarianceSqrt();
  }
  /// Sets the control points random effect standard deviation.
  void SetControlPointsRandomEffectVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ControlPoints"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the control points prior mean.
  MatrixType GetControlPointsPriorMean() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->GetMean().unvectorize(m_NumberOfControlPoints, Dimension);
  }
  /// Sets the control points prior mean.
  void SetControlPointsPriorMean(MatrixType const &cp) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->SetMean(cp);
  }
  /// Gets the control points prior variance.
  ScalarType GetControlPointsPriorVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->GetVariance();
  }
  /// Gets the control points prior standard deviation.
  ScalarType GetControlPointsPriorVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->GetVarianceSqrt();
  }
  /// Sets the control points random effect standard deviation.
  void SetControlPointsPriorVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ControlPoints"))->SetVarianceSqrt(varSqrt);
  }

  /// Returns the momenta.
  MatrixType GetMomenta() const {
    return recast<MatrixType>(Superclass::m_FixedEffects.at("Momenta"));
  }
  /// Sets the momenta to \e M.
  void SetMomenta(MatrixType const &mom) {
    Superclass::m_FixedEffects["Momenta"] = mom;
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("Momenta"))->SetMean(mom);
  }
  /// Sets the momenta by reading the file \e fn.
  void SetMomenta(std::string const &fn);

  /// Gets the momenta random effect variance.
  ScalarType GetMomentaRandomEffectVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("Momenta"))->GetVariance();
  }
  /// Gets the momenta random effect standard deviation.
  ScalarType GetMomentaRandomEffectVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("Momenta"))->GetVarianceSqrt();
  }
  /// Sets the momenta random effect standard deviation.
  void SetMomentaRandomEffectVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("Momenta"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the momenta prior mean.
  MatrixType GetMomentaPriorMean() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("Momenta"))->GetMean().unvectorize(m_NumberOfControlPoints, Dimension);
  }
  /// Sets the momenta prior mean.
  void SetMomentaPriorMean(MatrixType const &mom) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("Momenta"))->SetMean(mom);
  }
  /// Gets the momenta prior variance.
  ScalarType GetMomentaPriorVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("Momenta"))->GetVariance();
  }
  /// Gets the momenta prior standard deviation.
  ScalarType GetMomentaPriorVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("Momenta"))->GetVarianceSqrt();
  }
  /// Sets the momenta random effect standard deviation.
  void SetMomentaPriorVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("Momenta"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the modulation matrix.
  MatrixType GetModulationMatrix() const {
    return recast<MatrixType>(Superclass::m_FixedEffects.at("ModulationMatrix"));
  }
  /// Sets the modulation matrix.
  void SetModulationMatrix(MatrixType const &modMat) {
    Superclass::m_FixedEffects["ModulationMatrix"] = modMat;
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ModulationMatrix"))->SetMean(modMat);
  }

  /// Gets the modulation matrix random effect variance.
  ScalarType GetModulationMatrixRandomEffectVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ModulationMatrix"))->GetVariance();
  }
  /// Gets the modulation matrix random effect standard deviation.
  ScalarType GetModulationMatrixRandomEffectVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ModulationMatrix"))->GetVarianceSqrt();
  }
  /// Sets the modulation matrix random effect standard deviation.
  void SetModulationMatrixRandomEffectVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_PopulationRandomEffects.at("ModulationMatrix"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the modulation matrix prior mean.
  MatrixType GetModulationMatrixPriorMean() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(this->m_Priors.at("ModulationMatrix"))
        ->GetMean().unvectorize(m_DimensionOfTangentSpaces, m_NumberOfSources);
  }
  /// Sets the modulation matrix prior mean.
  void SetModulationMatrixPriorMean(MatrixType const &modMat) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ModulationMatrix"))->SetMean(modMat);
  }
  /// Gets the modulation matrix prior variance.
  ScalarType GetModulationMatrixPriorVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ModulationMatrix"))->GetVariance();
  }
  /// Gets the modulation matrix prior standard deviation.
  ScalarType GetModulationMatrixPriorVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ModulationMatrix"))->GetVarianceSqrt();
  }
  /// Sets the modulation matrix random effect standard deviation.
  void SetModulationMatrixPriorVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ModulationMatrix"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the reference time.
  ScalarType GetReferenceTime() const {
    return recast<ScalarType>(Superclass::m_FixedEffects.at("ReferenceTime"));
  }
  /// Sets the reference time.
  void SetReferenceTime(ScalarType const &refTime) {
    Superclass::m_FixedEffects["ReferenceTime"] = refTime;
    SetTimeShiftRandomEffectMean(refTime);
  }
  /// Gets the reference time prior mean.
  ScalarType GetReferenceTimePriorMean() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ReferenceTime"))->GetMean()[0];
  }
  /// Sets the reference time prior mean.
  void SetReferenceTimePriorMean(ScalarType const &refTime) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ReferenceTime"))->SetMean(refTime);
  }
  /// Gets the reference time prior variance.
  ScalarType GetReferenceTimePriorVariance() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ReferenceTime"))->GetVariance();
  }
  /// Gets the reference time prior standard deviation.
  ScalarType GetReferenceTimePriorVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ReferenceTime"))->GetVarianceSqrt();
  }
  /// Sets the reference time random effect standard deviation.
  void SetReferenceTimePriorVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_Priors.at("ReferenceTime"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the time shift variance.
  ScalarType GetTimeShiftVariance() const {
    return recast<ScalarType>(Superclass::m_FixedEffects.at("TimeShiftVariance"));
  }
  /// Sets the time shift variance.
  void SetTimeShiftVariance(ScalarType const &var) {
    Superclass::m_FixedEffects["TimeShiftVariance"] = var;
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("TimeShift"))->SetVariance(var);
  }
  /// Gets the time shift random effect mean.
  ScalarType GetTimeShiftRandomEffectMean() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("TimeShift"))->GetMean()(0);
  }
  /// Sets the time shift random effect mean.
  void SetTimeShiftRandomEffectMean(ScalarType const &tsm) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("TimeShift"))->SetMean(tsm);
  }
  /// Gets the time shift random effect standard deviation.
  ScalarType GetTimeShiftRandomEffectVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("TimeShift"))->GetVarianceSqrt();
  }
  /// Sets the time shift random effect standard deviation.
  void SetTimeShiftRandomEffectVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("TimeShift"))->SetVarianceSqrt(varSqrt);
  }
  /// Gets the time shift variance prior scale factor.
  ScalarType GetTimeShiftVariancePriorScaleFactor() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("TimeShiftVariance"))->GetScaleVector()[0];
  }
  /// Sets the time shift variance prior scale factor.
  void SetTimeShiftVariancePriorScaleFactor(ScalarType const &scale) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("TimeShiftVariance"))->SetScaleVector(VectorType(1, scale));
  }
  /// Gets the time shift variance prior standard deviation.
  ScalarType GetTimeShiftVariancePriorDegreeOfFreedom() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("TimeShiftVariance"))->GetDegreesOfFreedom()[0];
  }
  /// Sets the time shift variance random effect standard deviation.
  void SetTimeShiftVariancePriorDegreeOfFreedom(ScalarType const &dof) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("TimeShiftVariance"))->SetDegreesOfFreedom(VectorType(1, dof));
  }

  /// Gets the log acceleration variance.
  ScalarType GetLogAccelerationVariance() const {
    return recast<ScalarType>(Superclass::m_FixedEffects.at("LogAccelerationVariance"));
  }
  /// Sets the log acceleration variance.
  void SetLogAccelerationVariance(ScalarType const &var) {
    Superclass::m_FixedEffects["LogAccelerationVariance"] = var;
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("LogAcceleration"))->SetVariance(var);
  }

  /// Sets the log acceleration random effect mean.
  void SetLogAccelerationRandomEffectMean(VectorType const &logAcc) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("LogAcceleration"))->SetMean(logAcc);
  }
  /// Gets the log acceleration random effect standard deviation.
  ScalarType GetLogAccelerationRandomEffectVarianceSqrt() const {
    return std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("LogAcceleration"))->GetVarianceSqrt();
  }
  /// Sets the log acceleration random effect standard deviation.
  void SetLogAccelerationRandomEffectVarianceSqrt(ScalarType const &varSqrt) {
    std::static_pointer_cast<MultiScalarNormalDistributionType>(
        this->m_IndividualRandomEffects.at("LogAcceleration"))->SetVarianceSqrt(varSqrt);
  }

  /// Gets the log acceleration variance prior scale factor.
  ScalarType GetLogAccelerationVariancePriorScaleFactor() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("LogAccelerationVariance"))->GetScaleVector()[0];
  }
  /// Sets the log acceleration variance prior scale factor.
  void SetLogAccelerationVariancePriorScaleFactor(ScalarType const &scale) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("LogAccelerationVariance"))->SetScaleVector(VectorType(1, scale));
  }
  /// Gets the log acceleration variance prior standard deviation.
  ScalarType GetLogAccelerationVariancePriorDegreeOfFreedom() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("LogAccelerationVariance"))->GetDegreesOfFreedom()[0];
  }
  /// Sets the log acceleration variance random effect standard deviation.
  void SetLogAccelerationVariancePriorDegreeOfFreedom(ScalarType const &dof) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("LogAccelerationVariance"))->SetDegreesOfFreedom(VectorType(1, dof));
  }

  /// Gets the noise variance.
  VectorType GetNoiseVariance() const {
    return recast<VectorType>(Superclass::m_FixedEffects.at("NoiseVariance"));
  }
  /// Sets the noise variance.
  void SetNoiseVariance(VectorType const &var) {
    Superclass::m_FixedEffects["NoiseVariance"] = var;
  }

  /// Gets the noise variance prior scale vector.
  VectorType GetNoiseVariancePriorScaleVector() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->GetScaleVector();
  }
  /// Sets the noise variance prior scale factor.
  void SetNoiseVariancePriorScaleVector(VectorType const &scaleVec) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->SetScaleVector(scaleVec);
  }
  /// Gets the noise variance prior standard deviation.
  VectorType GetNoiseVariancePriorDegreesOfFreedom() const {
    return std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->GetDegreesOfFreedom();
  }
  /// Sets the noise variance random effect standard deviation.
  void SetNoiseVariancePriorDegreesOfFreedom(VectorType const &dofs) {
    std::static_pointer_cast<MultiScalarInverseWishartDistributionType>(
        this->m_Priors.at("NoiseVariance"))->SetDegreesOfFreedom(dofs);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates the bounding box and initializes the control points if needed.
  virtual void Update();

  /// Prints information about the LongitudinalAtlas model.
  virtual void Print() const;
  /// Saves the LongitudinalAtlas model.
  virtual void Write(const LongitudinalDataSetType *const dataSet,
                     LinearVariableMapType const &popRER,
                     LinearVariablesMapType const &indRER) const;

  /// Samples from the model, and stores the used parameters in \e popRER and \e indRER.
  virtual bool Sample(LongitudinalDataSetType *dataSet,
                      LinearVariableMapType &popRER,
                      LinearVariablesMapType &indRER) const;

  /// Computes the residuals.
  virtual bool ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                                const LinearVariableMapType &popRER,
                                const LinearVariablesMapType &indRER,
                                std::vector<std::vector<std::vector<ScalarType>>> &residuals);

  /// Evaluates how well the LongitudinalAtlas model explains a specific time series dataset.
  virtual ScalarType ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                  const LinearVariableMapType &popRER,
                                                  const LinearVariablesMapType &indRER);
  /// Computes the functional first mode, given an input random effects realization.
  virtual bool UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                                 const LinearVariableMapType &popRER,
                                                                 const LinearVariablesMapType &indRER,
                                                                 VectorType &logLikelihoodTerms) {
    logLikelihoodTerms.set_size(2);
    logLikelihoodTerms[0] = ComputeCompleteLogLikelihood(dataSet, popRER, indRER);
    logLikelihoodTerms[1] = 0.0;
    return false;
  }

  /// Computes the model log-likelihood, returning as well the detailed contribution of each subject.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               std::vector<ScalarType> &contributions) {
    return ComputeModelLogLikelihood(dataSet, popRER, indRER, 1.0, "All", contributions);
  };
  /// Optimized version, where only one population random effect realization (RER) has been modified since the last call.
  virtual ScalarType ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                               const LinearVariableMapType &popRER,
                                               const LinearVariablesMapType &indRER,
                                               const ScalarType &temperature,
                                               const std::string &modifiedVar,
                                               std::vector<ScalarType> &contributions);

  /// Computes the log-likelihood mode, ignoring contributions of priors and subjects other than \e i.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i) {
    return ComputeModelLogLikelihoodForSubject(dataSet, popRER, indRER, i, "All");
  }
  /// Optimized version, where only one individual random effect realization (RER) has been modified since the last call.
  virtual ScalarType ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                                         const LinearVariableMapType &popRER,
                                                         const LinearVariablesMapType &indRER,
                                                         const unsigned int &i,
                                                         const std::string &modifiedVar);

  /// Computes the gradient of the log-likelihood first mode.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad) {
    std::cerr << "Error : LongitudinalAtlas::ComputeCompleteLogLikelihoodGradient method is not available."
              << std::endl;
  }
  /// Computes the gradient of the log-likelihood first mode and the corresponding norms.
  virtual void ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    LinearVariableMapType &popGrad,
                                                    LinearVariablesMapType &indGrad,
                                                    VectorType &gradSquaredNorms) {
    std::cerr << "Error : LongitudinalAtlas::ComputeCompleteLogLikelihoodGradient method is not available."
              << std::endl;
  }
  /// Computes the marginal log-likelihood gradient with respect to the random effects of the model.
  virtual void ComputeLogLikelihoodGradientWrtRandomEffects(const LongitudinalDataSetType *const dataSet,
                                                            const LinearVariableMapType &popRER,
                                                            const LinearVariablesMapType &indRER,
                                                            LinearVariableMapType &popGrad,
                                                            LinearVariablesMapType &indGrad) {
    std::cerr << "Error : LongitudinalAtlas::ComputeLogLikelihoodGradientWrtRandomEffects is not available."
              << std::endl;
  }

  /// Computes the sufficient statistics of the LongitudinalAtlas model.
  virtual void ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                                           const LinearVariableMapType &popRER,
                                           const LinearVariablesMapType &indRER,
                                           LinearVariableMapType &sufficientStatistics);

  /// Updates the fixed effects based on the given sufficient statistics, maximizing the likelihood.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics) {
    return UpdateFixedEffects(dataSet, sufficientStatistics, 1);
  }
  /// Optional variation, featuring a temperature parameter.
  virtual void UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                                  const LinearVariableMapType &sufficientStatistics,
                                  const ScalarType &temperature);

  /// Recovers the memorized random effect realizations-based state.
  virtual void RecoverMemorizedState();

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Saves the current template shape(s).
  void WriteModelPredictions(const LongitudinalDataSetType *const dataSet,
                             LinearVariableMapType const &popRER,
                             LinearVariablesMapType const &indRER) const;
  /// Saves LongitudinalAtlas parameters.
  void WriteModelParameters(LinearVariablesMapType const &indRER) const;

  /// Computes the residuals for a specific subject.
  std::vector<std::vector<ScalarType>> ComputeResidualsForSubject(const LongitudinalDataSetType *const dataSet,
                                                                  const LinearVariableMapType &popRER,
                                                                  const LinearVariablesMapType &indRER,
                                                                  const unsigned int &i);

  /// Updates the memorized absolute times.
  void UpdateAbsoluteTimeIncrements(const LongitudinalDataSetType *const dataSet,
                                    const LinearVariableMapType &popRER,
                                    const LinearVariablesMapType &indRER);
  /// Updates the memorized reference geodesic.
  void UpdateReferenceGeodesic(const LinearVariableMapType &popRER);
  /// Updates the memorized projected modulation matrix.
  void UpdateProjectedModulationMatrix(const LinearVariableMapType &popRER);

  /// Returns the deformed object and associated control points for a given time.
  void GetDeformedObjectAndControlPointsAt(ScalarType const &time,
                                           std::shared_ptr<DeformableMultiObjectType> &shape,
                                           MatrixType &cp) const;

  /// Initializes the bounding box, based on the template data and the control points.
  void InitializeBoundingBox();

  /// Initializes the template data-related variables.
  void InitializeTemplateDataVariables();
  /// Initializes the control points-related variables.
  void InitializeControlPointsVariables();
  /// Initializes the momenta-related variables.
  void InitializeMomentaVariables();
  /// Initializes the modulation matrix-related variables.
  void InitializeModulationMatrixVariables();
  /// Initializes the reference time-related variables.
  void InitializeReferenceTimeVariables();
  /// Initializes the sources.
  void InitializeSources();
  /// Initializes the time shift-related variables.
  void InitializeTimeShiftVariables();
  /// Initializes the log acceleration-related variables.
  void InitializeLogAccelerationVariables();
  /// Initializes the noise-related variables.
  void InitializeNoiseVariables();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Template object.
  std::shared_ptr<DeformableMultiObjectType> m_Template;
  /// Template objects name.
  std::vector<std::string> m_TemplateObjectsName;
  /// Extension of template objects name.
  std::vector<std::string> m_TemplateObjectsNameExtension;
  /// Size parameters of the template data matrix list.
  std::vector<unsigned int> m_TemplateDataSizeParameters;
  /// Template bounding box: the union of the template bounding box and the bounding box around control points.
  MatrixType m_BoundingBox;

  /// Working deformation to compute LongitudinalAtlas deformation.
  std::shared_ptr<DiffeosType> m_Def;
  /// Number of time points per unit of time to be used for the reference geodesic computation.
  unsigned int m_ConcentrationOfTimePointsForReferenceGeodesic;
  /// Number of time points to be used for the Riemannian exponential computations.
  unsigned int m_NumberOfTimePointsForExponentiation;
  /// Unidirectional margin for the reference geodesic.
  ScalarType m_MarginOnGeodesicLength;

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

  /// CP spacing used if no set of control points is given.
  ScalarType m_CPSpacing;
  /// Flag which freezes the control points when true.
  bool m_FreezeControlPointsFlag;
  /// If the previous flag is active, stores the frozen control points.
  MatrixType m_FrozenControlPoints;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Memory attributes, for optimization of the log-likelihood computations.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Memorized absolute times.
  std::vector<std::vector<ScalarType>> m_AbsoluteTimeIncrements;
  std::vector<std::vector<ScalarType>> m_AbsoluteTimeIncrements_Memory;
  /// Memorized minimum absolute time.
  ScalarType m_MinimumAbsoluteTimeIncrement;
  ScalarType m_MinimumAbsoluteTimeIncrement_Memory;
  /// Memorized maximum absolute time.
  ScalarType m_MaximumAbsoluteTimeIncrement;
  ScalarType m_MaximumAbsoluteTimeIncrement_Memory;

  /// Forward part of the memorized reference geodesic.
  std::shared_ptr<DiffeosType> m_ForwardReferenceGeodesic;
  std::shared_ptr<DiffeosType> m_ForwardReferenceGeodesic_Memory;
  /// Backward part of the memorized reference geodesic.
  std::shared_ptr<DiffeosType> m_BackwardReferenceGeodesic;
  std::shared_ptr<DiffeosType> m_BackwardReferenceGeodesic_Memory;

  /// Memorized projected modulation matrix.
  MatrixType m_ProjectedModulationMatrix;
  MatrixType m_ProjectedModulationMatrix_Memory;

}; /* class LongitudinalAtlas */


