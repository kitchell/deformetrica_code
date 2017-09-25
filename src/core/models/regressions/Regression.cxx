/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "Regression.h"

using namespace def::algebra;


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
Regression<ScalarType, Dimension>
::Regression() : Superclass(), m_Def(NULL), m_SmoothingKernelWidth(0.0),
                 m_FreezeTemplateFlag(false), m_FreezeControlPointsFlag(false),
                 m_WriteFullTrajectoriesFlag(false) {
  Superclass::SetRegressionType();

  /// Initialization of the fixed effects.
  MatrixType controlPoints, initialMomenta;
  MatrixListType tempData;
  Superclass::m_FixedEffects["ControlPoints"] = controlPoints;
  Superclass::m_FixedEffects["InitialMomenta"] = initialMomenta;
  Superclass::m_FixedEffects["TemplateData"] = tempData;
}

template<class ScalarType, unsigned int Dimension>
Regression<ScalarType, Dimension>
::~Regression() {}

template<class ScalarType, unsigned int Dimension>
Regression<ScalarType, Dimension>
::Regression(const Regression &other) {
  m_Def = std::static_pointer_cast<DiffeosType>(other.m_Def->Clone());
  m_SmoothingKernelWidth = other.m_SmoothingKernelWidth;
  m_FreezeTemplateFlag = other.m_FreezeTemplateFlag;
  m_FreezeControlPointsFlag = other.m_FreezeControlPointsFlag;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::SetFixedEffects(const LinearVariableMapType &map) {
  if (m_FreezeTemplateFlag) {
    MatrixListType const aux = GetTemplateData();
    Superclass::SetFixedEffects(map);
    Superclass::m_FixedEffects["TemplateData"] = aux;
  } else {
    Superclass::SetFixedEffects(map);
    m_Template->UpdateImageIntensityAndLandmarkPointCoordinates(GetTemplateData());
    m_Template->Update();
  }
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::SetControlPoints(const std::string &fn) {
  if (strlen(fn.c_str())) {
    const MatrixType controlPoints = readMatrixDLM<ScalarType>(fn.c_str());
    Superclass::m_FixedEffects["ControlPoints"] = controlPoints;
    std::cout << "Using a set of " << controlPoints.rows() << " control points in file " << fn << std::endl;
  }
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::SetInitialMomenta(const std::string &fn) {
  if (strlen(fn.c_str())) {
    const MatrixType initialMomenta = readMatrixDLM<ScalarType>(fn.c_str());
    Superclass::m_FixedEffects["InitialMomenta"] = initialMomenta;
    std::cout << "Using an initial momenta matrix of size " << initialMomenta.rows() << " x "
              << initialMomenta.columns() << " from file: " << fn << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::Write(const LongitudinalDataSetType *const dataSet,
        LinearVariableMapType const &popRER,
        LinearVariablesMapType const &indRER) const {
  WriteTemplateFlow();
  WriteRegressionParameters();
  WriteTemplateData();
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::Update() {
  if (m_Def == NULL)
    throw std::runtime_error("A deformation should be set to the regression model");

  if (GetControlPoints().rows() == 0) {
    InitializeControlPoints();
  }

  MatrixType controlPoints = GetControlPoints();
  MatrixType initialMomenta = GetInitialMomenta();
  std::shared_ptr<DeformableMultiObjectType> temp = GetTemplate();

  temp->Update();
  m_NumberOfObjects = temp->GetNumberOfObjects();
  m_BoundingBox = temp->GetBoundingBox();

  for (unsigned int i = 0; i < controlPoints.rows(); i++) {
    for (unsigned int dim = 0; dim < Dimension; dim++) {
      m_BoundingBox(dim, 0) =
          (m_BoundingBox(dim, 0) < controlPoints(i, dim) ? m_BoundingBox(dim, 0) : controlPoints(i, dim));
      m_BoundingBox(dim, 1) =
          (m_BoundingBox(dim, 1) > controlPoints(i, dim) ? m_BoundingBox(dim, 1) : controlPoints(i, dim));
    }
  }

  if (!((initialMomenta.rows() == controlPoints.rows()) && (initialMomenta.columns() == Dimension))) {
    if (initialMomenta.rows() > 0)
      std::cout << "Warning: initial initialMomenta file has incompatible number of vectors. "
          "Initial momenta reset to zero." << std::endl;

    initialMomenta.set_size(controlPoints.rows(), Dimension);
    initialMomenta.fill(0.0);

    SetInitialMomenta(initialMomenta);
  }
}

template<class ScalarType, unsigned int Dimension>
bool
Regression<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<std::vector<ScalarType>>> &residuals) {
  std::vector<std::vector<ScalarType>> aux;
  bool oob = ComputeResiduals(dataSet, popRER, indRER, aux);
  residuals.resize(1);
  residuals[0] = aux;
  return oob;
}

template<class ScalarType, unsigned int Dimension>
bool
Regression<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   const LinearVariableMapType &popRER,
                   const LinearVariablesMapType &indRER,
                   std::vector<std::vector<ScalarType>> &residuals) {
  std::vector<std::shared_ptr<DeformableMultiObjectType>> targets;
  std::vector<unsigned int> timeIndices;
  const TimeSeriesDataSetType *const timeSeriesDataSet = static_cast<const TimeSeriesDataSetType *const>(dataSet);
  timeSeriesDataSet->GetDataTimeSeries(targets, timeIndices);

  if (targets.size() != timeIndices.size())
    throw std::runtime_error("Mismatch between number of observations and number of time indices.");

  // A copy of the deformation is needed since this method can be called by different threads (maybe?) at the same time
  UpdateDeformationAndKernelDataDomain(targets);
  std::shared_ptr<DiffeosType> subjectDef = m_Def->Clone();
  subjectDef->SetDeformableMultiObject(GetTemplate());
  subjectDef->SetStartPositions(GetControlPoints());
  subjectDef->SetStartMomentas(GetInitialMomenta());
  subjectDef->Update();

  if (subjectDef->OutOfBox())
    return true;

  unsigned int numberOfObservations = targets.size();
  residuals.resize(numberOfObservations);

  // Loop over the number of observations to compute the value of the data matching term
  for (unsigned int t = 0; t < numberOfObservations; ++t)
    residuals[t] = subjectDef->GetDeformedObjectAt(timeIndices[t])->ComputeMatch(targets[t]);

  return false;
}

template<class ScalarType, unsigned int Dimension>
bool
Regression<ScalarType, Dimension>
::UpdateFixedEffectsAndComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                                                    const LinearVariableMapType &popRER,
                                                    const LinearVariablesMapType &indRER,
                                                    VectorType &logLikelihoodTerms) {
  /// Initialization.
  logLikelihoodTerms.set_size(2);
  logLikelihoodTerms.fill(0.0);

  /// Data (residuals) term.
  std::vector<std::vector<ScalarType>> residuals;
  bool oob = ComputeResiduals(dataSet, popRER, indRER, residuals);

  for (unsigned int t = 0; t < residuals.size(); ++t)
    for (unsigned int i = 0; i < m_NumberOfObjects; i++)
      logLikelihoodTerms[0] -= 0.5 * residuals[t][i] / m_DataSigmaSquared[i];

  /// Regularity term.
  const MatrixType initialMomenta = GetInitialMomenta();
  const MatrixType controlPoints = GetControlPoints();
  const unsigned int nbControlPoints = controlPoints.rows();

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> momKernelObj = kfac->CreateKernelObject(m_Def->GetKernelType());
  momKernelObj->SetKernelWidth(m_Def->GetKernelWidth());
  momKernelObj->SetSources(controlPoints);
  momKernelObj->SetWeights(initialMomenta);

  MatrixType kMom = momKernelObj->Convolve(controlPoints);
  for (unsigned int i = 0; i < nbControlPoints; ++i)
    logLikelihoodTerms[1] -= dot_product(kMom.get_row(i), initialMomenta.get_row(i));
  logLikelihoodTerms[1] *= 0.5f;

  return oob; // Out of box flag.
}

template<class ScalarType, unsigned int Dimension>
ScalarType
Regression<ScalarType, Dimension>
::ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                            const LinearVariableMapType &popRER,
                            const LinearVariablesMapType &indRER,
                            std::vector<ScalarType> &contributions) {
  /// Data (residuals) term.
  std::vector<std::vector<ScalarType>> residuals;
  ComputeResiduals(dataSet, popRER, indRER, residuals);

  ScalarType out = 0.0;
  for (unsigned int t = 0; t < residuals.size(); ++t)
    for (unsigned int i = 0; i < m_NumberOfObjects; i++)
      out -= 0.5 * residuals[t][i] / m_DataSigmaSquared[i];

  contributions.resize(1);
  contributions[0] = out;

  return out;
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad) {
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target;
  std::vector<unsigned int> timeIndices;
  const TimeSeriesDataSetType *const timeSeriesDataSet = static_cast<const TimeSeriesDataSetType *const>(dataSet);
  timeSeriesDataSet->GetDataTimeSeries(target, timeIndices);
  UpdateDeformationAndKernelDataDomain(target);

  const MatrixType initialMomenta = GetInitialMomenta();
  MatrixType gradPos, gradMom;
  MatrixListType gradTempL_L2, gradTempL_Sob;

  ComputeDataTermGradient(initialMomenta, target, timeIndices, gradPos, gradMom, gradTempL_L2, gradTempL_Sob);
  AddGradientRegularityTerm(gradPos, gradMom, initialMomenta);

  if (!m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  popGrad["InitialMomenta"] = gradMom;
  if (!m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::ComputeCompleteLogLikelihoodGradient(const LongitudinalDataSetType *const dataSet,
                                       const LinearVariableMapType &popRER,
                                       const LinearVariablesMapType &indRER,
                                       LinearVariableMapType &popGrad,
                                       LinearVariablesMapType &indGrad,
                                       VectorType &gradSquaredNorms) {
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target;
  std::vector<unsigned int> timeIndices;
  const TimeSeriesDataSetType *const timeSeriesDataSet = static_cast<const TimeSeriesDataSetType *const>(dataSet);
  timeSeriesDataSet->GetDataTimeSeries(target, timeIndices);
  UpdateDeformationAndKernelDataDomain(target);

  const MatrixType initialMomenta = GetInitialMomenta();
  MatrixType gradPos, gradMom;
  MatrixListType gradTempL_L2, gradTempL_Sob;

  ComputeDataTermGradient(initialMomenta, target, timeIndices, gradPos, gradMom, gradTempL_L2, gradTempL_Sob);
  AddGradientRegularityTerm(gradPos, gradMom, initialMomenta);

  if (!m_FreezeControlPointsFlag) { popGrad["ControlPoints"] = gradPos; }
  popGrad["InitialMomenta"] = gradMom;
  if (!m_FreezeTemplateFlag) { popGrad["TemplateData"] = gradTempL_Sob; }

  // Must be in alphabetical order !!
  gradSquaredNorms.set_size(0);
  if (!m_FreezeControlPointsFlag) { gradSquaredNorms.push_back(gradPos.sum_of_squares()); }
  gradSquaredNorms.push_back(gradMom.sum_of_squares());
  if (!m_FreezeTemplateFlag) { gradSquaredNorms.push_back(dot_product(gradTempL_L2, gradTempL_Sob)); }
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                              const LinearVariableMapType &popRER,
                              const LinearVariablesMapType &indRER,
                              LinearVariableMapType &sufficientStatistics) {
//    this->UpdateDeformationAndKernelDataDomain(dataSet->GetDeformableMultiObjects());
//
//    /// Random control points.
//    sufficientStatistics["S1"] = popRER["ControlPoints"];
//
//    /// Applies the deformation to the template.
//	const std::vector<std::shared_ptr<DeformableMultiObjectType>> targets = dataSet->GetDeformableMultiObjects();
//	const std::vector<unsigned int> timeIndices = dataSet->GetTimeIndices();
//	const unsigned int imgIndex = m_Template->GetImageIndex();
//    const MatrixType initialMomenta = GetInitialMomenta();
//
//	std::shared_ptr<DiffeosType> subjectDef = m_Def->Clone();
//	subjectDef->SetDeformableMultiObject(GetTemplate());
//	subjectDef->SetStartPositions(GetControlPoints());
//	subjectDef->SetStartMomentas(initialMomenta);
//	subjectDef->Update();
//
//	/// Initialization with the first observation.
//	std::shared_ptr<const ParametricImageType> deformedTemplate
//		= std::static_pointer_cast<const ParametricImageType>(
//			subjectDef->GetDeformedObjectAt(timeIndices[0])->GetObjectList()[imgIndex]);
//	MatrixType deformedKernelField = deformedTemplate->GetPhotometricKernelField();
//	MatrixType Ii = targets[0]->GetImageIntensityAndLandmarkPointCoordinates()[imgIndex];
//
//	sufficientStatistics["S3"] = deformedKernelField.transpose() * Ii;
//	sufficientStatistics["S4"] = deformedKernelField.transpose() * deformedKernelField;
//
//	/// Loops and sums over the following observations.
//	for (unsigned int t = 1 ; t < timeIndices.size() ; ++t)
//	{
//		deformedTemplate = std::static_pointer_cast<const ParametricImageType>(
//				subjectDef->GetDeformedObjectAt(timeIndices[t])->GetObjectList()[imgIndex]);
//		deformedKernelField = deformedTemplate->GetPhotometricKernelField();
//		Ii = targets[t]->GetImageIntensityAndLandmarkPointCoordinates()[imgIndex];
//
//		sufficientStatistics["S3"] += deformedKernelField.transpose() * Ii;
//		sufficientStatistics["S4"] += deformedKernelField.transpose() * deformedKernelField;
//	}
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                     const LinearVariableMapType &sufficientStatistics) {
//    const MatrixType S3 = recast<MatrixType>(sufficientStatistics["S3"]);
//    const MatrixType S4 = recast<MatrixType>(sufficientStatistics["S4"]);
//
//    /// Control points random effect mean close-form update.
//    SetControlPoints(sufficientStatistics["S1"]);
//
//	/// Optional update of the parametric template photometric weights.
//	if (!m_FreezeTemplateFlag)
//	{
//		const ScalarType ss = m_DataSigmaSquared[m_Template->GetImageIndex()];
//		const MatrixType paramTempMean = MatrixType(GetParametricTemplatePriorMean());
//		const MatrixType paramTempCovInv = GetParametricTemplatePriorCovarianceInverse();
//
//		MatrixType w = inverse_sympd(S4 + ss * paramTempCovInv) * (S3 + ss * paramTempCovInv * paramTempMean);
//		SetTemplateData(MatrixListType(std::vector<MatrixType>(1, w)));
//	}
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::InitializeControlPoints(bool optimize) {
  const std::shared_ptr<const DeformableMultiObjectType> temp = GetTemplate();
  VectorType Xmin = temp->GetBoundingBox().get_column(0);
  VectorType Xmax = temp->GetBoundingBox().get_column(1);

  /// For the optimize==true case.
  ImageTypePointer img;
  ImageSpacingType spacing;
  if (optimize && temp->IsOfImageKind()[0]) {
    img = temp->GetImage();
    spacing = img->GetSpacing();
    std::cout << "We optimize the position of the control points." << std::endl;
  } else { optimize = false; }

  std::vector<VectorType> pointList;
  VectorType v(Dimension);
  switch (Dimension) {
    case 2: {
      ScalarType offsetX = 0.5 * (Xmax[0] - Xmin[0] - m_CPSpacing * floor((Xmax[0] - Xmin[0]) / m_CPSpacing));
      ScalarType offsetY = 0.5 * (Xmax[1] - Xmin[1] - m_CPSpacing * floor((Xmax[1] - Xmin[1]) / m_CPSpacing));

      for (ScalarType x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_CPSpacing) {
        for (ScalarType y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_CPSpacing) {
          bool getOut = false;
          v[0] = x;
          v[1] = y;
          if (optimize) {
            ImageIndexType ind;
            ImagePointType p;
            ScalarType value;
            ScalarType toSpan = 2 * m_Def->GetKernelWidth();
            p[0] = x;
            p[1] = y;
            img->TransformPhysicalPointToIndex(p, ind);
            value = img->GetPixel(ind);
            // Here we look : if all the intensities are equal in the neighborhood of the cp, we discard it.
            for (ScalarType xTest = x - toSpan; xTest < x + toSpan; xTest += spacing[0]) {
              if (getOut) {
                break;
              }
              p[0] = xTest;
              for (ScalarType yTest = y - toSpan; yTest < y + toSpan; yTest += spacing[1]) {
                p[1] = yTest;
                if (img->TransformPhysicalPointToIndex(p, ind)
                    && std::fabs(img->GetPixel(ind) - value) > 0.01) {
                  pointList.push_back(v);
                  getOut = true;
                  break;
                }
              }
            }
          } else // Else we don't have a criterion to discard the cp : we keep it.
            pointList.push_back(v);
        }
      }
      break;
    }
    case 3: {
      ScalarType offsetX = 0.5 * (Xmax[0] - Xmin[0] - m_CPSpacing * floor((Xmax[0] - Xmin[0]) / m_CPSpacing));
      ScalarType offsetY = 0.5 * (Xmax[1] - Xmin[1] - m_CPSpacing * floor((Xmax[1] - Xmin[1]) / m_CPSpacing));
      ScalarType offsetZ = 0.5 * (Xmax[2] - Xmin[2] - m_CPSpacing * floor((Xmax[2] - Xmin[2]) / m_CPSpacing));

      for (ScalarType x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_CPSpacing)
        for (ScalarType y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_CPSpacing)
          for (ScalarType z = Xmin[2] + offsetZ; z <= Xmax[2]; z += m_CPSpacing) {
            bool getOut = false;
            v[0] = x;
            v[1] = y;
            v[2] = z;
            if (optimize) {
              ImageIndexType ind;
              ImagePointType p;
              ScalarType value;
              ScalarType toSpan = 2 * m_Def->GetKernelWidth();

              p[0] = x;
              p[1] = y;
              p[2] = z;
              img->TransformPhysicalPointToIndex(p, ind);
              value = img->GetPixel(ind);
              // Here we look : if all the intensities are equal in the neighborhood of the cp, we discard it.
              for (ScalarType xTest = x - toSpan; xTest < x + toSpan; xTest += spacing[0]) {
                if (getOut) { break; }
                p[0] = xTest;
                for (ScalarType yTest = y - toSpan; yTest < y + toSpan; yTest += spacing[1]) {
                  if (getOut) { break; }
                  p[1] = yTest;
                  for (ScalarType zTest = z - toSpan; zTest < z + toSpan; zTest += spacing[2]) {
                    p[2] = zTest;
                    if (img->TransformPhysicalPointToIndex(p, ind)
                        && std::fabs(img->GetPixel(ind) - value) > 0.01) {
                      pointList.push_back(v);
                      getOut = true;
                      break;
                    }
                  }
                }
              }
            } else
              pointList.push_back(v);
          }
      break;
    }
    default:
      throw std::runtime_error(
          "Regression::InitializeControlPoints is not implemented in Dimensions other than 2 and 3.");
      break;
  }

  unsigned int NumCPs = pointList.size();
  MatrixType controlPoints;

  controlPoints.set_size(NumCPs, Dimension);
  for (unsigned int i = 0; i < NumCPs; i++)
    controlPoints.set_row(i, pointList[i]);

  Superclass::m_FixedEffects["ControlPoints"] = controlPoints;

  std::cout << "Set of " << NumCPs << " control points defined." << std::endl;
  if (optimize) {
    unsigned int unoptimizedNCP = 1;
    for (unsigned int d = 0; d < Dimension; ++d) { unoptimizedNCP *= 1 + floor((Xmax[d] - Xmin[d]) / m_CPSpacing); }
    std::cout << "Without optimization, " << unoptimizedNCP << " control points would have been defined." << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::WriteTemplateFlow() const {
  const std::string outputDir = def::utils::settings.output_dir;

  m_Def->SetDeformableMultiObject(GetTemplate());
  m_Def->SetStartPositions(GetControlPoints());
  m_Def->SetStartMomentas(GetInitialMomenta());
  m_Def->Update();

  if (m_Def->OutOfBox())
    throw std::runtime_error("Deformation out of box in Regression::WriteTemplateFlow: this should not be");

  std::vector<std::string> names = m_TemplateObjectsName;
  for (unsigned int i = 0; i < names.size(); i++) {
    std::ostringstream oss;
    oss << outputDir << Superclass::m_Name << "_baseline_" << names[i] << "_trajectory_";
    names[i] = oss.str();
  }

  m_Def->WriteFlow(names, m_TemplateObjectsNameExtension);

  // Output deformation field
  DeformationFieldIO<ScalarType, Dimension> defFieldIO;
  defFieldIO.SetDiffeos(m_Def);
  defFieldIO.SetAnatomicalCoordinateSystemLabel("LPS");
  defFieldIO.WriteDeformationField(outputDir + Superclass::m_Name, true);

  if (m_WriteFullTrajectoriesFlag) {
    // Write trajectories of control points and momenta.
    MatrixListType TrajectoryPositions = m_Def->GetTrajectoryPositions();
    MatrixListType TrajectoryMomentas = m_Def->GetTrajectoryMomentas();
    for (unsigned int t = 0; t < m_Def->GetNumberOfTimePoints(); ++t) {
      std::ostringstream oss;
      oss << outputDir << Superclass::m_Name << "__CP_t_" << t << ".txt" << std::ends;
      writeMatrixDLM<ScalarType>(oss.str().c_str(), TrajectoryPositions[t]);

      std::ostringstream oss2;
      oss2 << outputDir << Superclass::m_Name << "__MOM_t_" << t << ".txt" << std::ends;
      writeMatrixDLM<ScalarType>(oss2.str().c_str(), TrajectoryMomentas[t]);
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::WriteRegressionParameters() const {
  const std::string outputDir = def::utils::settings.output_dir;

  // Write control points.
  std::ostringstream oss1;
  oss1 << outputDir << Superclass::m_Name << "_ControlPoints.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss1.str().c_str(), GetControlPoints());

  // Write initialMomenta.
  std::ostringstream oss2;
  oss2 << outputDir << Superclass::m_Name << "_InitialMomenta.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss2.str().c_str(), GetInitialMomenta());
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::WriteTemplateData() const {
  const std::string outputDir = def::utils::settings.output_dir;

  // Writes template objects.
  std::vector<std::string> FullNames(m_TemplateObjectsName.size());
  for (int i = 0; i < m_TemplateObjectsName.size(); i++) {
    std::ostringstream oss;
    oss << outputDir << Superclass::m_Name << "_" << m_TemplateObjectsName[i] << m_TemplateObjectsNameExtension[i]
        << std::ends;
    FullNames[i] = oss.str();
  }

  GetTemplate()->WriteMultiObject(FullNames);
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::ComputeDataTermGradient(const MatrixType &initialMomenta,
                          const std::vector<std::shared_ptr<DeformableMultiObjectType>> targets,
                          std::vector<unsigned int> timeIndices,
                          MatrixType &gradPos, MatrixType &gradMom,
                          MatrixListType &gradTempL_L2, MatrixListType &gradTempL_Sob) {
  const unsigned int numberOfObjects = m_NumberOfObjects;

  if (targets.size() != timeIndices.size())
    throw std::runtime_error("Mismatch between number of observations and number of time indices.");

  gradTempL_L2 = GetTemplateData(); // to get the correct matrices size
  gradTempL_L2.fill(0.0);

  // A copy of the deformation is needed since this method can be called by different threads (maybe?) at the same time
  std::shared_ptr<DiffeosType> subjectDef = m_Def->Clone(); // shoot control points and initialMomenta and flow baseline

  subjectDef->SetDeformableMultiObject(GetTemplate());
  subjectDef->SetStartPositions(GetControlPoints());
  subjectDef->SetStartMomentas(initialMomenta);
  subjectDef->Update();

  if (subjectDef->OutOfBox())
    throw std::runtime_error("Out of box in Regression::ComputeDataTermGradient (this should not be...)");

  // Compute gradient of data matching term
  unsigned int numberOfObservations = targets.size();
  MatrixListType GradientDataTermOfLandmarkTypes(numberOfObservations);
  MatrixListType GradientDataTermOfImageTypes(numberOfObservations);
  for (unsigned int s = 0; s < numberOfObservations; s++) {
    std::shared_ptr<DeformableMultiObjectType> deformedTemplateObjects
        = subjectDef->GetDeformedObjectAt(timeIndices[s]);

    // Get the gradient of the similarity metric between deformed template and target
    MatrixListType currGradDataTi = deformedTemplateObjects->ComputeMatchGradient(targets[s]);

    for (unsigned int i = 0; i < numberOfObjects; i++)
      currGradDataTi[i] /= 2.0 * m_DataSigmaSquared[i];

    GetTemplate()->ListToMatrices(
        currGradDataTi, GradientDataTermOfLandmarkTypes[s], GradientDataTermOfImageTypes[s]);
  }

  subjectDef->IntegrateAdjointEquations(GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes, timeIndices);

  /// Get the gradient w.r.t. deformation parameters and landmark points positions
  gradPos = -subjectDef->GetAdjointPosAt0();
  gradMom = -subjectDef->GetAdjointMomAt0();

  /// If needed, computes the gradient w.r.t. the template.
  if (!m_FreezeTemplateFlag) {
    MatrixType dTempLLandmarkTypes = -subjectDef->GetAdjointLandmarkPointsAt0();
    MatrixType dTempLImageTypes = -subjectDef->SplatResidualImages(targets, timeIndices);
    if (m_Template->GetNumberOfImageKindObjects()) {
      dTempLImageTypes /= m_DataSigmaSquared(m_Template->GetImageIndex());
    }

    GetTemplate()->MatricesToList(gradTempL_L2, dTempLLandmarkTypes, dTempLImageTypes);
    gradTempL_Sob = ConvolveGradTemplate(gradTempL_L2);
  }
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::AddGradientRegularityTerm(MatrixType gradPos, MatrixType gradMom, const MatrixType &initialMomenta) {
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();

  std::shared_ptr<KernelType> momKernelObj = kfac->CreateKernelObject(this->m_Def->GetKernelType());
  momKernelObj->SetKernelWidth(this->m_Def->GetKernelWidth());
  momKernelObj->SetSources(GetControlPoints());

  momKernelObj->SetWeights(initialMomenta);
  gradMom -= momKernelObj->Convolve(GetControlPoints());
  gradPos -= momKernelObj->ConvolveGradient(GetControlPoints(), initialMomenta);
}

template<class ScalarType, unsigned int Dimension>
MatrixListType
Regression<ScalarType, Dimension>
::ConvolveGradTemplate(MatrixListType &gradTemplate_L2) {
  const unsigned int numberOfObjects = m_NumberOfObjects;

  if (m_SmoothingKernelWidth < 1e-20) {
    return gradTemplate_L2;
  }

  MatrixListType gradTemplate_Sob(numberOfObjects);

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> momKernelObj = kfac->CreateKernelObject(m_Def->GetKernelType());
  momKernelObj->SetKernelWidth(m_SmoothingKernelWidth);

  const MatrixListType templateData = GetTemplateData();
  for (unsigned int i = 0; i < numberOfObjects; i++) {
    if (GetTemplate()->GetObjectList()[i]->IsOfLandmarkKind()) {
      momKernelObj->SetSources(templateData[i]);
      momKernelObj->SetWeights(gradTemplate_L2[i]);

      gradTemplate_Sob[i] = momKernelObj->Convolve(templateData[i]);
    } else
      gradTemplate_Sob[i] = gradTemplate_L2[i];
  }

  return gradTemplate_Sob;
}

template<class ScalarType, unsigned int Dimension>
void
Regression<ScalarType, Dimension>
::UpdateDeformationAndKernelDataDomain(const std::vector<std::shared_ptr<DeformableMultiObjectType>> target) {
  // There are 3 boxes:
  // boundingBox: tighlty enclosing template data and control points
  // a padded version with padding equal to 0.5f*PaddingFactor*DeformationKernelWidth: points should never move outside this box, it used in Diffeos->CheckBoundingBox
  // a padded version with padding equal to      PaddingFactor*DeformationKernelWidth: used to define p3MKernel grid to account for periodic boundary conditions

  MatrixType dataDomain = m_BoundingBox;
  for (unsigned int t = 0; t < target.size(); ++t) {
    MatrixType BB = target[t]->GetBoundingBox();
    for (int d = 0; d < Dimension; d++) {
      dataDomain(d, 0) = (dataDomain(d, 0) < BB(d, 0) ? dataDomain(d, 0) : BB(d, 0));
      dataDomain(d, 1) = (dataDomain(d, 1) > BB(d, 1) ? dataDomain(d, 1) : BB(d, 1));
    }
  }

  m_Def->SetDataDomain(dataDomain); // will be used to define bounding box with padding factor of 0.5*PaddingFactor*DeformationKernelWidth

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetDataDomain(dataDomain);  // will be used to define bounding box with padding factor of    PaddingFactor*DeformationKernelWidth
}

template
class Regression<double, 2>;
template
class Regression<double, 3>;
