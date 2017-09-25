/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "AbstractAtlas.h"
#include "LinearAlgebra.h"

using namespace def::algebra;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class ScalarType, unsigned int Dimension>
AbstractAtlas<ScalarType, Dimension>
::AbstractAtlas() : Superclass(), m_Def(NULL), m_SmoothingKernelWidth(0.0), m_NumberOfThreads(1),
                    m_FreezeTemplateFlag(0), m_FreezeControlPointsFlag(0) {
  MatrixType controlPoints;
  MatrixListType tempData;
  Superclass::m_FixedEffects["ControlPoints"] = controlPoints;
  Superclass::m_FixedEffects["TemplateData"] = tempData;
}

template<class ScalarType, unsigned int Dimension>
AbstractAtlas<ScalarType, Dimension>
::~AbstractAtlas() {}

template<class ScalarType, unsigned int Dimension>
AbstractAtlas<ScalarType, Dimension>
::AbstractAtlas(const AbstractAtlas &other) : Superclass(other) {
  m_Def = std::static_pointer_cast<DiffeosType>(other.m_Def->Clone());
  m_SmoothingKernelWidth = other.m_SmoothingKernelWidth;
  m_NumberOfThreads = other.m_NumberOfThreads;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::SetFixedEffects(const LinearVariableMapType &map) {
  const MatrixListType templateData_memory = GetTemplateData();
  const MatrixType controlPoints_memory = GetControlPoints();

  Superclass::SetFixedEffects(map);

  if (m_FreezeTemplateFlag) Superclass::m_FixedEffects["TemplateData"] = templateData_memory;
  else {
    m_Template->UpdateImageIntensityAndLandmarkPointCoordinates(templateData_memory);
    m_Template->Update();
  }

  if (m_FreezeControlPointsFlag) SetControlPoints(controlPoints_memory);
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::SetControlPoints(const std::string &fn) {
  if (strlen(fn.c_str())) {
    const MatrixType controlPoints = readMatrixDLM<ScalarType>(fn.c_str());
    Superclass::m_FixedEffects["ControlPoints"] = controlPoints;
    std::cout << "Using a set of " << controlPoints.rows() << " control points in file " << fn << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::Update() {
  if (m_Def == NULL)
    throw std::runtime_error("A deformation should be set to the atlas model");

  /// Template.
  m_Template->Update();
  m_NumberOfObjects = m_Template->GetNumberOfObjects();
  m_BoundingBox = m_Template->GetBoundingBox();

  /// Control points.
  if (GetControlPoints().rows() == 0) {
    if (m_CPSpacing == 0.0) {
      m_CPSpacing = m_Def->GetKernelWidth();
      std::cout << "No initial CP spacing given: using diffeo kernel width of " << m_CPSpacing << std::endl;
    }

    InitializeControlPoints();
  }

  /// Bounding box.
  const MatrixType controlPoints = GetControlPoints();
  for (unsigned int i = 0; i < controlPoints.rows(); i++) {
    for (unsigned int dim = 0; dim < Dimension; dim++) {
      m_BoundingBox(dim, 0) =
          (m_BoundingBox(dim, 0) < controlPoints(i, dim) ? m_BoundingBox(dim, 0) : controlPoints(i, dim));
      m_BoundingBox(dim, 1) =
          (m_BoundingBox(dim, 1) > controlPoints(i, dim) ? m_BoundingBox(dim, 1) : controlPoints(i, dim));
    }
  }
}

template<class ScalarType, unsigned int Dimension>
bool
AbstractAtlas<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   LinearVariableMapType const &popRER,
                   LinearVariablesMapType const &indRER,
                   std::vector<std::vector<std::vector<ScalarType>>> &residuals) {
  std::cout << "Abstract atlas compute residuals" << std::endl;

  std::vector<std::vector<ScalarType>> aux;


  bool oob = ComputeResiduals(dataSet, popRER, indRER, aux);
  for (unsigned int i = 0; i < aux.size(); ++i) {
    residuals[i].resize(1);
    residuals[i][0] = aux[i];
  }

  return oob;
}

template<class ScalarType, unsigned int Dimension>
bool
AbstractAtlas<ScalarType, Dimension>
::ComputeResiduals(const MatrixType &controlPoints,
                   const std::vector<MatrixType> &momentas,
                   const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                   std::vector<std::vector<ScalarType>> &residuals) {
  assert(momentas.size() == target.size());

  const unsigned int numberOfSubjects = momentas.size();
  UpdateDeformationAndKernelDataDomain(target);

  //
  // Version without multi-threading :
  //
  if (m_NumberOfThreads < 2) {
    residuals.resize(numberOfSubjects);

    bool oob = false;
    for (int s = 0; s < numberOfSubjects; s++) {
      oob = this->ComputeResidualsSubject(controlPoints, momentas[s], target[s], residuals[s]);

      if (oob)
        return true;
    }

    return false;
  }

  //
  // Multi-threaded version :
  //
  m_MT_SubjectCounter = 0;
  m_MT_NumberOfSubjects = numberOfSubjects;

  m_MT_ControlPoints = controlPoints;
  m_MT_Momentas = momentas;

  if (m_MT_Target != target)
    m_MT_Target = target;

  m_MT_OutOfBox = false;

  m_MT_Residuals.resize(numberOfSubjects);
  for (int s = 0; s < numberOfSubjects; s++)
    m_MT_Residuals[s].resize(m_NumberOfObjects);

  itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

  //unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
  unsigned int numThreads = m_NumberOfThreads;
  if (numThreads > m_MT_NumberOfSubjects)
    numThreads = m_MT_NumberOfSubjects;
  if (numThreads < 2)
    numThreads = 2;

  threader->SetNumberOfThreads(numThreads);
  threader->SetSingleMethod(&AbstractAtlas::_residualsThread, (void *) this);
  threader->SingleMethodExecute();

//	Residuals = m_MT_Residuals;

  residuals.resize(numberOfSubjects);
  for (int s = 0; s < numberOfSubjects; s++)
    residuals[s] = m_MT_Residuals[s];

  return m_MT_OutOfBox;
}

template<class ScalarType, unsigned int Dimension>
bool
AbstractAtlas<ScalarType, Dimension>
::ComputeResidualsSubject(const MatrixType &controlPoints,
                          const MatrixType &momenta,
                          const std::shared_ptr<DeformableMultiObjectType> target,
                          std::vector<ScalarType> &residuals) const {
  if (momenta.rows() != controlPoints.rows())
    throw std::runtime_error("Number of Momentas and Control Points mismatch");

  // A copy of the deformation is needed since this method can be called by different threads at the same time
  std::shared_ptr<DiffeosType> subjectDef = m_Def->Clone();

  subjectDef->SetDeformableMultiObject(GetTemplate());
  subjectDef->SetStartPositions(controlPoints);
  subjectDef->SetStartMomentas(momenta);
  subjectDef->Update();

  if (subjectDef->OutOfBox())
    return true;

  residuals = subjectDef->GetDeformedObject()->ComputeMatch(target);

  // no normalization in the residuals
//	for (int i = 0; i < m_NumberOfObjects; i++)
//		Residuals[i] /= ( 2.0 * m_DataSigmaSquared(i) );

  return false;
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
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
AbstractAtlas<ScalarType, Dimension>
::Write(std::vector<MatrixType> const &momentas, MatrixType const &covarianceMomentaInverse,
        VectorType const &dataSigmaSquared, const LongitudinalDataSetType *const dataSet) const {
  if ((m_TemplateObjectsName.size() == 0) || (m_TemplateObjectsNameExtension.size()) == 0)
    throw std::runtime_error("No file names given");

  if ((m_TemplateObjectsName.size() != m_NumberOfObjects)
      || (m_TemplateObjectsNameExtension.size() != m_NumberOfObjects))
    throw std::runtime_error("Number of objects and objects' name mismatch in AbstractAtlas::WriteOutput");

  WriteAtlasToSubjectDeformations(momentas, dataSet);
  WriteAtlasParameters(momentas, covarianceMomentaInverse, dataSigmaSquared);
  WriteTemplateData();
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::WriteAtlasToSubjectDeformations(std::vector<MatrixType> const &momentas, const LongitudinalDataSetType *const dataSet) const {

  const MatrixType controlPoints = GetControlPoints();
  MatrixType initialMomenta;
  MatrixType residuals(momentas.size(), this->GetTemplate()->GetNumberOfObjects(), 0.);
  const CrossSectionalDataSetType *const crossSectionalDataSet = static_cast<const CrossSectionalDataSetType *const>(dataSet);
  std::vector<std::shared_ptr<DeformableMultiObjectType>> target = crossSectionalDataSet->GetDeformableMultiObjects();
  const std::string outputDir = def::utils::settings.output_dir;

  for (unsigned int s = 0; s < momentas.size(); s++) {
    initialMomenta = momentas[s];

    if (initialMomenta.rows() != controlPoints.rows())
      throw std::runtime_error("Number of Momentas and Control Points mismatch");

    m_Def->SetDeformableMultiObject(GetTemplate());
    m_Def->SetStartPositions(controlPoints);
    m_Def->SetStartMomentas(initialMomenta);
    m_Def->Update();

    if (m_Def->OutOfBox())
      throw std::runtime_error("Deformation out of box in AbstractAtlas::WriteAtlasDeformation: this should not be.");

    std::vector<std::string> names(m_TemplateObjectsName.size());
    for (unsigned int i = 0; i < m_TemplateObjectsName.size(); i++) {
      std::ostringstream oss;
      oss << outputDir << Superclass::m_Name << "_" << m_TemplateObjectsName[i] << "_to_subject_" << s;
      names[i] = oss.str();
    }

    m_Def->WriteFlow(names, m_TemplateObjectsNameExtension);

    std::shared_ptr<DeformableMultiObjectType> deformedTemplate = m_Def->GetDeformedObject();
    std::vector<ScalarType> aux = deformedTemplate->ComputeMatch(target[s]);
    for (unsigned int i=0 ; i < aux.size(); ++i)
      residuals(s,i) = aux[i];

    // Output deformation field.
    DeformationFieldIO<ScalarType, Dimension> defFieldIO;
    defFieldIO.SetDiffeos(m_Def);
    defFieldIO.SetAnatomicalCoordinateSystemLabel("LPS");
    defFieldIO.WriteDeformationField(outputDir + Superclass::m_Name, true);
  }

  writeMatrixDLM<ScalarType>(outputDir+"Residuals.txt", residuals);

}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::WriteAtlasParameters(std::vector<MatrixType> const &momentas, MatrixType const &covarianceMomentaInverse,
                       VectorType const &dataSigmaSquared) const {
  const std::string outputDir = def::utils::settings.output_dir;

  /// Write control points.
  std::ostringstream oss1;
  oss1 << outputDir << Superclass::m_Name << "_ControlPoints.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss1.str().c_str(), GetControlPoints());

  /// Write initial momentas.
  std::ostringstream oss2;
  oss2 << outputDir << Superclass::m_Name << "_Momentas.txt" << std::ends;
  writeMultipleMatrixDLM<ScalarType>(oss2.str().c_str(), momentas);

  /// Write inverse of optimal covariance momenta matrix.
  std::ostringstream ossCMI;
  ossCMI << outputDir << this->m_Name << "_CovarianceMomentaInverse.txt" << std::ends;
  writeMatrixDLM<ScalarType>(ossCMI.str().c_str(), covarianceMomentaInverse);

  /// Write data sigma.
  MatrixType dataSigma(1, dataSigmaSquared.size());
  for (int i = 0; i < dataSigmaSquared.size(); i++)
    dataSigma(0, i) = sqrt(dataSigmaSquared(i));

  std::ostringstream ossDSS;
  ossDSS << outputDir << this->m_Name << "_DataSigma.txt" << std::ends;
  writeMatrixDLM<ScalarType>(ossDSS.str().c_str(), dataSigma);
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
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
AbstractAtlas<ScalarType, Dimension>
::ComputeDataTermGradient(const MatrixType &controlPoints,
                          const std::vector<MatrixType> &momentas,
                          const VectorType &dataSigmaSquared,
                          const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                          MatrixType &gradPos, std::vector<MatrixType> &gradMom, MatrixListType &gradTempL_L2) {

  if (momentas.size() != target.size())
    throw std::runtime_error("Number of subjects in momentas and target lists mismatch");

  int nbSubjects = momentas.size();
  int nbObjects = m_NumberOfObjects;

  gradPos.set_size(controlPoints.rows(), Dimension);
  gradPos.fill(0.0);
  gradMom.resize(nbSubjects);
  gradTempL_L2 = this->GetTemplateData(); // to get the correct matrices size
  for (int i = 0; i < nbObjects; i++)
    gradTempL_L2[i].fill(0.0);

  /// Version without multi-threading.
  if (m_NumberOfThreads < 2) {
    for (int s = 0; s < nbSubjects; s++) {
      MatrixType dPos;
      MatrixType dMom;
      MatrixListType dTempL(nbObjects);
      ComputeDataTermGradientSubject(controlPoints, momentas[s], dataSigmaSquared, target[s],
                                     dPos, dMom, dTempL);

      gradPos += dPos;
      gradMom[s] = dMom;
      for (unsigned int i = 0; i < nbObjects; i++) { gradTempL_L2[i] += dTempL[i]; }
    }
    return;
  }

  /// Multi-threading.
  m_MT_SubjectCounter = 0;
  m_MT_NumberOfSubjects = nbSubjects;

  m_MT_ControlPoints = controlPoints;
  m_MT_Momentas = momentas;
  m_MT_DataSigmaSquared = dataSigmaSquared;

  if (m_MT_Target != target) { m_MT_Target = target; }

  m_MT_GradPos = gradPos;
  m_MT_GradMom = gradMom;
  m_MT_GradTempL_L2 = gradTempL_L2;

  itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

  //unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
  unsigned int numThreads = m_NumberOfThreads;
  if (numThreads > m_MT_NumberOfSubjects)
    numThreads = m_MT_NumberOfSubjects;
  if (numThreads < 2)
    numThreads = 2;

  threader->SetNumberOfThreads(numThreads);
  threader->SetSingleMethod(&AbstractAtlas::_gradientResidualsThread, (void *) this);
  threader->SingleMethodExecute();

  gradPos = m_MT_GradPos;
  gradMom = m_MT_GradMom;
  gradTempL_L2 = m_MT_GradTempL_L2;
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::ComputeDataTermGradientSubject(const MatrixType &controlPoints,
                                 const MatrixType &momenta,
                                 const VectorType &dataSigmaSquared,
                                 const std::shared_ptr<DeformableMultiObjectType> target,
                                 MatrixType &dPos, MatrixType &dMom, MatrixListType &dTempL) {

  if (momenta.rows() != controlPoints.rows())
    throw std::runtime_error("Number of momenta and Control Points mismatch");

  // A copy of the deformation is needed since this method can be called by different threads at the same time
  std::shared_ptr<DiffeosType> subjectDef = m_Def->Clone();
  subjectDef->SetDeformableMultiObject(GetTemplate());
  subjectDef->SetStartPositions(controlPoints);
  subjectDef->SetStartMomentas(momenta);
  subjectDef->Update();

  if (subjectDef->OutOfBox())
    throw std::runtime_error("Out of box in AbstractAtlas::ComputeResidualGradient (this should not be).");

  /// Get the deformed template
  std::shared_ptr<DeformableMultiObjectType> deformedTemplateObjects = subjectDef->GetDeformedObject();

  /// Get the gradient of the similarity metric between deformed template and target
  MatrixListType GradientSimilarityMetric = deformedTemplateObjects->ComputeMatchGradient(target);

  /// Divide each gradient of the data term by 1/(2*DataSigmaSquared)
  for (int i = 0; i < m_NumberOfObjects; i++)
    GradientSimilarityMetric[i] /= (2.0f * dataSigmaSquared(i));

  /// Pool the gradients into big matrices for all objects of type landmarks (and children types) and image (and children types).
  MatrixType GradientDataTermOfLandmarkTypes;
  MatrixType GradientDataTermOfImageTypes;
  GetTemplate()->ListToMatrices(GradientSimilarityMetric, GradientDataTermOfLandmarkTypes,
                                GradientDataTermOfImageTypes);

  /// Integrate the adjoint equations.
  subjectDef->IntegrateAdjointEquations(GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes);

  /// Get the gradient w.r.t. deformation parameters and landmark points positions.
  dPos = -subjectDef->GetAdjointPosAt0();
  dMom = -subjectDef->GetAdjointMomAt0();

  /// If needed, computes the gradient w.r.t. the template.
  if (!m_FreezeTemplateFlag) {
    MatrixType dTempLLandmarkTypes = -subjectDef->GetAdjointLandmarkPointsAt0();
    MatrixType dTempLImageTypes = -subjectDef->SplatResidualImage(target);
    if (m_Template->GetNumberOfImageKindObjects()) {
      dTempLImageTypes /= dataSigmaSquared(m_Template->GetImageIndex());
    }

    GetTemplate()->MatricesToList(dTempL, dTempLLandmarkTypes, dTempLImageTypes);
  }
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::ComputeDataTermGradient(const MatrixType &controlPoints,
                          const std::vector<MatrixType> &momentas,
                          const VectorType &dataSigmaSquared,
                          const std::vector<std::shared_ptr<DeformableMultiObjectType>> target,
                          MatrixType &gradPos, std::vector<MatrixType> &gradMom) {

  if (momentas.size() != target.size())
    throw std::runtime_error("Number of subjects in momentas and target lists mismatch");

  int nbSubjects = momentas.size();
  int nbObjects = m_NumberOfObjects;

  gradPos.set_size(controlPoints.rows(), Dimension);
  gradPos.fill(0.0);
  gradMom.resize(nbSubjects);

  /// Version without multi-threading.
  if (m_NumberOfThreads < 2) {
    for (int s = 0; s < nbSubjects; s++) {
      MatrixType dPos;
      MatrixType dMom;
      MatrixListType dTempL(nbObjects);
      ComputeDataTermGradientSubject(controlPoints, momentas[s], dataSigmaSquared, target[s], dPos, dMom);

      gradPos += dPos;
      gradMom[s] = dMom;
    }
    return;
  }

  /// Multi-threading.
  m_MT_SubjectCounter = 0;
  m_MT_NumberOfSubjects = nbSubjects;

  m_MT_ControlPoints = controlPoints;
  m_MT_Momentas = momentas;
  m_MT_DataSigmaSquared = dataSigmaSquared;

  if (m_MT_Target != target)
    m_MT_Target = target;

  m_MT_GradPos = gradPos;
  m_MT_GradMom = gradMom;

  itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

  //unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
  unsigned int numThreads = m_NumberOfThreads;
  if (numThreads > m_MT_NumberOfSubjects)
    numThreads = m_MT_NumberOfSubjects;
  if (numThreads < 2)
    numThreads = 2;

  threader->SetNumberOfThreads(numThreads);
  threader->SetSingleMethod(&AbstractAtlas::_gradientWithoutTemplateResidualsThread, (void *) this);
  threader->SingleMethodExecute();

  gradPos = m_MT_GradPos;
  gradMom = m_MT_GradMom;
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::ComputeDataTermGradientSubject(const MatrixType &controlPoints,
                                 const MatrixType &momenta,
                                 const VectorType &dataSigmaSquared,
                                 const std::shared_ptr<DeformableMultiObjectType> target,
                                 MatrixType &dPos, MatrixType &dMom) {
  assert(momenta.rows() == controlPoints.rows());

  // A copy of the deformation is needed since this method can be called by different threads at the same time
  std::shared_ptr<DiffeosType> subjectDef = m_Def->Clone();
  subjectDef->SetDeformableMultiObject(GetTemplate());
  subjectDef->SetStartPositions(controlPoints);
  subjectDef->SetStartMomentas(momenta);
  subjectDef->Update();

  if (subjectDef->OutOfBox())
    throw std::runtime_error("Out of box in AbstractAtlas::ComputeResidualGradient (this should not be).");

  /// Get the deformed template
  std::shared_ptr<DeformableMultiObjectType> deformedTemplateObjects = subjectDef->GetDeformedObject();

  /// Get the gradient of the similarity metric between deformed template and target
  MatrixListType GradientSimilarityMetric = deformedTemplateObjects->ComputeMatchGradient(target);

  /// Divide each gradient of the data term by 1/(2*DataSigmaSquared)
  for (int i = 0; i < m_NumberOfObjects; i++)
    GradientSimilarityMetric[i] /= (2.0 * dataSigmaSquared(i));

  /// Pool the gradients into big matrices for all objects of type landmarks (and children types) and image (and children types).
  MatrixType GradientDataTermOfLandmarkTypes;
  MatrixType GradientDataTermOfImageTypes;
  GetTemplate()->ListToMatrices(GradientSimilarityMetric, GradientDataTermOfLandmarkTypes,
                                GradientDataTermOfImageTypes);

  /// Integrate the adjoint equations.
  subjectDef->IntegrateAdjointEquations(GradientDataTermOfLandmarkTypes, GradientDataTermOfImageTypes);

  /// Get the gradient w.r.t. deformation parameters and landmark points positions.
  dPos = -subjectDef->GetAdjointPosAt0();
  dMom = -subjectDef->GetAdjointMomAt0();
}

template<class ScalarType, unsigned int Dimension>
MatrixListType
AbstractAtlas<ScalarType, Dimension>
::ConvolveGradTemplate(MatrixListType &GradTemplate_L2) {

  if (m_SmoothingKernelWidth < 1e-20) {
    return GradTemplate_L2;
  }

  MatrixListType GradTemplate_Sob(m_NumberOfObjects);

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> momKernelObj = kfac->CreateKernelObject(m_Def->GetKernelType());
  momKernelObj->SetKernelWidth(m_SmoothingKernelWidth);

  const MatrixListType tempData = GetTemplateData();
  for (unsigned int i = 0; i < m_NumberOfObjects; i++) {
    if (GetTemplate()->GetObjectList()[i]->IsOfLandmarkKind()) {
      momKernelObj->SetSources(tempData[i]);
      momKernelObj->SetWeights(GradTemplate_L2[i]);

      GradTemplate_Sob[i] = momKernelObj->Convolve(tempData[i]);
    } else
      GradTemplate_Sob[i] = GradTemplate_L2[i];
  }

  return GradTemplate_Sob;
}

template<class ScalarType, unsigned int Dimension>
void
AbstractAtlas<ScalarType, Dimension>
::UpdateDeformationAndKernelDataDomain(const std::vector<std::shared_ptr<DeformableMultiObjectType>> target) {
  // There are 3 boxes:
  // Superclass::m_FixedEffects->GetBoundingBox() : tightly enclosing template data and control points
  // a padded version with padding equal to 0.5f*PaddingFactor*DeformationKernelWidth: points should never move outside this box, it used in Diffeos->CheckBoundingBox
  // a padded version with padding equal to     PaddingFactor*DeformationKernelWidth: used to define p3MKernel grid to account for periodic boundary conditions

  MatrixType DataDomain = m_BoundingBox;
  for (unsigned int s = 0; s < target.size(); s++) {
    MatrixType BB = target[s]->GetBoundingBox();
    for (int d = 0; d < Dimension; d++) {
      DataDomain(d, 0) = (DataDomain(d, 0) < BB(d, 0) ? DataDomain(d, 0) : BB(d, 0));
      DataDomain(d, 1) = (DataDomain(d, 1) > BB(d, 1) ? DataDomain(d, 1) : BB(d, 1));
    }
  }

  m_Def->SetDataDomain(DataDomain); // will be used to define bounding box with padding factor of 0.5*PaddingFactor*DeformationKernelWidth

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetDataDomain(DataDomain);  // will be used to define bounding box with padding factor of    PaddingFactor*DeformationKernelWidth
}


//
// Multi-threaded components
//

template<class ScalarType, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
AbstractAtlas<ScalarType, Dimension>
::_residualsThread(void *arg) {
  typedef itk::MultiThreader::ThreadInfoStruct ThreadInfoType;
  ThreadInfoType *infoStruct = static_cast<ThreadInfoType *>(arg);
  AbstractAtlas *obj = static_cast<AbstractAtlas *>(infoStruct->UserData);

  unsigned int s = obj->m_MT_NumberOfSubjects;

  while (true) {
    obj->m_Mutex.Lock();
    s = obj->m_MT_SubjectCounter++;
    obj->m_Mutex.Unlock();

    if (s >= obj->m_MT_NumberOfSubjects)
      break;

    std::vector<ScalarType> values;

    bool oob = obj->ComputeResidualsSubject(obj->m_MT_ControlPoints, obj->m_MT_Momentas[s],
                                            obj->m_MT_Target[s], values);

    obj->m_Mutex.Lock();
    obj->m_MT_OutOfBox = oob;
    obj->m_Mutex.Unlock();

    if (oob)
      break;

    for (int i = 0; i < obj->m_NumberOfObjects; i++) {
      obj->m_Mutex.Lock();
      obj->m_MT_Residuals[s][i] = values[i];
      obj->m_Mutex.Unlock();
    }

  }

  return ITK_THREAD_RETURN_VALUE;
}

template<class ScalarType, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
AbstractAtlas<ScalarType, Dimension>
::_gradientResidualsThread(void *arg) {
  typedef itk::MultiThreader::ThreadInfoStruct ThreadInfoType;
  ThreadInfoType *infoStruct = static_cast<ThreadInfoType *>( arg );
  AbstractAtlas *obj = static_cast<AbstractAtlas *>(infoStruct->UserData);

  unsigned int s = obj->m_MT_NumberOfSubjects;

  while (true) {
    obj->m_Mutex.Lock();
    s = obj->m_MT_SubjectCounter++;
    obj->m_Mutex.Unlock();

    if (s >= obj->m_MT_NumberOfSubjects)
      break;

    MatrixType dPos;
    MatrixType dMom;
    MatrixListType dTempL(obj->m_NumberOfObjects);

    obj->ComputeDataTermGradientSubject(obj->m_MT_ControlPoints, obj->m_MT_Momentas[s],
                                        obj->m_MT_DataSigmaSquared, obj->m_MT_Target[s],
                                        dPos, dMom, dTempL);

    obj->m_Mutex.Lock();
    obj->m_MT_GradPos += dPos;
    obj->m_Mutex.Unlock();

    obj->m_Mutex.Lock();
    obj->m_MT_GradMom[s] = dMom;
    obj->m_Mutex.Unlock();

    for (unsigned int i = 0; i < obj->m_NumberOfObjects; i++) {
      obj->m_Mutex.Lock();
      obj->m_MT_GradTempL_L2[i] += dTempL[i];
      obj->m_Mutex.Unlock();
    }
  }

  return ITK_THREAD_RETURN_VALUE;
}

template<class ScalarType, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
AbstractAtlas<ScalarType, Dimension>
::_gradientWithoutTemplateResidualsThread(void *arg) {
  typedef itk::MultiThreader::ThreadInfoStruct ThreadInfoType;
  ThreadInfoType *infoStruct = static_cast<ThreadInfoType *>( arg );
  AbstractAtlas *obj = static_cast<AbstractAtlas *>(infoStruct->UserData);

  unsigned int s = obj->m_MT_NumberOfSubjects;

  while (true) {
    obj->m_Mutex.Lock();
    s = obj->m_MT_SubjectCounter++;
    obj->m_Mutex.Unlock();

    if (s >= obj->m_MT_NumberOfSubjects)
      break;

    MatrixType dPos;
    MatrixType dMom;

    obj->ComputeDataTermGradientSubject(obj->m_MT_ControlPoints, obj->m_MT_Momentas[s],
                                        obj->m_MT_DataSigmaSquared, obj->m_MT_Target[s],
                                        dPos, dMom);

    obj->m_Mutex.Lock();
    obj->m_MT_GradPos += dPos;
    obj->m_Mutex.Unlock();

    obj->m_Mutex.Lock();
    obj->m_MT_GradMom[s] = dMom;
    obj->m_Mutex.Unlock();
  }

  return ITK_THREAD_RETURN_VALUE;
}

template
class AbstractAtlas<double, 2>;
template
class AbstractAtlas<double, 3>;

