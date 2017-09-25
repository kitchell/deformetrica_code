/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "LongitudinalRegistration.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
LongitudinalRegistration<ScalarType, Dimension>
::LongitudinalRegistration() : Superclass(),
                           m_Template(NULL),
                           m_ForwardReferenceGeodesic(NULL),
                           m_BackwardReferenceGeodesic(NULL),
                           m_PerpendicularDeformation(NULL),
                           m_MarginOnGeodesicLength(2.0) {

  Superclass::SetLongitudinalMatchingType();

  /// Fixed effects.
  Superclass::m_FixedEffects["TimeShift"] = 0.0;
  Superclass::m_FixedEffects["LogAcceleration"] = 0.0;
  Superclass::m_FixedEffects["Sources"] = VectorType(1, 0.0);

  /// Scalar hyperparameters.
  m_ReferenceTime = std::sqrt(-1);
  m_TimeShiftVariance = 0.0;
  m_LogAccelerationVariance = 0.0;
}

template<class ScalarType, unsigned int Dimension>
LongitudinalRegistration<ScalarType, Dimension>
::~LongitudinalRegistration() {}

template<class ScalarType, unsigned int Dimension>
LongitudinalRegistration<ScalarType, Dimension>
::LongitudinalRegistration(const LongitudinalRegistration &other) : Superclass(other) {
  m_Template = other.m_Template->Clone();
  m_TemplateObjectsName = other.m_TemplateObjectsName;
  m_TemplateObjectsNameExtension = other.m_TemplateObjectsNameExtension;
  m_BoundingBox = other.m_BoundingBox;

  m_ControlPoints = other.m_ControlPoints;
  m_Momenta = other.m_Momenta;
  m_ModulationMatrix = other.m_ModulationMatrix;
  m_ProjectedModulationMatrix = other.m_ProjectedModulationMatrix;
  m_ReferenceTime = other.m_ReferenceTime;

  m_TimeShiftVariance = other.m_TimeShiftVariance;
  m_LogAccelerationVariance = other.m_LogAccelerationVariance;
  m_NoiseVariance = other.m_NoiseVariance;

  m_ConcentrationOfTimePointsForReferenceGeodesic = other.m_ConcentrationOfTimePointsForReferenceGeodesic;
  m_NumberOfTimePointsForExponentiation = other.m_NumberOfTimePointsForExponentiation;

  m_NumberOfObjects = other.m_NumberOfObjects;
  m_DimensionOfDiscretizedObjects = other.m_DimensionOfDiscretizedObjects;
  m_NumberOfControlPoints = other.m_NumberOfControlPoints;
  m_DimensionOfTangentSpaces = other.m_DimensionOfTangentSpaces;
  m_NumberOfSources = other.m_NumberOfSources;

  m_AbsoluteTimeIncrements = other.m_AbsoluteTimeIncrements;
  m_MinimumAbsoluteTimeIncrement = other.m_MinimumAbsoluteTimeIncrement;
  m_MaximumAbsoluteTimeIncrement = other.m_MaximumAbsoluteTimeIncrement;

  m_ForwardReferenceGeodesic = other.m_ForwardReferenceGeodesic->Clone();
  m_BackwardReferenceGeodesic = other.m_BackwardReferenceGeodesic->Clone();
  m_PerpendicularDeformation = other.m_PerpendicularDeformation->Clone();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::SetControlPoints(const std::string &fn) {
  if (strlen(fn.c_str())) {
    m_ControlPoints = readMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using a set of " << m_ControlPoints.rows() << " control points in file " << fn << std::endl;
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::SetMomenta(const std::string &fn) {
  if (strlen(fn.c_str())) {
    m_Momenta = readMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using a momenta matrix of size " << m_Momenta.rows() << " x "
              << m_Momenta.columns() << " from file: " << fn << std::endl;
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::SetModulationMatrix(const std::string &fn) {
  if (strlen(fn.c_str())) {
    m_ModulationMatrix = readMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using a modulation matrix with " << m_ModulationMatrix.cols() << " sources "
              << " from file: " << fn << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::Update() {
  /// Working deformation.
  if (m_PerpendicularDeformation == NULL) {
    throw std::runtime_error("A deformation should be set to the LongitudinalRegistration model.");
  }

  /// Template.
  if (m_Template == 0) { throw std::runtime_error("A template should be set to the LongitudinalRegistration model."); }
  m_Template->Update();
  m_NumberOfObjects = m_Template->GetNumberOfObjects();
  m_DimensionOfDiscretizedObjects = m_Template->GetDimensionOfDiscretizedObjects();
  InitializeBoundingBox();

  /// Control points.
  if (m_ControlPoints.rows() == 0) {
    throw std::runtime_error("Control points must be specified to the LongitudinalRegistration model.");
  }
  m_NumberOfControlPoints = m_ControlPoints.rows();
  m_DimensionOfTangentSpaces = m_ControlPoints.n_elem();

  /// Momenta.
  if ((m_Momenta.rows() != m_ControlPoints.rows()) || (m_Momenta.cols() != m_ControlPoints.cols())) {
    throw std::runtime_error("Valid momenta must be specified to the LongitudinalRegistration model.");
  }

  /// Modulation matrix.
  ProjectModulationMatrix();
  m_NumberOfSources = m_ModulationMatrix.cols();

  /// Reference time.
  if (m_ReferenceTime != m_ReferenceTime) { // i.e. isnan.
    throw std::runtime_error("A reference time must be specified to the LongitudinalRegistration model.");
  }

  /// Time-shift variance.
  if (m_TimeShiftVariance == 0) {
    throw std::runtime_error("A time-shift variance must be specified to the LongitudinalRegistration model.");
  }

  /// Log-acceleration variance.
  if (m_LogAccelerationVariance == 0) {
    throw std::runtime_error("A log-acceleration variance must be specified to the LongitudinalRegistration model.");
  }

  /// Noise variance.
  if (m_NoiseVariance.size() == 0) {
    throw std::runtime_error("A noise variance must be specified to the LongitudinalRegistration model.");
  }

  /// Fixed effects.
  Superclass::m_FixedEffects["TimeShift"] = m_ReferenceTime;
  Superclass::m_FixedEffects["Sources"] = VectorType(m_NumberOfSources, 0.0);

  /// Initializes the reference geodesic.
  UpdateReferenceGeodesic();
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::Print() const {
  std::cout << "[ model print - TODO ]" << std::endl;

//  const ScalarType referenceTime = Superclass::m_FixedEffects["ReferenceTime"].vectorize()[0];
//  const ScalarType timeShiftVariance = Superclass::m_FixedEffects["TimeShiftVariance"].vectorize()[0];
//  const ScalarType logAccelerationVariance = Superclass::m_FixedEffects["LogAccelerationVariance"].vectorize()[0];
//  const VectorType noiseVariance = Superclass::m_FixedEffects["NoiseVariance"].vectorize();
//
//  std::cout << ">> Model fixed effects :" << std::endl;
//  std::cout << "\t\tReference time = " << referenceTime << std::endl;
//  std::cout << "\t\tTime-shift standard deviation = " << std::sqrt(timeShiftVariance) << std::endl;
//  std::cout << "\t\tLog-acceleration standard deviation = " << std::sqrt(logAccelerationVariance) << std::endl;
//  for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
//    std::cout << "\t\tNoise standard deviation ";
//    if (m_NumberOfObjects > 1) { std::cout << "(" << m_TemplateObjectsName[k] << ")"; }
//    std::cout << " = " << std::sqrt(noiseVariance[k]) << std::endl;
//  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::Write(const LongitudinalDataSetType *const dataSet,
        LinearVariableMapType const &popRER,
        LinearVariablesMapType const &indRER) const {
  WriteModelPredictions(dataSet);
  WriteModelParameters();
}

template<class ScalarType, unsigned int Dimension>
bool
LongitudinalRegistration<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   std::vector<std::vector<ScalarType>> &residuals) {
  /// Extract needed variables.
  const unsigned int numberOfSubjects = dataSet->GetNumberOfSubjects();
  const unsigned int nbObservations = m_AbsoluteTimeIncrements.size();

  const TimeSeriesDataSetType *const timeSeriesDataSet = static_cast<const TimeSeriesDataSetType *const>(dataSet);
  const std::vector<std::shared_ptr<DeformableMultiObjectType>>
      targets = timeSeriesDataSet->GetDeformableMultiObjects();

  const VectorType sources = recast<VectorType>(Superclass::m_FixedEffects.at("Sources"));

  /// Transport long the reference geodesic.
  MatrixType spaceShift = (m_ProjectedModulationMatrix * sources).unvectorize(m_NumberOfControlPoints, Dimension);

  // Divide the target absolute times into backward and forward targets.
  std::vector<ScalarType> forwardAbsoluteTimeIncrements, backwardAbsoluteTimeIncrements;
  for (unsigned int t = 0; t < nbObservations; ++t) {
    if (m_AbsoluteTimeIncrements[t] < 0.0) {
      backwardAbsoluteTimeIncrements.push_back(-m_AbsoluteTimeIncrements[t]);
    } else { forwardAbsoluteTimeIncrements.push_back(m_AbsoluteTimeIncrements[t]); }
  }
  std::reverse(backwardAbsoluteTimeIncrements.begin(), backwardAbsoluteTimeIncrements.end());

  // Perform the transport.
  MatrixListType backwardTransportedSpaceShifts = m_BackwardReferenceGeodesic->ParallelTransport(
      spaceShift, backwardAbsoluteTimeIncrements).reverse();
  MatrixListType forwardTransportedSpaceShifts = m_ForwardReferenceGeodesic->ParallelTransport(
      spaceShift, forwardAbsoluteTimeIncrements);

  // Concatenate the results.
  MatrixListType transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
  assert(transportedSpaceShifts.size() == nbObservations);

  /// Exponentiation : shoot at each target time-point. Then compute the residual.
  residuals.resize(nbObservations);
  {
    ThreadPool pool(def::utils::settings.number_of_threads);
    for (unsigned int t = 0; t < nbObservations; ++t) {
      pool.enqueue([&, t]() {
        std::shared_ptr<DeformableMultiObjectType> referenceShape;
        MatrixType referenceControlPoints;
        GetDeformedObjectAndControlPointsAt(m_AbsoluteTimeIncrements[t], referenceShape, referenceControlPoints);

        std::shared_ptr<DiffeosType> expDef = m_PerpendicularDeformation->Clone();
        expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
        expDef->SetDeformableMultiObject(referenceShape);
        expDef->SetStartPositions(referenceControlPoints);
        expDef->SetStartMomentas(transportedSpaceShifts[t]);
        expDef->Update();

        std::shared_ptr<DeformableMultiObjectType> predictedShape = expDef->GetDeformedObject();
        residuals[t] = predictedShape->ComputeMatch(targets[t]);
      });
    }
  }
  return false; // TODO : return true if an out-of-box is detected.
}

template<class ScalarType, unsigned int Dimension>
ScalarType
LongitudinalRegistration<ScalarType, Dimension>
::ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                               const LinearVariableMapType &popRER,
                               const LinearVariablesMapType &indRER) {
  /// Initialization.
  UpdateAbsoluteTimeIncrements(dataSet);
  if ((m_MaximumAbsoluteTimeIncrement > 0 &&
      m_MaximumAbsoluteTimeIncrement > m_ForwardReferenceGeodesic->GetTN()) ||
      (m_MinimumAbsoluteTimeIncrement < 0 &&
          m_MinimumAbsoluteTimeIncrement < (-m_BackwardReferenceGeodesic->GetTN()))) {
    UpdateReferenceGeodesic(); // SHOULD ONLY BE EXTENDED. TODO.
  }

  /// Data (residuals) term.
  std::vector<std::vector<ScalarType>> residuals;
  ComputeResiduals(dataSet, residuals);

  ScalarType out = 0.0;
  for (unsigned int t = 0; t < residuals.size(); ++t) {
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      out -= 0.5 * (residuals[t][k] / m_NoiseVariance(k));
    }
  }

  /// Regularization terms. Actually prior terms, but not coded with the dedicated formalism here, for simplicity.
  // Time-shift.
  const ScalarType timeShift = recast<ScalarType>(Superclass::m_FixedEffects.at("TimeShift"));
  out -= 0.5 * pow(timeShift - m_ReferenceTime, 2) / m_TimeShiftVariance;
  // Log-acceleration.
  out -= 0.5 * Superclass::m_FixedEffects.at("LogAcceleration").sum_of_squares() / m_LogAccelerationVariance;
  // Source.
  out -= 0.5 * Superclass::m_FixedEffects.at("TimeShift").sum_of_squares();

  return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::WriteModelPredictions(const LongitudinalDataSetType *const dataSet) const {
  const std::string outputDir = def::utils::settings.output_dir;

  /// Remove pre-existing output files.
  namespace bfs=boost::filesystem;
  bfs::path path_to_remove(outputDir);
  for (bfs::directory_iterator end_dir_it, it(path_to_remove); it != end_dir_it; ++it) {
    if (it->path().string()[outputDir.size()] != '.' &&
        it->path().string().find(def::utils::settings.output_state_filename) == std::string::npos) {
      bfs::remove_all(it->path());
    }
  }

  /// Get the needed fixed effects.
  const ScalarType timeShift = recast<ScalarType>(Superclass::m_FixedEffects.at("TimeShift"));
  const ScalarType logAcceleration = recast<ScalarType>(Superclass::m_FixedEffects.at("LogAcceleration"));
  const VectorType sources = recast<VectorType>(Superclass::m_FixedEffects.at("Sources"));

  /// Compute the absolute times (cf. this->UpdateAbsoluteTimeIncrements, here a const version).
  const TimeSeriesDataSetType *const timeSeriesDataSet = static_cast<const TimeSeriesDataSetType *const>(dataSet);
  const std::vector<ScalarType> times = timeSeriesDataSet->GetTimes();
  const ScalarType acceleration = std::exp(logAcceleration);
  const unsigned long nbObservations = times.size();

  std::vector<ScalarType> absoluteTimeIncrements(nbObservations);
  for (unsigned int t = 0; t < nbObservations; ++t) {
    absoluteTimeIncrements[t] = acceleration * (times[t] - timeShift);
  }
  const ScalarType minimumAbsoluteTimeIncrement = absoluteTimeIncrements[0];
  const ScalarType maximumAbsoluteTimeIncrement = absoluteTimeIncrements[nbObservations - 1];

  /// Compute the reference geodesic (cf. this->UpdateReferenceGeodesic, here a const version).
  // Forward part.
  std::shared_ptr<DiffeosType> forwardDef = m_ForwardReferenceGeodesic->Clone();
  const ScalarType forwardGeodesicLength = maximumAbsoluteTimeIncrement;
  const int forwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * forwardGeodesicLength + 1.5;
  forwardDef->SetT0(0.0);
  forwardDef->SetTN(forwardGeodesicLength);
  forwardDef->SetNumberOfTimePoints(forwardNumberOfTimePoints);
  forwardDef->SetDeformableMultiObject(m_Template);
  forwardDef->SetStartPositions(m_ControlPoints);
  forwardDef->SetStartMomentas(m_Momenta);
  forwardDef->Update();

  // Backward part.
  std::shared_ptr<DiffeosType> backwardDef = m_BackwardReferenceGeodesic->Clone();
  const ScalarType backwardGeodesicLength = -minimumAbsoluteTimeIncrement;
  const int backwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * backwardGeodesicLength + 1.5;
  backwardDef->SetT0(0.0);
  backwardDef->SetTN(backwardGeodesicLength);
  backwardDef->SetNumberOfTimePoints(backwardNumberOfTimePoints);
  backwardDef->SetDeformableMultiObject(m_Template);
  backwardDef->SetStartPositions(m_ControlPoints);
  backwardDef->SetStartMomentas(-m_Momenta);
  backwardDef->Update();

  // Auxiliary variables.
  unsigned int clampedForwardNumberOfTimePoints = 0;
  if (forwardNumberOfTimePoints > 0) { clampedForwardNumberOfTimePoints = forwardNumberOfTimePoints; }
  unsigned int clampedBackwardNumberOfTimePoints = 1;
  if (backwardNumberOfTimePoints > 1) { clampedBackwardNumberOfTimePoints = backwardNumberOfTimePoints; }

  /// Write the reference geodesic, and optionally the transported control points and momenta.
  // Optional initialization.
  const bool writeFullTrajectoriesFlag = false;
  MatrixListType forwardPositions, forwardMomentas, backwardPositions, backwardMomenta;
  if (writeFullTrajectoriesFlag) {
    forwardPositions = forwardDef->GetTrajectoryPositions();
    forwardMomentas = forwardDef->GetTrajectoryMomentas();
    backwardPositions = backwardDef->GetTrajectoryPositions();
    backwardMomenta = -backwardDef->GetTrajectoryMomentas();
  }

  // Forward part.
  const ScalarType forwardStepSize = forwardGeodesicLength / (forwardNumberOfTimePoints - 1);
  for (unsigned int t = 0; t < clampedForwardNumberOfTimePoints; ++t) {
    const ScalarType time = m_ReferenceTime + t * forwardStepSize;
    const std::string roundedTime = std::to_string((unsigned int) (time * 100.0 + 0.5));
    const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
    const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

    std::vector<std::string> filenames(m_NumberOfObjects);
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      std::ostringstream filename;
      filename << outputDir << Superclass::m_Name << "_";
      if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
      filename << "_ReferenceGeodesic__tp_" << backwardNumberOfTimePoints + t - 1 << "__age_";
      if (time < 0) { filename << "-"; }
      filename << floor << "." << decimals << m_TemplateObjectsNameExtension[k];
      filenames[k] = filename.str();

      if (writeFullTrajectoriesFlag) {
        std::ostringstream oss1, oss2;
        oss1 << outputDir << Superclass::m_Name;
        oss2 << outputDir << Superclass::m_Name;
        oss1 << "__ReferenceGeodesic__CP__tp_" << backwardNumberOfTimePoints + t - 1 << "__age_";
        oss2 << "__ReferenceGeodesic__MOM__tp_" << backwardNumberOfTimePoints + t - 1 << "__age_";
        if (time < 0) {
          oss1 << "-";
          oss2 << "-";
        }
        oss1 << floor << "." << decimals << ".txt" << std::ends;
        oss2 << floor << "." << decimals << ".txt" << std::ends;
        writeMatrixDLM<ScalarType>(oss1.str().c_str(), forwardPositions[t]);
        writeMatrixDLM<ScalarType>(oss2.str().c_str(), forwardMomentas[t]);
      }
    }
    std::shared_ptr<DeformableMultiObjectType> deformedShape = forwardDef->GetDeformedObjectAt(t);
    deformedShape->WriteMultiObject(filenames);
  }

  // Backward part.
  const ScalarType backwardStepSize = backwardGeodesicLength / (backwardNumberOfTimePoints - 1);
  for (unsigned int t = 1; t < clampedBackwardNumberOfTimePoints; ++t) {
    const ScalarType time = m_ReferenceTime - t * backwardStepSize;
    const std::string roundedTime = std::to_string((unsigned int) (time * 100.0 + 0.5));
    const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
    const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

    std::vector<std::string> filenames(m_NumberOfObjects);
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      std::ostringstream filename;
      filename << outputDir << Superclass::m_Name << "_";
      if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
      filename << "_ReferenceGeodesic__tp_" << backwardNumberOfTimePoints - 1 - t << "__age_";
      if (time < 0) { filename << "-"; }
      filename << floor << "." << decimals << m_TemplateObjectsNameExtension[k];
      filenames[k] = filename.str();
    }
    std::shared_ptr<DeformableMultiObjectType> deformedShape = backwardDef->GetDeformedObjectAt(t);
    deformedShape->WriteMultiObject(filenames);

    if (writeFullTrajectoriesFlag) {
      std::ostringstream oss1, oss2;
      oss1 << outputDir << Superclass::m_Name;
      oss2 << outputDir << Superclass::m_Name;
      oss1 << "__ReferenceGeodesic__CP__tp_" << backwardNumberOfTimePoints - 1 - t << "__age_";
      oss2 << "__ReferenceGeodesic__MOM__tp_" << backwardNumberOfTimePoints - 1 - t << "__age_";
      if (time < 0) {
        oss1 << "-";
        oss2 << "-";
      }
      oss1 << floor << "." << decimals << ".txt" << std::ends;
      oss2 << floor << "." << decimals << ".txt" << std::ends;
      writeMatrixDLM<ScalarType>(oss1.str().c_str(), backwardPositions[t]);
      writeMatrixDLM<ScalarType>(oss2.str().c_str(), backwardMomenta[t]);
    }
  }

  /// Compute and write the model-based reconstruction of the observations.
  const std::string subjectId = timeSeriesDataSet->GetSubjectId();

  // Transport.
  MatrixType spaceShift = (m_ProjectedModulationMatrix * sources).unvectorize(m_NumberOfControlPoints, Dimension);

  // Divide the target absolute times into backward and forward targets.
  std::vector<ScalarType> forwardAbsoluteTimeIncrements, backwardAbsoluteTimeIncrements;
  for (unsigned int t = 0; t < nbObservations; ++t) {
    if (absoluteTimeIncrements[t] < 0.0) {
      backwardAbsoluteTimeIncrements.push_back(-absoluteTimeIncrements[t]);
    } else { forwardAbsoluteTimeIncrements.push_back(absoluteTimeIncrements[t]); }
  }
  std::reverse(backwardAbsoluteTimeIncrements.begin(), backwardAbsoluteTimeIncrements.end());

  // Perform the transport.
  MatrixListType backwardTransportedSpaceShifts = backwardDef->ParallelTransport(
      spaceShift, backwardAbsoluteTimeIncrements).reverse();
  MatrixListType forwardTransportedSpaceShifts = forwardDef->ParallelTransport(
      spaceShift, forwardAbsoluteTimeIncrements);

  // Concatenate the results.
  MatrixListType transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
  assert(transportedSpaceShifts.size() == nbObservations);

  // Exponentiation.
  for (unsigned int t = 0; t < nbObservations; ++t) {
    std::shared_ptr<DeformableMultiObjectType> referenceShape;
    MatrixType referenceControlPoints;

    // Cf. this->GetDeformedObjectAndControlPointsAt, here a const version.
    if (absoluteTimeIncrements[t] >= 0.0) {
      ScalarType continuousIndex = absoluteTimeIncrements[t] / forwardStepSize;
      unsigned int index = continuousIndex;
      if ((continuousIndex - index) >= 0.5) { ++index; }
      referenceShape = forwardDef->GetDeformedObjectAt(index);
      referenceControlPoints = forwardDef->GetDeformedControlPointsAt(index);
    } else {
      ScalarType continuousIndex = -absoluteTimeIncrements[t] / backwardStepSize;
      unsigned int index = continuousIndex;
      if ((continuousIndex - index) >= 0.5) { ++index; }
      referenceShape = backwardDef->GetDeformedObjectAt(index);
      referenceControlPoints = backwardDef->GetDeformedControlPointsAt(index);
    }

    // Riemannian exponential.
    std::shared_ptr<DiffeosType> expDef = m_PerpendicularDeformation->Clone();
    expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
    expDef->SetDeformableMultiObject(referenceShape);
    expDef->SetStartPositions(referenceControlPoints);
    expDef->SetStartMomentas(transportedSpaceShifts[t]);
    expDef->Update();

    std::shared_ptr<DeformableMultiObjectType> predictedShape = expDef->GetDeformedObject();

    // Write.
    const std::string roundedTime = std::to_string((unsigned int) (times[t] * 100.0 + 0.5));
    const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
    const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

    std::vector<std::string> filenames(m_NumberOfObjects);
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      std::ostringstream filename;
      filename << outputDir << Superclass::m_Name << "_";
      if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
      filename << "_Prediction__subject_" << subjectId
               << "__tp_" << t
               << "__age_" << floor << "." << decimals
               << m_TemplateObjectsNameExtension[k];
      filenames[k] = filename.str();
    }
    predictedShape->WriteMultiObject(filenames);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::WriteModelParameters() const {
  const std::string outputDir = def::utils::settings.output_dir;

  // Write time-shift.
  const ScalarType timeShift = recast<ScalarType>(Superclass::m_FixedEffects.at("TimeShift"));
  std::ostringstream oss1;
  oss1 << outputDir << Superclass::m_Name << "__Parameters__TimeShift.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss1.str().c_str(), MatrixType(1, 1, timeShift));

  // Write log-acceleration.
  const ScalarType logAcceleration = recast<ScalarType>(Superclass::m_FixedEffects.at("LogAcceleration"));
  std::ostringstream oss2;
  oss2 << outputDir << Superclass::m_Name << "__Parameters__LogAcceleration.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss2.str().c_str(), MatrixType(1, 1, timeShift));

  // Write sources.
  const VectorType sources = recast<VectorType>(Superclass::m_FixedEffects.at("Sources"));
  std::ostringstream oss3;
  oss3 << outputDir << Superclass::m_Name << "__Parameters__Sources.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss3.str().c_str(), MatrixType(sources));
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::UpdateAbsoluteTimeIncrements(const LongitudinalDataSetType *const dataSet) {
  /// Extract the needed fixed effects.
  const ScalarType timeShift = recast<ScalarType>(Superclass::m_FixedEffects.at("TimeShift"));
  const ScalarType logAcceleration = recast<ScalarType>(Superclass::m_FixedEffects.at("LogAcceleration"));

  /// Initialization.
  const TimeSeriesDataSetType *const timeSeriesDataSet = static_cast<const TimeSeriesDataSetType *const>(dataSet);
  const std::vector<ScalarType> times = timeSeriesDataSet->GetTimes();

  const ScalarType acceleration = std::exp(logAcceleration);
  const unsigned long nbObservations = times.size();

  /// Looping over the observations.
  m_AbsoluteTimeIncrements.resize(nbObservations);
  for (unsigned int t = 0; t < nbObservations; ++t) {
    m_AbsoluteTimeIncrements[t] = acceleration * (times[t] - timeShift);
  }

  m_MinimumAbsoluteTimeIncrement = m_AbsoluteTimeIncrements[0];
  m_MaximumAbsoluteTimeIncrement = m_AbsoluteTimeIncrements[nbObservations - 1];
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::UpdateReferenceGeodesic() {
  /// Shoot to obtain the reference geodesic. Flow to obtain the continuously deformed template objects.
  ScalarType forwardGeodesicLength = m_MaximumAbsoluteTimeIncrement;
  ScalarType backwardGeodesicLength = - m_MinimumAbsoluteTimeIncrement;
  const ScalarType totalGeodesicLength = forwardGeodesicLength + backwardGeodesicLength;

  const ScalarType totalMargin =
      (m_MarginOnGeodesicLength - 1.0) * (m_MaximumAbsoluteTimeIncrement - m_MinimumAbsoluteTimeIncrement);
  const ScalarType forwardMargin = totalMargin * forwardGeodesicLength / totalGeodesicLength;
  const ScalarType backwardMargin = totalMargin * backwardGeodesicLength / totalGeodesicLength;

  // Forward.
  if (forwardGeodesicLength > 0) { forwardGeodesicLength += forwardMargin; }
  const int forwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * forwardGeodesicLength + 1.5;
  m_ForwardReferenceGeodesic->SetT0(0.0);
  m_ForwardReferenceGeodesic->SetTN(forwardGeodesicLength);
  m_ForwardReferenceGeodesic->SetNumberOfTimePoints(forwardNumberOfTimePoints);
  m_ForwardReferenceGeodesic->SetDeformableMultiObject(m_Template);
  m_ForwardReferenceGeodesic->SetStartPositions(m_ControlPoints);
  m_ForwardReferenceGeodesic->SetStartMomentas(m_Momenta);
  if (forwardNumberOfTimePoints > 1) {
    m_ForwardReferenceGeodesic->Update();
    assert(!m_ForwardReferenceGeodesic->OutOfBox());
  }

  // Backward.
  if (backwardGeodesicLength > 0) { backwardGeodesicLength += backwardMargin; }
  const int backwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * backwardGeodesicLength + 1.5;
  m_BackwardReferenceGeodesic->SetT0(0.0);
  m_BackwardReferenceGeodesic->SetTN(backwardGeodesicLength);
  m_BackwardReferenceGeodesic->SetNumberOfTimePoints(backwardNumberOfTimePoints);
  m_BackwardReferenceGeodesic->SetDeformableMultiObject(m_Template);
  m_BackwardReferenceGeodesic->SetStartPositions(m_ControlPoints);
  m_BackwardReferenceGeodesic->SetStartMomentas(-m_Momenta);
  if (backwardNumberOfTimePoints > 1) {
    m_BackwardReferenceGeodesic->Update();
    assert(!m_BackwardReferenceGeodesic->OutOfBox());
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::GetDeformedObjectAndControlPointsAt(ScalarType const &timeIncrement,
                                      std::shared_ptr<DeformableMultiObjectType> &shape,
                                      MatrixType &cp) const {
  /// Get the backward and forward number of time points i.e. the step size.
  const int forwardNumberOfSteps = m_ForwardReferenceGeodesic->GetNumberOfTimePoints();
  const int backwardNumberOfSteps = m_BackwardReferenceGeodesic->GetNumberOfTimePoints();
  const ScalarType forwardStepSize = m_ForwardReferenceGeodesic->GetTN() / (forwardNumberOfSteps - 1);
  const ScalarType backwardStepSize = m_BackwardReferenceGeodesic->GetTN() / (backwardNumberOfSteps - 1);

  /// Find the discrete time that is the closest possible to the target.
  if (timeIncrement >= 0.0) {
    if (forwardNumberOfSteps > 1) {
      ScalarType continuousIndex = timeIncrement / forwardStepSize;
      unsigned int index = continuousIndex;
      if ((continuousIndex - index) >= 0.5) { ++index; }

      shape = m_ForwardReferenceGeodesic->GetDeformedObjectAt(index);
      cp = m_ForwardReferenceGeodesic->GetDeformedControlPointsAt(index);
    } else {
      shape = m_ForwardReferenceGeodesic->GetDeformableMultiObject();
      cp = m_ForwardReferenceGeodesic->GetStartPositions();
    }

  } else {

    if (backwardNumberOfSteps > 1) {
      ScalarType continuousIndex = -timeIncrement / backwardStepSize;
      unsigned int index = continuousIndex;
      if ((continuousIndex - index) >= 0.5) { ++index; }

      shape = m_BackwardReferenceGeodesic->GetDeformedObjectAt(index);
      cp = m_BackwardReferenceGeodesic->GetDeformedControlPointsAt(index);
    } else {
      shape = m_BackwardReferenceGeodesic->GetDeformableMultiObject();
      cp = m_BackwardReferenceGeodesic->GetStartPositions();
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::InitializeBoundingBox() {
  /// Initialization.
  m_BoundingBox = m_Template->GetBoundingBox();
  MatrixType controlPoints = GetControlPoints();

  /// Compute the model bounding box.
  for (unsigned int i = 0; i < controlPoints.rows(); i++) {
    for (unsigned int dim = 0; dim < Dimension; dim++) {
      m_BoundingBox(dim, 0) =
          (m_BoundingBox(dim, 0) < controlPoints(i, dim) ? m_BoundingBox(dim, 0) : controlPoints(i, dim));
      m_BoundingBox(dim, 1) =
          (m_BoundingBox(dim, 1) > controlPoints(i, dim) ? m_BoundingBox(dim, 1) : controlPoints(i, dim));
    }
  }

  /// Set the related bounding boxes.
  m_ForwardReferenceGeodesic->SetDataDomain(m_BoundingBox);
  m_BackwardReferenceGeodesic->SetDataDomain(m_BoundingBox);
  m_PerpendicularDeformation->SetDataDomain(m_BoundingBox);
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetDataDomain(m_BoundingBox);
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalRegistration<ScalarType, Dimension>
::ProjectModulationMatrix() {

  m_ProjectedModulationMatrix = m_ModulationMatrix;
  if (m_Momenta.sum_of_squares() > 1e-20) {
    /// Instantiate the kernel object.
    KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernel = kFactory->CreateKernelObject(m_PerpendicularDeformation->GetKernelType());

    /// Compute the product metric matrix times the momenta vector.
    kernel->SetKernelWidth(m_PerpendicularDeformation->GetKernelWidth());
    kernel->SetSources(m_ControlPoints);
    kernel->SetWeights(m_Momenta);
    const MatrixType Km = kernel->Convolve(m_ControlPoints);

    /// Compute the momenta squared norm.
    const ScalarType mKm = dot_product(m_Momenta, Km);

    /// Vectorized version of the momenta.
    const VectorType KmVectorized = Km.vectorize();
    const VectorType momVectorized = m_Momenta.vectorize();

    /// Loop on the columns.
    for (unsigned int k = 0; k < m_ProjectedModulationMatrix.cols(); ++k) {
      ScalarType cKm = dot_product(m_ProjectedModulationMatrix.get_column(k), KmVectorized);
      m_ProjectedModulationMatrix.set_column(k, m_ProjectedModulationMatrix.get_column(k) - momVectorized * cKm / mKm);
    }
  } else {
    std::cout << "Warning : in LongitudinalRegistration::ProjectModulationMatrix, the given momenta seem to be null."
              << std::endl;
  }
}

template
class LongitudinalRegistration<double, 2>;
template
class LongitudinalRegistration<double, 3>;
