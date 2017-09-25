/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "LongitudinalAtlas.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
LongitudinalAtlas<ScalarType, Dimension>
::LongitudinalAtlas() : Superclass(),
                        m_Template(NULL),
                        m_Def(NULL),
                        m_ForwardReferenceGeodesic(NULL),
                        m_BackwardReferenceGeodesic(NULL) {
  Superclass::SetLongitudinalAtlasType();

  /// Fixed effects.
  ScalarType referenceTimeFE(sqrt(-1.0)), logAccelerationVarianceFE(-1.0), timeShiftVarianceFE(-1.0);
  VectorType noiseVarianceFE(1, -1.0);
  MatrixType controlPointsFE, momentaFE, modulationMatrixFE;
  MatrixListType templateDataFE;
  Superclass::m_FixedEffects["TemplateData"] = templateDataFE;
  Superclass::m_FixedEffects["ControlPoints"] = controlPointsFE;
  Superclass::m_FixedEffects["Momenta"] = momentaFE;
  Superclass::m_FixedEffects["ModulationMatrix"] = modulationMatrixFE;
  Superclass::m_FixedEffects["ReferenceTime"] = referenceTimeFE;
  Superclass::m_FixedEffects["TimeShiftVariance"] = timeShiftVarianceFE;
  Superclass::m_FixedEffects["LogAccelerationVariance"] = logAccelerationVarianceFE;
  Superclass::m_FixedEffects["NoiseVariance"] = noiseVarianceFE;

  /// Random effects.
  auto templateDataRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_PopulationRandomEffects["TemplateData"] = templateDataRE;
  auto controlPointsRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_PopulationRandomEffects["ControlPoints"] = controlPointsRE;
  auto momentaRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_PopulationRandomEffects["Momenta"] = momentaRE;
  auto modulationMatrixRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_PopulationRandomEffects["ModulationMatrix"] = modulationMatrixRE;

  auto sourcesRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_IndividualRandomEffects["Sources"] = sourcesRE;
  auto timeShiftRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_IndividualRandomEffects["TimeShift"] = timeShiftRE;
  auto logAccelerationRE = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_IndividualRandomEffects["LogAcceleration"] = logAccelerationRE;

  /// Bayesian priors.
  auto templateDataPrior = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_Priors["TemplateData"] = templateDataPrior;
  auto controlPointsPrior = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_Priors["ControlPoints"] = controlPointsPrior;
  auto momentaPrior = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_Priors["Momenta"] = momentaPrior;
  auto modulationMatrixPrior = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_Priors["ModulationMatrix"] = modulationMatrixPrior;
  auto referenceTimePrior = std::make_shared<MultiScalarNormalDistributionType>();
  this->m_Priors["ReferenceTime"] = referenceTimePrior;

  auto timeShiftVariancePrior = std::make_shared<MultiScalarInverseWishartDistributionType>();
  this->m_Priors["TimeShiftVariance"] = timeShiftVariancePrior;
  auto logAccelerationVariancePrior = std::make_shared<MultiScalarInverseWishartDistributionType>();
  this->m_Priors["LogAccelerationVariance"] = logAccelerationVariancePrior;
  auto noiseVariancePrior = std::make_shared<MultiScalarInverseWishartDistributionType>();
  this->m_Priors["NoiseVariance"] = noiseVariancePrior;
}

template<class ScalarType, unsigned int Dimension>
LongitudinalAtlas<ScalarType, Dimension>
::~LongitudinalAtlas() {}

template<class ScalarType, unsigned int Dimension>
LongitudinalAtlas<ScalarType, Dimension>
::LongitudinalAtlas(const LongitudinalAtlas &other) : Superclass(other) {
  m_Template = other.m_Template->Clone();
  m_TemplateObjectsName = other.m_TemplateObjectsName;
  m_TemplateObjectsNameExtension = other.m_TemplateObjectsNameExtension;
  m_TemplateDataSizeParameters = other.m_TemplateDataSizeParameters;
  m_BoundingBox = other.m_BoundingBox;

  m_Def = other.m_Def->Clone();
  m_ConcentrationOfTimePointsForReferenceGeodesic = other.m_ConcentrationOfTimePointsForReferenceGeodesic;
  m_NumberOfTimePointsForExponentiation = other.m_NumberOfTimePointsForExponentiation;
  m_MarginOnGeodesicLength = other.m_MarginOnGeodesicLength;

  m_NumberOfObjects = other.m_NumberOfObjects;
  m_DimensionOfDiscretizedObjects = other.m_DimensionOfDiscretizedObjects;
  m_NumberOfControlPoints = other.m_NumberOfControlPoints;
  m_DimensionOfTangentSpaces = other.m_DimensionOfTangentSpaces;
  m_NumberOfSources = other.m_NumberOfSources;

  m_CPSpacing = other.m_CPSpacing;
  m_FreezeControlPointsFlag = other.m_FreezeControlPointsFlag;
  m_FrozenControlPoints = other.m_FrozenControlPoints;

  m_AbsoluteTimeIncrements = other.m_AbsoluteTimeIncrements;
  m_AbsoluteTimeIncrements_Memory = other.m_AbsoluteTimeIncrements_Memory;
  m_MinimumAbsoluteTimeIncrement = other.m_MinimumAbsoluteTimeIncrement;
  m_MinimumAbsoluteTimeIncrement_Memory = other.m_MinimumAbsoluteTimeIncrement_Memory;
  m_MaximumAbsoluteTimeIncrement = other.m_MaximumAbsoluteTimeIncrement;
  m_MaximumAbsoluteTimeIncrement_Memory = other.m_MaximumAbsoluteTimeIncrement_Memory;

  m_ForwardReferenceGeodesic = other.m_ForwardReferenceGeodesic->Clone();
  m_ForwardReferenceGeodesic_Memory = other.m_ForwardReferenceGeodesic_Memory->Clone();
  m_BackwardReferenceGeodesic = other.m_BackwardReferenceGeodesic->Clone();
  m_BackwardReferenceGeodesic_Memory = other.m_BackwardReferenceGeodesic_Memory->Clone();

  m_ProjectedModulationMatrix = other.m_ProjectedModulationMatrix;
  m_ProjectedModulationMatrix_Memory = other.m_ProjectedModulationMatrix_Memory;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::SetFixedEffects(const LinearVariableMapType &map) {
  Superclass::SetFixedEffects(map);

  SetTemplateData(recast<MatrixListType>(map.at("TemplateData")));
  if (!m_FreezeControlPointsFlag) { SetControlPoints(recast<MatrixType>(map.at("ControlPoints"))); }
  SetMomenta(recast<MatrixType>(map.at("Momenta")));
  SetModulationMatrix(recast<MatrixType>(map.at("ModulationMatrix")));
  SetReferenceTime(recast<ScalarType>(map.at("ReferenceTime")));

  SetTimeShiftVariance(recast<ScalarType>(map.at("TimeShiftVariance")));
  SetLogAccelerationVariance(recast<ScalarType>(map.at("LogAccelerationVariance")));
  SetNoiseVariance(recast<VectorType>(map.at("NoiseVariance")));
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::SetControlPoints(const std::string &fn) {
  if (strlen(fn.c_str())) {
    const MatrixType controlPoints = readMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using a set of " << controlPoints.rows() << " control points in file " << fn << std::endl;

    SetControlPoints(controlPoints);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::SetMomenta(const std::string &fn) {
  if (strlen(fn.c_str())) {
    const MatrixType momenta = readMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using a momenta matrix of size " << momenta.rows() << " x "
              << momenta.columns() << " from file: " << fn << std::endl;

    SetMomenta(momenta);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::Update() {
  if (m_Def == NULL)
    throw std::runtime_error("A deformation should be set to the LongitudinalAtlas model.");

  InitializeTemplateDataVariables();
  InitializeBoundingBox();
  InitializeControlPointsVariables();
  InitializeMomentaVariables();
  InitializeModulationMatrixVariables();
  InitializeReferenceTimeVariables();

  InitializeSources();
  InitializeTimeShiftVariables();
  InitializeLogAccelerationVariables();
  InitializeNoiseVariables();
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::Print() const {
  const ScalarType referenceTime = Superclass::m_FixedEffects["ReferenceTime"].vectorize()[0];
  const ScalarType timeShiftVariance = Superclass::m_FixedEffects["TimeShiftVariance"].vectorize()[0];
  const ScalarType logAccelerationVariance = Superclass::m_FixedEffects["LogAccelerationVariance"].vectorize()[0];
  const VectorType noiseVariance = Superclass::m_FixedEffects["NoiseVariance"].vectorize();

  std::cout << ">> Model fixed effects :" << std::endl;
  std::cout << "\t\tReference time = " << referenceTime << std::endl;
  std::cout << "\t\tTime-shift standard deviation = " << std::sqrt(timeShiftVariance) << std::endl;
  std::cout << "\t\tLog-acceleration standard deviation = " << std::sqrt(logAccelerationVariance) << std::endl;
  for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
    std::cout << "\t\tNoise standard deviation ";
    if (m_NumberOfObjects > 1) { std::cout << "(" << m_TemplateObjectsName[k] << ")"; }
    std::cout << " = " << std::sqrt(noiseVariance[k]) << std::endl;
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::Write(const LongitudinalDataSetType *const dataSet,
        LinearVariableMapType const &popRER,
        LinearVariablesMapType const &indRER) const {
  WriteModelPredictions(dataSet, popRER, indRER);
  WriteModelParameters(indRER);
}

template<class ScalarType, unsigned int Dimension>
bool
LongitudinalAtlas<ScalarType, Dimension>
::Sample(LongitudinalDataSetType *dataSet,
         LinearVariableMapType &popRER,
         LinearVariablesMapType &indRER) const {
  /// Initialization.
  const std::vector<std::vector<ScalarType>> times = dataSet->GetTimes();
  const unsigned int numberOfSubjects = times.size();

  /// Get the needed population fixed effects.
  const MatrixListType templateData = GetTemplateData();
  const std::shared_ptr<DeformableMultiObjectType> temp = m_Template->Clone();
  temp->UpdateImageIntensityAndLandmarkPointCoordinates(templateData);
  temp->Update();

  const MatrixType controlPoints = GetControlPoints();
  const MatrixType momenta = GetMomenta();
  const MatrixType modulationMatrix = GetModulationMatrix();

  /// Simulate the needed individual random effects.
  std::vector<VectorType> sourcesRERs(numberOfSubjects);
  std::vector<ScalarType> logAccelerationRERs(numberOfSubjects);
  std::vector<ScalarType> timeShiftRERs(numberOfSubjects);

  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    sourcesRERs[i] = std::static_pointer_cast<MultiScalarNormalDistributionType>(
        Superclass::m_IndividualRandomEffects.at("Sources"))->Sample();
    timeShiftRERs[i] = std::static_pointer_cast<MultiScalarNormalDistributionType>(
        Superclass::m_IndividualRandomEffects.at("TimeShift"))->Sample()(0);
    logAccelerationRERs[i] = std::static_pointer_cast<MultiScalarNormalDistributionType>(
        Superclass::m_IndividualRandomEffects.at("LogAcceleration"))->Sample()(0);
  }

  /// Compute the absolute times (cf. this->UpdateAbsoluteTimeIncrements, here a const version).
  ScalarType minimumAbsoluteTimeIncrement = 0.0;
  ScalarType maximumAbsoluteTimeIncrement = 0.0;

  std::vector<std::vector<ScalarType>> absoluteTimeIncrements(numberOfSubjects);
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    const unsigned long nbObservations_i = times[i].size();
    const ScalarType acceleration_i = std::exp(logAccelerationRERs[i]);

    absoluteTimeIncrements[i].resize(nbObservations_i);
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      absoluteTimeIncrements[i][t] = acceleration_i * (times[i][t] - timeShiftRERs[i]);
    }

    if (absoluteTimeIncrements[i][0] < minimumAbsoluteTimeIncrement) {
      minimumAbsoluteTimeIncrement = absoluteTimeIncrements[i][0];
    }
    if (absoluteTimeIncrements[i][nbObservations_i - 1] > maximumAbsoluteTimeIncrement) {
      maximumAbsoluteTimeIncrement = absoluteTimeIncrements[i][nbObservations_i - 1];
    }
  }

  /// Compute the reference geodesic (cf. this->UpdateReferenceGeodesic, here a const version).
  // Forward part.
  std::shared_ptr<DiffeosType> forwardDef = m_ForwardReferenceGeodesic->Clone();
  const ScalarType forwardGeodesicLength = maximumAbsoluteTimeIncrement;
  const int forwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * forwardGeodesicLength + 1.5;
  forwardDef->SetT0(0.0);
  forwardDef->SetTN(forwardGeodesicLength);
  forwardDef->SetNumberOfTimePoints(forwardNumberOfTimePoints);
  forwardDef->SetDeformableMultiObject(temp);
  forwardDef->SetStartPositions(controlPoints);
  forwardDef->SetStartMomentas(momenta);
  forwardDef->Update();

  // Backward part.
  std::shared_ptr<DiffeosType> backwardDef = m_BackwardReferenceGeodesic->Clone();
  const ScalarType backwardGeodesicLength = -minimumAbsoluteTimeIncrement;
  const int backwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * backwardGeodesicLength + 1.5;
  backwardDef->SetT0(0.0);
  backwardDef->SetTN(backwardGeodesicLength);
  backwardDef->SetNumberOfTimePoints(backwardNumberOfTimePoints);
  backwardDef->SetDeformableMultiObject(temp);
  backwardDef->SetStartPositions(controlPoints);
  backwardDef->SetStartMomentas(-momenta);
  backwardDef->Update();

  /// Compute the projected modulation matrix (cf. this->UpdateProjectedModulationMatrix, here a const version).
  MatrixType projectedModulationMatrixRER = modulationMatrix;
  if (momenta.sum_of_squares() > 1e-20) {
    // Instantiate the kernel object.
    KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernel = kFactory->CreateKernelObject(m_Def->GetKernelType());

    // Compute the product metric matrix times the momenta vector.
    kernel->SetKernelWidth(m_Def->GetKernelWidth());
    kernel->SetSources(controlPoints);
    kernel->SetWeights(momenta);
    const MatrixType Km = kernel->Convolve(controlPoints);

    // Compute the momenta squared norm.
    const ScalarType mKm = dot_product(momenta, Km);

    // Vectorized version of the momenta.
    const VectorType KmVectorized = Km.vectorize();
    const VectorType momVectorized = momenta.vectorize();

    // Loop on the columns.
    for (unsigned int k = 0; k < projectedModulationMatrixRER.cols(); ++k) {
      ScalarType cKm = dot_product(projectedModulationMatrixRER.get_column(k), KmVectorized);
      projectedModulationMatrixRER.set_column(
          k, projectedModulationMatrixRER.get_column(k) - momVectorized * cKm / mKm);
    }
  }

  /// For each subject, transport along the reference geodesic and then shoot at the target time-points.
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType >>> samples(numberOfSubjects);
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    const unsigned int nbObservations_i = absoluteTimeIncrements[i].size();

    /// Transport.
    MatrixType spaceShift = (projectedModulationMatrixRER * sourcesRERs[i])
        .unvectorize(m_NumberOfControlPoints, Dimension);

    // Divide the target absolute times into backward and forward targets.
    std::vector<ScalarType> forwardabsoluteTimeIncrements, backwardabsoluteTimeIncrements;
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      if (absoluteTimeIncrements[i][t] < 0.0) {
        backwardabsoluteTimeIncrements.push_back(-absoluteTimeIncrements[i][t]);
      } else { forwardabsoluteTimeIncrements.push_back(absoluteTimeIncrements[i][t]); }
    }
    std::reverse(backwardabsoluteTimeIncrements.begin(), backwardabsoluteTimeIncrements.end());

    // Perform the transport.
    MatrixListType backwardTransportedSpaceShifts = backwardDef->ParallelTransport(
        spaceShift, backwardabsoluteTimeIncrements).reverse();
    MatrixListType forwardTransportedSpaceShifts = forwardDef->ParallelTransport(
        spaceShift, forwardabsoluteTimeIncrements);

    // Concatenate the results.
    MatrixListType transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
    assert(transportedSpaceShifts.size() == nbObservations_i);

    /// Exponentiation.
    const ScalarType forwardStepSize = forwardGeodesicLength / (forwardNumberOfTimePoints - 1);
    const ScalarType backwardStepSize = backwardGeodesicLength / (backwardNumberOfTimePoints - 1);

    samples[i].resize(nbObservations_i);
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      std::shared_ptr<DeformableMultiObjectType> referenceShape;
      MatrixType referenceControlPoints;

      // Cf. this->GetDeformedObjectAndControlPointsAt, here a const version.
      if (absoluteTimeIncrements[i][t] >= 0.0) {
        ScalarType continuousIndex = absoluteTimeIncrements[i][t] / forwardStepSize;
        unsigned int index = continuousIndex;
        if ((continuousIndex - index) >= 0.5) { ++index; }
        referenceShape = forwardDef->GetDeformedObjectAt(index);
        referenceControlPoints = forwardDef->GetDeformedControlPointsAt(index);
      } else {
        ScalarType continuousIndex = -absoluteTimeIncrements[i][t] / backwardStepSize;
        unsigned int index = continuousIndex;
        if ((continuousIndex - index) >= 0.5) { ++index; }
        referenceShape = backwardDef->GetDeformedObjectAt(index);
        referenceControlPoints = backwardDef->GetDeformedControlPointsAt(index);
      }

      std::shared_ptr<DiffeosType> expDef = m_Def->Clone();
      expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
      expDef->SetDeformableMultiObject(referenceShape);
      expDef->SetStartPositions(referenceControlPoints);
      expDef->SetStartMomentas(transportedSpaceShifts[t]);
      expDef->Update();

      samples[i][t] = expDef->GetDeformedObject();
    }
  }

  /// Exports the results.
  // Population parameters.
  popRER["TemplateData"] = templateData;
  if (!m_FreezeControlPointsFlag) { popRER["ControlPoints"] = controlPoints; }
  popRER["Momenta"] = momenta;
  popRER["ModulationMatrix"] = modulationMatrix;

  // Individual parameters.
  indRER["Sources"] = sourcesRERs;
  indRER["LogAcceleration"] = logAccelerationRERs;
  indRER["TimeShift"] = timeShiftRERs;

  // Generated dataset.
  dataSet->SetDeformableMultiObjects(samples);

  return false; // TODO : return true if an out-of-box is detected.
}

template<class ScalarType, unsigned int Dimension>
bool
LongitudinalAtlas<ScalarType, Dimension>
::ComputeResiduals(const LongitudinalDataSetType *const dataSet,
                   const LinearVariableMapType &popRER,
                   const LinearVariablesMapType &indRER,
                   std::vector<std::vector<std::vector<ScalarType>>> &residuals) {
  /// Extract needed variables.
  const unsigned int numberOfSubjects = dataSet->GetNumberOfSubjects();
  const std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>>
      targets = dataSet->GetDeformableMultiObjects();
  const std::vector<VectorType> sourcesRERs = recast<VectorType>(indRER.at("Sources"));

  /// For each subject, transport along the reference geodesic and then shoot at the target time-points.
  residuals.resize(numberOfSubjects);
  {
    ThreadPool pool(def::utils::settings.number_of_threads);
    for (unsigned int i = 0; i < numberOfSubjects; ++i) {
      pool.enqueue([&, i]() {
        const unsigned int nbObservations_i = m_AbsoluteTimeIncrements[i].size();

        /// Transport.
        MatrixType spaceShift = (m_ProjectedModulationMatrix * sourcesRERs[i])
            .unvectorize(m_NumberOfControlPoints, Dimension);

        // Divide the target absolute times into backward and forward targets.
        std::vector<ScalarType> forwardabsoluteTimeIncrements, backwardabsoluteTimeIncrements;
        for (unsigned int t = 0; t < nbObservations_i; ++t) {
          if (m_AbsoluteTimeIncrements[i][t] < 0.0) {
            backwardabsoluteTimeIncrements.push_back(-m_AbsoluteTimeIncrements[i][t]);
          } else { forwardabsoluteTimeIncrements.push_back(m_AbsoluteTimeIncrements[i][t]); }
        }
        std::reverse(backwardabsoluteTimeIncrements.begin(), backwardabsoluteTimeIncrements.end());

        // Perform the transport.
        MatrixListType backwardTransportedSpaceShifts = m_BackwardReferenceGeodesic->ParallelTransport(
            spaceShift, backwardabsoluteTimeIncrements).reverse();
        MatrixListType forwardTransportedSpaceShifts = m_ForwardReferenceGeodesic->ParallelTransport(
            spaceShift, forwardabsoluteTimeIncrements);

        // Concatenate the results.
        MatrixListType
            transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
        assert(transportedSpaceShifts.size() == nbObservations_i);

        /// Exponentiation.
        residuals[i].resize(nbObservations_i);
        for (unsigned int t = 0; t < nbObservations_i; ++t) {
          std::shared_ptr<DeformableMultiObjectType> referenceShape;
          MatrixType referenceControlPoints;
          GetDeformedObjectAndControlPointsAt(m_AbsoluteTimeIncrements[i][t], referenceShape, referenceControlPoints);

          std::shared_ptr<DiffeosType> expDef = m_Def->Clone();
          expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
          expDef->SetDeformableMultiObject(referenceShape);
          expDef->SetStartPositions(referenceControlPoints);
          expDef->SetStartMomentas(transportedSpaceShifts[t]);
          expDef->Update();

          std::shared_ptr<DeformableMultiObjectType> predictedShape = expDef->GetDeformedObject();
          residuals[i][t] = predictedShape->ComputeMatch(targets[i][t]);
        }
      });
    }
  }
  return false; // TODO : return true if an out-of-box is detected.
}

template<class ScalarType, unsigned int Dimension>
ScalarType
LongitudinalAtlas<ScalarType, Dimension>
::ComputeCompleteLogLikelihood(const LongitudinalDataSetType *const dataSet,
                               const LinearVariableMapType &popRER,
                               const LinearVariablesMapType &indRER) {
  /// Initialization.
  UpdateAbsoluteTimeIncrements(dataSet, popRER, indRER);
  UpdateReferenceGeodesic(popRER);
  UpdateProjectedModulationMatrix(popRER);

  /// Data (residuals) term.
  std::vector<std::vector<std::vector<ScalarType>>> residuals;
  ComputeResiduals(dataSet, popRER, indRER, residuals);

  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const unsigned int numberOfSubjects = residuals.size();

  ScalarType out = 0.0;
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    for (unsigned int t = 0; t < residuals[i].size(); ++t) {
      for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
        out -= 0.5 * (residuals[i][t][k] / noiseVariances(k)
            + m_DimensionOfDiscretizedObjects[k] * log(noiseVariances(k)));
      }
    }
  }

  /// Population random effects terms.
  // Template term.
  out += this->m_PopulationRandomEffects["TemplateData"]->ComputeLogLikelihood(popRER.at("TemplateData"));

  // Control points term.
  if (!m_FreezeControlPointsFlag) {
    out += this->m_PopulationRandomEffects["ControlPoints"]->ComputeLogLikelihood(popRER.at("ControlPoints"));
  }

  // Momenta term.
  out += this->m_PopulationRandomEffects["Momenta"]->ComputeLogLikelihood(popRER.at("Momenta"));

  // Modulation matrix term.
  out += this->m_PopulationRandomEffects["ModulationMatrix"]->ComputeLogLikelihood(popRER.at("ModulationMatrix"));

  /// Individual random effects terms.
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    // Sources term.
    out += this->m_IndividualRandomEffects["Sources"]->ComputeLogLikelihood(indRER.at("Sources")[i]);

    // Time shift term.
    out += this->m_IndividualRandomEffects["TimeShift"]->ComputeLogLikelihood(indRER.at("TimeShift")[i]);

    // Log-acceleration term.
    out += this->m_IndividualRandomEffects["LogAcceleration"]->ComputeLogLikelihood(indRER.at("LogAcceleration")[i]);
  }

  /// Prior terms.
  // Template term.
  out += this->m_Priors["TemplateData"]->ComputeLogLikelihood(this->m_FixedEffects["TemplateData"]);

  // Control points term.
  if (!m_FreezeControlPointsFlag) {
    out += this->m_Priors["ControlPoints"]->ComputeLogLikelihood(this->m_FixedEffects["ControlPoints"]);
  }

  // Momenta term.
  out += this->m_Priors["Momenta"]->ComputeLogLikelihood(this->m_FixedEffects["Momenta"]);

  // Modulation matrix term.
  out += this->m_Priors["ModulationMatrix"]->ComputeLogLikelihood(this->m_FixedEffects["ModulationMatrix"]);

  // Reference time term.
  out += this->m_Priors["ReferenceTime"]->ComputeLogLikelihood(this->m_FixedEffects["ReferenceTime"]);

  // Time-shift variance term.
  out += this->m_Priors["TimeShiftVariance"]->ComputeLogLikelihood(this->m_FixedEffects["TimeShiftVariance"]);

  // Log-acceleration variance term.
  out +=
      this->m_Priors["LogAccelerationVariance"]->ComputeLogLikelihood(this->m_FixedEffects["LogAccelerationVariance"]);

  // Noise variance term.
  out += this->m_Priors["NoiseVariance"]->ComputeLogLikelihood(this->m_FixedEffects["NoiseVariance"]);

  return out;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
LongitudinalAtlas<ScalarType, Dimension>
::ComputeModelLogLikelihood(const LongitudinalDataSetType *const dataSet,
                            const LinearVariableMapType &popRER,
                            const LinearVariablesMapType &indRER,
                            const ScalarType &temperature,
                            const std::string &modifiedVar,
                            std::vector<ScalarType> &contributions) {
  /// Update the relevant memorized intermediate variables based on the random effects realizations.
  if (modifiedVar == "TemplateData") { // NOTE : ONLY FLOW() IS ACTUALLY NECESSARY. TODO.
    UpdateReferenceGeodesic(popRER);
  } else if (!m_FreezeControlPointsFlag and modifiedVar == "ControlPoints") {
    UpdateReferenceGeodesic(popRER);
    UpdateProjectedModulationMatrix(popRER);
  } else if (modifiedVar == "Momenta") {
    UpdateReferenceGeodesic(popRER);
    UpdateProjectedModulationMatrix(popRER);
  } else if (modifiedVar == "ModulationMatrix") {
    UpdateProjectedModulationMatrix(popRER);
  } else if (modifiedVar == "Sources") {
  } else if (modifiedVar == "TimeShift" || modifiedVar == "LogAcceleration") {
    UpdateAbsoluteTimeIncrements(dataSet, popRER, indRER);

    // Memorize the current state. Necessary if the geodesic has been updated at the previous call.
    *m_ForwardReferenceGeodesic_Memory = *m_ForwardReferenceGeodesic;
    *m_BackwardReferenceGeodesic_Memory = *m_BackwardReferenceGeodesic;

    if ((m_MinimumAbsoluteTimeIncrement < 0.0 &&
        m_MinimumAbsoluteTimeIncrement < (-m_BackwardReferenceGeodesic->GetTN())) ||
        (m_MaximumAbsoluteTimeIncrement > 0.0 &&
            m_MaximumAbsoluteTimeIncrement > m_ForwardReferenceGeodesic->GetTN())) {
      UpdateReferenceGeodesic(popRER); // SHOULD ONLY BE EXTENDED. TODO.
    }
  } else {
    UpdateAbsoluteTimeIncrements(dataSet, popRER, indRER);
    UpdateReferenceGeodesic(popRER);
    UpdateProjectedModulationMatrix(popRER);
  }

  /// Data (residuals) term.
  std::vector<std::vector<std::vector<ScalarType>>> residuals;
  ComputeResiduals(dataSet, popRER, indRER, residuals);

  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);
  const unsigned int numberOfSubjects = residuals.size();

  contributions.resize(numberOfSubjects);
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    contributions[i] = 0.0;
    for (unsigned int t = 0; t < residuals[i].size(); ++t) {
      for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
        contributions[i] -= 0.5 * (residuals[i][t][k] / (temperature * noiseVariances(k)));
      }
    }
  }

  ScalarType out = std::accumulate(contributions.begin(), contributions.end(), 0.0);
  return out;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
LongitudinalAtlas<ScalarType, Dimension>
::ComputeModelLogLikelihoodForSubject(const LongitudinalDataSetType *const dataSet,
                                      const LinearVariableMapType &popRER,
                                      const LinearVariablesMapType &indRER,
                                      const unsigned int &i,
                                      const std::string &modifiedVar) {
  /// Update the relevant memorized intermediate variables based on the random effects realizations.
  if (modifiedVar != "Sources") { // I.e. LogAcceleration, TimeShift or default.
    UpdateAbsoluteTimeIncrements(dataSet, popRER, indRER);

    // Memorize the current state. Necessary if the geodesic has been updated at the previous call.
    *m_ForwardReferenceGeodesic_Memory = *m_ForwardReferenceGeodesic;
    *m_BackwardReferenceGeodesic_Memory = *m_BackwardReferenceGeodesic;

    if ((m_MinimumAbsoluteTimeIncrement < 0.0 &&
        m_MinimumAbsoluteTimeIncrement < (-m_BackwardReferenceGeodesic->GetTN())) ||
        (m_MaximumAbsoluteTimeIncrement > 0.0 &&
            m_MaximumAbsoluteTimeIncrement > m_ForwardReferenceGeodesic->GetTN())) {
//      std::cout << ">> Warning : the margin on the geodesic length seems unsufficient." << std::endl;
      UpdateReferenceGeodesic(popRER);
    }
  }

  /// Data (residuals) term for a specific subject.
  std::vector<std::vector<ScalarType>> residuals = ComputeResidualsForSubject(dataSet, popRER, indRER, i);
  const VectorType noiseVariances = recast<VectorType>(this->m_FixedEffects["NoiseVariance"]);

  ScalarType out = 0.0;
  for (unsigned int t = 0; t < residuals.size(); ++t) {
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      out -= 0.5 * (residuals[t][k] / noiseVariances(k));
    }
  }

  return out;
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::ComputeSufficientStatistics(const LongitudinalDataSetType *const dataSet,
                              const LinearVariableMapType &popRER,
                              const LinearVariablesMapType &indRER,
                              LinearVariableMapType &sufficientStatistics) {
  /// Initialization.
  sufficientStatistics.clear();
  UpdateAbsoluteTimeIncrements(dataSet, popRER, indRER);
  UpdateReferenceGeodesic(popRER);
  UpdateProjectedModulationMatrix(popRER);

  std::vector<std::vector<std::vector<ScalarType>>> residuals;
  ComputeResiduals(dataSet, popRER, indRER, residuals);

  /// Population sufficient statistics.
  sufficientStatistics["S1"] = recast<MatrixListType>(popRER.at("TemplateData"));
  if (!m_FreezeControlPointsFlag) { sufficientStatistics["S2"] = recast<MatrixType>(popRER.at("ControlPoints")); }
  sufficientStatistics["S3"] = recast<MatrixType>(popRER.at("Momenta"));
  sufficientStatistics["S4"] = recast<MatrixType>(popRER.at("ModulationMatrix"));

  /// Individual sufficient statistics.
  sufficientStatistics["S5"] = indRER.at("TimeShift").sum();
  sufficientStatistics["S6"] = indRER.at("TimeShift").sum_of_squares();
  sufficientStatistics["S7"] = indRER.at("LogAcceleration").sum_of_squares();

  // Special case of the residuals. Initialization.
  unsigned int numberOfSubjects = dataSet->GetNumberOfSubjects();
  sufficientStatistics["S8"] = VectorType(residuals[0][0]);
  for (unsigned int j = 1; j < residuals[0].size(); ++j) {
    sufficientStatistics["S8"] += VectorType(residuals[0][j]);
  }

  // Summing over all subjects.
  for (unsigned int i = 1; i < numberOfSubjects; ++i) {
    for (unsigned int j = 0; j < residuals[i].size(); ++j) {
      sufficientStatistics["S8"] += VectorType(residuals[i][j]);
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::UpdateFixedEffects(const LongitudinalDataSetType *const dataSet,
                     const LinearVariableMapType &sufficientStatistics,
                     const ScalarType &temperature) {

  const unsigned int numberOfSubjects = dataSet->GetNumberOfSubjects();
  const MatrixListType S1 = recast<MatrixListType>(sufficientStatistics["S1"]);
  const MatrixType S3 = recast<MatrixType>(sufficientStatistics["S3"]);
  const MatrixType S4 = recast<MatrixType>(sufficientStatistics["S4"]);
  const ScalarType S5 = recast<ScalarType>(sufficientStatistics["S5"]);
  const ScalarType S6 = recast<ScalarType>(sufficientStatistics["S6"]);
  const ScalarType S7 = recast<ScalarType>(sufficientStatistics["S7"]);
  const VectorType S8 = recast<VectorType>(sufficientStatistics["S8"]);

  /// Template data update.
  const ScalarType templateDataRandomEffectVariance = GetTemplateDataRandomEffectVariance();
  const MatrixListType templateDataPriorMean = GetTemplateDataPriorMean();
  const ScalarType templateDataPriorVariance = GetTemplateDataPriorVariance();
  const MatrixListType newTemplateData =
      (templateDataPriorVariance * S1 + temperature * templateDataRandomEffectVariance * templateDataPriorMean)
          / (templateDataPriorVariance + temperature * templateDataRandomEffectVariance);
  SetTemplateData(newTemplateData);

  /// Control points update.
  if (!m_FreezeControlPointsFlag) {
    const MatrixType S2 = recast<MatrixType>(sufficientStatistics["S2"]);
    const ScalarType controlPointsRandomEffectVariance = GetControlPointsRandomEffectVariance();
    const MatrixType controlPointsPriorMean = GetControlPointsPriorMean();
    const ScalarType controlPointsPriorVariance = GetControlPointsPriorVariance();
    const MatrixType newControlPoints =
        (controlPointsPriorVariance * S2 + temperature * controlPointsRandomEffectVariance * controlPointsPriorMean)
            / (controlPointsPriorVariance + temperature * controlPointsRandomEffectVariance);
    SetControlPoints(newControlPoints);
  }

  /// Momenta update.
  const ScalarType momentaRandomEffectVariance = GetMomentaRandomEffectVariance();
  const MatrixType momentaPriorMean = GetMomentaPriorMean();
  const ScalarType momentaPriorVariance = GetMomentaPriorVariance();
  const MatrixType newMomenta =
      (momentaPriorVariance * S3 + temperature * momentaRandomEffectVariance * momentaPriorMean)
          / (momentaPriorVariance + temperature * momentaRandomEffectVariance);
  SetMomenta(newMomenta);

  /// Modulation matrix update.
  const ScalarType modulationMatrixRandomEffectVariance = GetModulationMatrixRandomEffectVariance();
  const MatrixType modulationMatrixPriorMean = GetModulationMatrixPriorMean();
  const ScalarType modulationMatrixPriorVariance = GetModulationMatrixPriorVariance();
  const MatrixType newModulationMatrix =
      (modulationMatrixPriorVariance * S4
          + temperature * modulationMatrixRandomEffectVariance * modulationMatrixPriorMean)
          / (modulationMatrixPriorVariance + temperature * modulationMatrixRandomEffectVariance);
  SetModulationMatrix(newModulationMatrix);

  /// Intricated update of the reference time and the time-shift variance.
  const ScalarType referenceTimePriorMean = GetReferenceTimePriorMean();
  const ScalarType referenceTimePriorVariance = GetReferenceTimePriorVariance();
  const ScalarType timeShiftVariancePriorDof = GetTimeShiftVariancePriorDegreeOfFreedom();
  const ScalarType timeShiftVariancePriorScale = GetTimeShiftVariancePriorScaleFactor();

  ScalarType referenceTime_old = GetReferenceTime(), referenceTime_new = referenceTime_old;
  ScalarType timeShiftVariance_old = GetTimeShiftVariance(), timeShiftVariance_new = timeShiftVariance_old;

  const unsigned int maxNbIterations = 20;
  const ScalarType tolerance = 1e-5;
  unsigned int k = 0;
  for (; k < maxNbIterations; ++k) {
    referenceTime_new = (referenceTimePriorVariance * S5 + temperature * timeShiftVariance_new * referenceTimePriorMean)
        / (numberOfSubjects * referenceTimePriorVariance + temperature * timeShiftVariance_new);
    timeShiftVariance_new = (S6 - 2 * referenceTime_new * S5 + numberOfSubjects * pow(referenceTime_new, 2)
        + timeShiftVariancePriorDof * timeShiftVariancePriorScale)
        / (temperature * (numberOfSubjects + timeShiftVariancePriorDof));

    if (std::fabs(referenceTime_new - referenceTime_old) < tolerance &&
        std::fabs(timeShiftVariance_new - timeShiftVariance_old) < tolerance) { break; }
    else {
      referenceTime_old = referenceTime_new;
      timeShiftVariance_old = timeShiftVariance_new;
    }
  }
  assert(k < maxNbIterations - 1);

  SetReferenceTime(referenceTime_new);
  SetTimeShiftVariance(timeShiftVariance_new);

  /// Log-acceleration variance update.
  const ScalarType logAccelerationVariancePriorDof = GetLogAccelerationVariancePriorDegreeOfFreedom();
  const ScalarType logAccelerationVariancePriorScale = GetLogAccelerationVariancePriorScaleFactor();
  const ScalarType
      newLogAccelerationVariance = (S7 + logAccelerationVariancePriorDof * logAccelerationVariancePriorScale)
      / (temperature * (numberOfSubjects + logAccelerationVariancePriorDof));
  SetLogAccelerationVariance(newLogAccelerationVariance);

  /// Noise variance update.
  const unsigned int totalNumberOfObservations = dataSet->GetTotalNumberOfObservations();
  const VectorType noiseVariancePriorDofs = GetNoiseVariancePriorDegreesOfFreedom();
  const VectorType noiseVariancePriorScales = GetNoiseVariancePriorScaleVector();
  VectorType newNoiseVariance(m_NumberOfObjects, 0.0);
  for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
    newNoiseVariance[k] = (S8[k] + noiseVariancePriorDofs[k] * noiseVariancePriorScales[k])
        / (temperature * (m_DimensionOfDiscretizedObjects[k] * totalNumberOfObservations + noiseVariancePriorDofs[k]));
  }
  SetNoiseVariance(newNoiseVariance);
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::RecoverMemorizedState() {
  m_AbsoluteTimeIncrements = m_AbsoluteTimeIncrements_Memory;
  m_MinimumAbsoluteTimeIncrement = m_MinimumAbsoluteTimeIncrement_Memory;
  m_MaximumAbsoluteTimeIncrement = m_MaximumAbsoluteTimeIncrement_Memory;

  *m_ForwardReferenceGeodesic = *m_ForwardReferenceGeodesic_Memory;
  *m_BackwardReferenceGeodesic = *m_BackwardReferenceGeodesic_Memory;

  m_ProjectedModulationMatrix = m_ProjectedModulationMatrix_Memory;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::WriteModelPredictions(const LongitudinalDataSetType *const dataSet,
                        LinearVariableMapType const &popRER,
                        LinearVariablesMapType const &indRER) const {
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

  /// Getting the needed population fixed effects.
  const MatrixListType templateData = GetTemplateData();
  const std::shared_ptr<DeformableMultiObjectType> temp = m_Template->Clone();
  temp->UpdateImageIntensityAndLandmarkPointCoordinates(templateData);
  temp->Update();

  const MatrixType controlPoints = GetControlPoints();
  const MatrixType momenta = GetMomenta();
  const MatrixType modulationMatrix = GetModulationMatrix();
  const ScalarType referenceTime = GetReferenceTime();

  /// Get the needed individual random effects.
  const std::vector<VectorType> sourcesRERs = recast<VectorType>(indRER.at("Sources"));
  const std::vector<ScalarType> timeShiftRERs = recast<ScalarType>(indRER.at("TimeShift"));
  const std::vector<ScalarType> logAccelerationRERs = recast<ScalarType>(indRER.at("LogAcceleration"));

  /// Compute the absolute times (cf. this->UpdateAbsoluteTimeIncrements, here a const version).
  const std::vector<std::vector<ScalarType>> times = dataSet->GetTimes();
  const unsigned long numberOfSubjects = times.size();
  ScalarType minimumAbsoluteTimeIncrement = 0.0;
  ScalarType maximumAbsoluteTimeIncrement = 0.0;

  std::vector<std::vector<ScalarType>> absoluteTimeIncrements(numberOfSubjects);
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    const unsigned long nbObservations_i = times[i].size();
    const ScalarType acceleration_i = std::exp(logAccelerationRERs[i]);

    absoluteTimeIncrements[i].resize(nbObservations_i);
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      absoluteTimeIncrements[i][t] = acceleration_i * (times[i][t] - timeShiftRERs[i]);
    }

    if (absoluteTimeIncrements[i][0] < minimumAbsoluteTimeIncrement) {
      minimumAbsoluteTimeIncrement = absoluteTimeIncrements[i][0];
    }
    if (absoluteTimeIncrements[i][nbObservations_i - 1] > maximumAbsoluteTimeIncrement) {
      maximumAbsoluteTimeIncrement = absoluteTimeIncrements[i][nbObservations_i - 1];
    }
  }

  /// Compute the reference geodesic (cf. this->UpdateReferenceGeodesic, here a const version).
  // Forward part.
  std::shared_ptr<DiffeosType> forwardDef = m_ForwardReferenceGeodesic->Clone();
  const ScalarType forwardGeodesicLength = maximumAbsoluteTimeIncrement;
  const int forwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * forwardGeodesicLength + 1.5;
  forwardDef->SetT0(0.0);
  forwardDef->SetTN(forwardGeodesicLength);
  forwardDef->SetNumberOfTimePoints(forwardNumberOfTimePoints);
  forwardDef->SetDeformableMultiObject(temp);
  forwardDef->SetStartPositions(controlPoints);
  forwardDef->SetStartMomentas(momenta);
  forwardDef->Update();

  // Backward part.
  std::shared_ptr<DiffeosType> backwardDef = m_BackwardReferenceGeodesic->Clone();
  const ScalarType backwardGeodesicLength = -minimumAbsoluteTimeIncrement;
  const int backwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * backwardGeodesicLength + 1.5;
  backwardDef->SetT0(0.0);
  backwardDef->SetTN(backwardGeodesicLength);
  backwardDef->SetNumberOfTimePoints(backwardNumberOfTimePoints);
  backwardDef->SetDeformableMultiObject(temp);
  backwardDef->SetStartPositions(controlPoints);
  backwardDef->SetStartMomentas(-momenta);
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
    const ScalarType time = referenceTime + t * forwardStepSize;
    const std::string roundedTime = std::to_string((unsigned int) (time * 100.0 + 0.5));
    const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
    const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

    std::vector<std::string> filenames(m_NumberOfObjects);
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      std::ostringstream filename;
      filename << outputDir << Superclass::m_Name << "_";
      if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
      filename << "_ReferenceGeodesic__tp_" << clampedBackwardNumberOfTimePoints + t - 1 << "__age_";
      if (time < 0) { filename << "-"; }
      filename << floor << "." << decimals << m_TemplateObjectsNameExtension[k];
      filenames[k] = filename.str();

      if (writeFullTrajectoriesFlag) {
        std::ostringstream oss1, oss2;
        oss1 << outputDir << Superclass::m_Name;
        oss2 << outputDir << Superclass::m_Name;
        oss1 << "__ReferenceGeodesic__CP__tp_" << clampedBackwardNumberOfTimePoints + t - 1 << "__age_";
        oss2 << "__ReferenceGeodesic__MOM__tp_" << clampedBackwardNumberOfTimePoints + t - 1 << "__age_";
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
    const ScalarType time = referenceTime - t * backwardStepSize;
    const std::string roundedTime = std::to_string((unsigned int) (time * 100.0 + 0.5));
    const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
    const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

    std::vector<std::string> filenames(m_NumberOfObjects);
    for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
      std::ostringstream filename;
      filename << outputDir << Superclass::m_Name << "_";
      if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
      filename << "_ReferenceGeodesic__tp_" << clampedBackwardNumberOfTimePoints - 1 - t << "__age_";
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
      oss1 << "__ReferenceGeodesic__CP__tp_" << clampedBackwardNumberOfTimePoints - 1 - t << "__age_";
      oss2 << "__ReferenceGeodesic__MOM__tp_" << clampedBackwardNumberOfTimePoints - 1 - t << "__age_";
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

  /// Compute the projected modulation matrix (cf. this->UpdateProjectedModulationMatrix, here a const version).
  MatrixType projectedModulationMatrix = modulationMatrix;
  if (momenta.sum_of_squares() > 1e-20) {
    // Instantiate the kernel object.
    KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernel = kFactory->CreateKernelObject(m_Def->GetKernelType());

    // Compute the product metric matrix times the momenta vector.
    kernel->SetKernelWidth(m_Def->GetKernelWidth());
    kernel->SetSources(controlPoints);
    kernel->SetWeights(momenta);
    const MatrixType Km = kernel->Convolve(controlPoints);

    // Compute the momenta squared norm.
    const ScalarType mKm = dot_product(momenta, Km);

    // Vectorized parameters.
    const VectorType KmVectorized = Km.vectorize();
    const VectorType momVectorized = momenta.vectorize();

    // Loop on the columns.
    for (unsigned int k = 0; k < projectedModulationMatrix.cols(); ++k) {
      ScalarType cKm = dot_product(projectedModulationMatrix.get_column(k), KmVectorized);
      projectedModulationMatrix.set_column(k, projectedModulationMatrix.get_column(k) - momVectorized * cKm / mKm);
    }
  }

  /// Write the projected modulation matrix.
  std::ostringstream oss;
  oss << outputDir << Superclass::m_Name << "__Parameters__ProjectedModulationMatrix.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss.str().c_str(), projectedModulationMatrix);

  /// Compute and write the model-based reconstruction of the observations.
  const std::vector<std::string> subjectIds = dataSet->GetSubjectIds();
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    const unsigned int nbObservations_i = absoluteTimeIncrements[i].size();

    // Transport.
    MatrixType spaceShift = (projectedModulationMatrix * sourcesRERs[i])
        .unvectorize(m_NumberOfControlPoints, Dimension);

    // Divide the target absolute times into backward and forward targets.
    std::vector<ScalarType> forwardabsoluteTimeIncrements, backwardabsoluteTimeIncrements;
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      if (absoluteTimeIncrements[i][t] < 0.0) {
        backwardabsoluteTimeIncrements.push_back(-absoluteTimeIncrements[i][t]);
      } else { forwardabsoluteTimeIncrements.push_back(absoluteTimeIncrements[i][t]); }
    }
    std::reverse(backwardabsoluteTimeIncrements.begin(), backwardabsoluteTimeIncrements.end());

    // Perform the transport.
    MatrixListType backwardTransportedSpaceShifts = backwardDef->ParallelTransport(
        spaceShift, backwardabsoluteTimeIncrements).reverse();
    MatrixListType forwardTransportedSpaceShifts = forwardDef->ParallelTransport(
        spaceShift, forwardabsoluteTimeIncrements);

    // Concatenate the results.
    MatrixListType transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
    assert(transportedSpaceShifts.size() == nbObservations_i);

    // Exponentiation.
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      std::shared_ptr<DeformableMultiObjectType> referenceShape;
      MatrixType referenceControlPoints;

      // Cf. this->GetDeformedObjectAndControlPointsAt, here a const version.
      if (absoluteTimeIncrements[i][t] >= 0.0) {
        ScalarType continuousIndex = absoluteTimeIncrements[i][t] / forwardStepSize;
        unsigned int index = continuousIndex;
        if ((continuousIndex - index) >= 0.5) { ++index; }
        referenceShape = forwardDef->GetDeformedObjectAt(index);
        referenceControlPoints = forwardDef->GetDeformedControlPointsAt(index);
      } else {
        ScalarType continuousIndex = -absoluteTimeIncrements[i][t] / backwardStepSize;
        unsigned int index = continuousIndex;
        if ((continuousIndex - index) >= 0.5) { ++index; }
        referenceShape = backwardDef->GetDeformedObjectAt(index);
        referenceControlPoints = backwardDef->GetDeformedControlPointsAt(index);
      }

      // Riemannian exponential.
      std::shared_ptr<DiffeosType> expDef = m_Def->Clone();
      expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
      expDef->SetDeformableMultiObject(referenceShape);
      expDef->SetStartPositions(referenceControlPoints);
      expDef->SetStartMomentas(transportedSpaceShifts[t]);
      expDef->Update();

      std::shared_ptr<DeformableMultiObjectType> predictedShape = expDef->GetDeformedObject();

      // Write.
      const std::string roundedTime = std::to_string((unsigned int) (times[i][t] * 100.0 + 0.5));
      const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
      const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

      std::vector<std::string> filenames(m_NumberOfObjects);
      for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
        std::ostringstream filename;
        filename << outputDir << Superclass::m_Name << "_";
        if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
        filename << "_Prediction__subject_" << subjectIds[i]
                 << "__tp_" << t
                 << "__age_" << floor << "." << decimals
                 << m_TemplateObjectsNameExtension[k];
        filenames[k] = filename.str();
      }
      predictedShape->WriteMultiObject(filenames);
    }
  }

  /// Compute and writes the independant components.
  for (unsigned int i = 0; i < m_NumberOfSources; ++i) {
    // Construct the source vector and the associated space-shift.
    VectorType source(m_NumberOfSources, 0.0);
    source(i) = 1.0;
    MatrixType spaceShift = (projectedModulationMatrix * source).unvectorize(m_NumberOfControlPoints, Dimension);

    // Perform the transport.
    MatrixListType backwardTransportedSpaceShifts(0), forwardTransportedSpaceShifts(0);
    if (backwardNumberOfTimePoints > 1) {
      backwardTransportedSpaceShifts = backwardDef->ParallelTransport(spaceShift)
          .sublist(backwardNumberOfTimePoints - 1, 1);
    }
    if (forwardNumberOfTimePoints > 0) {
      forwardTransportedSpaceShifts = forwardDef->ParallelTransport(spaceShift);
    }

    // Concatenate the results.
    MatrixListType transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
    assert(transportedSpaceShifts.size() == clampedForwardNumberOfTimePoints + clampedBackwardNumberOfTimePoints - 1);

    // Exponentiation.
    for (unsigned int t = 0; t < clampedForwardNumberOfTimePoints + clampedBackwardNumberOfTimePoints - 1; ++t) {
      std::shared_ptr<DeformableMultiObjectType> referenceShape;
      MatrixType referenceControlPoints;
      ScalarType time;

      if (t < clampedBackwardNumberOfTimePoints - 1) {
        referenceShape = backwardDef->GetDeformedObjectAt(clampedBackwardNumberOfTimePoints - t - 1);
        referenceControlPoints = backwardDef->GetDeformedControlPointsAt(clampedBackwardNumberOfTimePoints - t - 1);
        time = referenceTime - backwardDef->GetTN() + t * backwardStepSize;
      } else {
        referenceShape = forwardDef->GetDeformedObjectAt(t - clampedBackwardNumberOfTimePoints + 1);
        referenceControlPoints = forwardDef->GetDeformedControlPointsAt(t - clampedBackwardNumberOfTimePoints + 1);
        time = referenceTime + (t - clampedBackwardNumberOfTimePoints + 1) * forwardStepSize;
      }

      // Riemannian exponential.
      std::shared_ptr<DiffeosType> expDef = m_Def->Clone();
      expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
      expDef->SetDeformableMultiObject(referenceShape);
      expDef->SetStartPositions(referenceControlPoints);
      expDef->SetStartMomentas(transportedSpaceShifts[t]);
      expDef->Update();

      std::shared_ptr<DeformableMultiObjectType> predictedShape = expDef->GetDeformedObject();

      // Write.
      const std::string roundedTime = std::to_string((unsigned int) (time * 100.0 + 0.5));
      const std::string floor = roundedTime.substr(0, roundedTime.size() - 2);
      const std::string decimals = roundedTime.substr(roundedTime.size() - 2);

      std::vector<std::string> filenames(m_NumberOfObjects);
      for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
        std::ostringstream filename;
        filename << outputDir << Superclass::m_Name << "_";
        if (m_NumberOfObjects > 1) { filename << "_" << m_TemplateObjectsName[k] << "_"; }
        filename << "_IndependentComponent_" << i << "__tp_" << t << "__age_";
        if (time < 0) { filename << "-"; }
        filename << floor << "." << decimals
                 << m_TemplateObjectsNameExtension[k];
        filenames[k] = filename.str();
      }
      predictedShape->WriteMultiObject(filenames);

      // Optionnally writes the transported momentas along the geodesic.
      if (writeFullTrajectoriesFlag) {
        std::ostringstream oss;
        oss << outputDir << Superclass::m_Name;
        oss << "__IndependentComponent_" << i;
        oss << "__TransportedMomenta__tp_" << t << "__age_";
        if (time < 0) { oss << "-"; }
        oss << floor << "." << decimals << ".txt" << std::ends;
        writeMatrixDLM<ScalarType>(oss.str().c_str(), transportedSpaceShifts[t]);
      }
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::WriteModelParameters(LinearVariablesMapType const &indRER) const {
  const std::string outputDir = def::utils::settings.output_dir;

  // Write template.
  std::vector<std::string> filenames(m_NumberOfObjects);
  for (unsigned int k = 0; k < m_NumberOfObjects; ++k) {
    std::ostringstream filename;
    filename << outputDir << Superclass::m_Name << "__Parameters__Template";
    if (m_NumberOfObjects > 1) { filename << "__" << m_TemplateObjectsName[k]; }
    filename << m_TemplateObjectsNameExtension[k];
    filenames[k] = filename.str();
  }
  GetTemplate()->WriteMultiObject(filenames);

  // Write template data.
  std::ostringstream oss1;
  oss1 << outputDir << Superclass::m_Name << "__Parameters__TemplateData.txt" << std::ends;
  writeMultipleMatrixDLM<ScalarType>(oss1.str().c_str(), GetTemplateData());

  // Write control points.
  std::ostringstream oss2;
  oss2 << outputDir << Superclass::m_Name << "__Parameters__ControlPoints.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss2.str().c_str(), GetControlPoints());

  // Write momenta.
  std::ostringstream oss3;
  oss3 << outputDir << Superclass::m_Name << "__Parameters__Momenta.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss3.str().c_str(), GetMomenta());

  // Write modulation matrix.
  std::ostringstream oss4;
  oss4 << outputDir << Superclass::m_Name << "__Parameters__ModulationMatrix.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss4.str().c_str(), GetModulationMatrix());

  // Write reference time.
  std::ostringstream oss5;
  oss5 << outputDir << Superclass::m_Name << "__Parameters__ReferenceTime.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss5.str().c_str(), MatrixType(1, 1, GetReferenceTime()));

  // Write time-shift standard deviation.
  std::ostringstream oss6;
  oss6 << outputDir << Superclass::m_Name << "__Parameters__TimeShiftStd.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss6.str().c_str(), MatrixType(1, 1, std::sqrt(GetTimeShiftVariance())));

  // Write log-acceleration standard deviation.
  std::ostringstream oss7;
  oss7 << outputDir << Superclass::m_Name << "__Parameters__LogAccelerationStd.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss7.str().c_str(), MatrixType(1, 1, std::sqrt(GetLogAccelerationVariance())));

  // Write noise standard deviation.
  std::ostringstream oss8;
  oss8 << outputDir << Superclass::m_Name << "__Parameters__NoiseStd.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss8.str().c_str(), MatrixType(GetNoiseVariance().sqrt()));

  // Write sources.
  MatrixType sourcesRERs = recast<VectorType>(indRER.at("Sources"));
  std::ostringstream oss9;
  oss9 << outputDir << Superclass::m_Name << "__Parameters__Sources.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss9.str().c_str(), sourcesRERs);

  // Write log-acceleration factors.
  MatrixType logAccelerationRERs = recast<ScalarType>(indRER.at("LogAcceleration"));
  std::ostringstream oss10;
  oss10 << outputDir << Superclass::m_Name << "__Parameters__LogAccelerations.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss10.str().c_str(), logAccelerationRERs);

  // Write log-acceleration factors.
  MatrixType timeShiftRERs = recast<ScalarType>(indRER.at("TimeShift"));
  std::ostringstream oss11;
  oss11 << outputDir << Superclass::m_Name << "__Parameters__TimeShifts.txt" << std::ends;
  writeMatrixDLM<ScalarType>(oss11.str().c_str(), timeShiftRERs);
}

template<class ScalarType, unsigned int Dimension>
std::vector<std::vector<ScalarType>>
LongitudinalAtlas<ScalarType, Dimension>
::ComputeResidualsForSubject(const LongitudinalDataSetType *const dataSet,
                             const LinearVariableMapType &popRER,
                             const LinearVariablesMapType &indRER,
                             const unsigned int &i) {
  /// Extract needed variables.
  std::vector<std::shared_ptr<DeformableMultiObjectType>> targets;
  std::vector<ScalarType> times;
  dataSet->GetDataForSubject(i, targets, times);
  const VectorType sourcesRER = recast<VectorType>(indRER.at("Sources")[i]);
  const unsigned int nbObservations_i = m_AbsoluteTimeIncrements[i].size();

  /// For the chosen subject, transport along the reference geodesic and then shoot at the target time-points.
  /// Transport.
  MatrixType spaceShift = (m_ProjectedModulationMatrix * sourcesRER).unvectorize(m_NumberOfControlPoints, Dimension);

  // Divide the target absolute times into backward and forward targets.
  std::vector<ScalarType> forwardabsoluteTimeIncrements, backwardabsoluteTimeIncrements;
  for (unsigned int t = 0; t < nbObservations_i; ++t) {
    if (m_AbsoluteTimeIncrements[i][t] < 0.0) {
      backwardabsoluteTimeIncrements.push_back(-m_AbsoluteTimeIncrements[i][t]);
    } else { forwardabsoluteTimeIncrements.push_back(m_AbsoluteTimeIncrements[i][t]); }
  }
  std::reverse(backwardabsoluteTimeIncrements.begin(), backwardabsoluteTimeIncrements.end());

  // Perform the transport.
  MatrixListType backwardTransportedSpaceShifts = m_BackwardReferenceGeodesic->ParallelTransport(
      spaceShift, backwardabsoluteTimeIncrements).reverse();
  MatrixListType forwardTransportedSpaceShifts = m_ForwardReferenceGeodesic->ParallelTransport(
      spaceShift, forwardabsoluteTimeIncrements);

  // Concatenate the results.
  MatrixListType transportedSpaceShifts = concatenate(backwardTransportedSpaceShifts, forwardTransportedSpaceShifts);
  assert(transportedSpaceShifts.size() == nbObservations_i);

  /// Exponentiation.
  std::vector<std::vector<ScalarType>> residuals(nbObservations_i);
  for (unsigned int t = 0; t < nbObservations_i; ++t) {
    std::shared_ptr<DeformableMultiObjectType> referenceShape;
    MatrixType referenceControlPoints;
    GetDeformedObjectAndControlPointsAt(m_AbsoluteTimeIncrements[i][t], referenceShape, referenceControlPoints);

    std::shared_ptr<DiffeosType> expDef = m_Def->Clone();
    expDef->SetNumberOfTimePoints(m_NumberOfTimePointsForExponentiation);
    expDef->SetDeformableMultiObject(referenceShape);
    expDef->SetStartPositions(referenceControlPoints);
    expDef->SetStartMomentas(transportedSpaceShifts[t]);
    expDef->Update();

    std::shared_ptr<DeformableMultiObjectType> predictedShape = expDef->GetDeformedObject();
    residuals[t] = predictedShape->ComputeMatch(targets[t]);
  }

  return residuals;
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::UpdateAbsoluteTimeIncrements(const LongitudinalDataSetType *const dataSet,
                               const LinearVariableMapType &popRER,
                               const LinearVariablesMapType &indRER) {
  /// Memorize the current state.
  m_AbsoluteTimeIncrements_Memory = m_AbsoluteTimeIncrements;
  m_MinimumAbsoluteTimeIncrement_Memory = m_MinimumAbsoluteTimeIncrement;
  m_MaximumAbsoluteTimeIncrement_Memory = m_MaximumAbsoluteTimeIncrement;

  /// Extract the needed population random effects realizations (RER).
  const std::vector<ScalarType> timeShiftRERs = recast<ScalarType>(indRER.at("TimeShift"));
  const std::vector<ScalarType> logAccelerationRERs = recast<ScalarType>(indRER.at("LogAcceleration"));

  /// Initialization.
  const std::vector<std::vector<ScalarType>> times = dataSet->GetTimes();
  const unsigned long numberOfSubjects = times.size();

  m_MinimumAbsoluteTimeIncrement = 0.0;
  m_MaximumAbsoluteTimeIncrement = 0.0;

  /// Looping over the subjects.
  m_AbsoluteTimeIncrements.resize(numberOfSubjects);
  for (unsigned int i = 0; i < numberOfSubjects; ++i) {
    const unsigned long nbObservations_i = times[i].size();
    const ScalarType acceleration_i = std::exp(logAccelerationRERs[i]);

    m_AbsoluteTimeIncrements[i].resize(nbObservations_i);
    for (unsigned int t = 0; t < nbObservations_i; ++t) {
      m_AbsoluteTimeIncrements[i][t] = acceleration_i * (times[i][t] - timeShiftRERs[i]);
    }

    if (m_AbsoluteTimeIncrements[i][0] < m_MinimumAbsoluteTimeIncrement) {
      m_MinimumAbsoluteTimeIncrement = m_AbsoluteTimeIncrements[i][0];
    }
    if (m_AbsoluteTimeIncrements[i][nbObservations_i - 1] > m_MaximumAbsoluteTimeIncrement) {
      m_MaximumAbsoluteTimeIncrement = m_AbsoluteTimeIncrements[i][nbObservations_i - 1];
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::UpdateReferenceGeodesic(const LinearVariableMapType &popRER) {
  /// Memorize the current state.
  *m_ForwardReferenceGeodesic_Memory = *m_ForwardReferenceGeodesic;
  *m_BackwardReferenceGeodesic_Memory = *m_BackwardReferenceGeodesic;

  /// Extract the needed population random effects realizations (RER).
  const MatrixListType templateDataRER = recast<MatrixListType>(popRER.at("TemplateData"));
  const std::shared_ptr<DeformableMultiObjectType> templateRER = m_Template->Clone();
  templateRER->UpdateImageIntensityAndLandmarkPointCoordinates(templateDataRER);
  templateRER->Update();

  MatrixType controlPoints;
  if (!m_FreezeControlPointsFlag) { controlPoints = recast<MatrixType>(popRER.at("ControlPoints")); }
  else { controlPoints = m_FrozenControlPoints; }

  const MatrixType momentaRER = recast<MatrixType>(popRER.at("Momenta"));

  /// Shoot to obtain the reference geodesic. Flow to obtain the continuously deformed template objects.
  const ScalarType margin =
      (m_MarginOnGeodesicLength - 1.0) * (m_MaximumAbsoluteTimeIncrement - m_MinimumAbsoluteTimeIncrement) / 2.;

  // Forward.
  ScalarType forwardGeodesicLength = m_MaximumAbsoluteTimeIncrement;
  if (forwardGeodesicLength > 0) { forwardGeodesicLength += margin; }
  const int forwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * forwardGeodesicLength + 1.5;
  m_ForwardReferenceGeodesic->SetT0(0.0);
  m_ForwardReferenceGeodesic->SetTN(forwardGeodesicLength);
  m_ForwardReferenceGeodesic->SetNumberOfTimePoints(forwardNumberOfTimePoints);
  m_ForwardReferenceGeodesic->SetDeformableMultiObject(templateRER);
  m_ForwardReferenceGeodesic->SetStartPositions(controlPoints);
  m_ForwardReferenceGeodesic->SetStartMomentas(momentaRER);
  if (forwardNumberOfTimePoints > 1) {
    m_ForwardReferenceGeodesic->Update();
    assert(!m_ForwardReferenceGeodesic->OutOfBox());
  }

  // Backward.
  ScalarType backwardGeodesicLength = -m_MinimumAbsoluteTimeIncrement;
  if (backwardGeodesicLength > 0) { backwardGeodesicLength += margin; }
  const int backwardNumberOfTimePoints = m_ConcentrationOfTimePointsForReferenceGeodesic * backwardGeodesicLength + 1.5;
  m_BackwardReferenceGeodesic->SetT0(0.0);
  m_BackwardReferenceGeodesic->SetTN(backwardGeodesicLength);
  m_BackwardReferenceGeodesic->SetNumberOfTimePoints(backwardNumberOfTimePoints);
  m_BackwardReferenceGeodesic->SetDeformableMultiObject(templateRER);
  m_BackwardReferenceGeodesic->SetStartPositions(controlPoints);
  m_BackwardReferenceGeodesic->SetStartMomentas(-momentaRER);
  if (backwardNumberOfTimePoints > 1) {
    m_BackwardReferenceGeodesic->Update();
    assert(!m_BackwardReferenceGeodesic->OutOfBox());
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::UpdateProjectedModulationMatrix(const LinearVariableMapType &popRER) {
  /// Memorize the current state.
  m_ProjectedModulationMatrix_Memory = m_ProjectedModulationMatrix;

  /// Project.
  m_ProjectedModulationMatrix = recast<MatrixType>(popRER.at("ModulationMatrix"));
  if (popRER.at("Momenta").sum_of_squares() > 1e-20) {

    /// Extract the needed population random effects realizations (RER).
    MatrixType controlPoints;
    if (!m_FreezeControlPointsFlag) { controlPoints = recast<MatrixType>(popRER.at("ControlPoints")); }
    else { controlPoints = m_FrozenControlPoints; }

    const MatrixType momentaRER = recast<MatrixType>(popRER.at("Momenta"));

    /// Instantiate the kernel object.
    KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernel = kFactory->CreateKernelObject(m_Def->GetKernelType());

    /// Compute the product metric matrix times the momenta vector.
    kernel->SetKernelWidth(m_Def->GetKernelWidth());
    kernel->SetSources(controlPoints);
    kernel->SetWeights(momentaRER);
    const MatrixType Km = kernel->Convolve(controlPoints);

    /// Compute the momenta squared norm.
    const ScalarType mKm = dot_product(momentaRER, Km);

    /// Vectorized version of the momenta.
    const VectorType KmVectorized = Km.vectorize();
    const VectorType momVectorized = momentaRER.vectorize();

    /// Loop on the columns.
    for (unsigned int k = 0; k < m_ProjectedModulationMatrix.cols(); ++k) {
      ScalarType cKm = dot_product(m_ProjectedModulationMatrix.get_column(k), KmVectorized);
      m_ProjectedModulationMatrix.set_column(k, m_ProjectedModulationMatrix.get_column(k) - momVectorized * cKm / mKm);
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
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
LongitudinalAtlas<ScalarType, Dimension>
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
  m_Def->SetDataDomain(m_BoundingBox);
  m_ForwardReferenceGeodesic->SetDataDomain(m_BoundingBox);
  m_BackwardReferenceGeodesic->SetDataDomain(m_BoundingBox);
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetDataDomain(m_BoundingBox);
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeTemplateDataVariables() {
  if (m_Template == 0)
    throw std::runtime_error("A template should be set to the LongitudinalAtlas model.");

  m_Template->Update();
  m_NumberOfObjects = m_Template->GetNumberOfObjects();

  /// Sets the template data fixed effect, as well as the template data random effect mean.
  MatrixListType templateData = m_Template->GetImageIntensityAndLandmarkPointCoordinates();
  SetTemplateData(templateData);
  m_DimensionOfDiscretizedObjects = m_Template->GetDimensionOfDiscretizedObjects();
  m_TemplateDataSizeParameters = templateData.GetSizeParameters();

  /// If needed, sets the template data random effect standard deviation.
  ScalarType templateDataRandomEffectVarianceSqrt = GetTemplateDataRandomEffectVarianceSqrt();
  if (templateDataRandomEffectVarianceSqrt < 0.0) {
    ScalarType kernelWidth = m_Def->GetKernelWidth();
    std::cout << "No standard deviation given for the template data random effect : "
        "defaulting to 1/50th of the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
    templateDataRandomEffectVarianceSqrt = kernelWidth / 50.0;
    SetTemplateDataRandomEffectVarianceSqrt(templateDataRandomEffectVarianceSqrt);
  }

  /// Sets the template data prior mean, as the initial template data.
  SetTemplateDataPriorMean(templateData);

  /// If needed, sets the template data prior standard deviation.
  ScalarType templateDataPriorVarianceSqrt = GetTemplateDataPriorVarianceSqrt();
  if (templateDataPriorVarianceSqrt < 0.0) {
    ScalarType kernelWidth = m_Def->GetKernelWidth();
    std::cout << "No standard deviation given for the template data prior : "
        "defaulting to the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
    templateDataPriorVarianceSqrt = kernelWidth;
    SetTemplateDataPriorVarianceSqrt(templateDataPriorVarianceSqrt);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeControlPointsVariables() {

  /// If needed, generates a regular lattice of control points.
  MatrixType controlPoints = GetControlPoints();
  if (controlPoints.rows() == 0) {

    /// Sets the spatial step.
    if (m_CPSpacing == 0.0) {
      m_CPSpacing = m_Def->GetKernelWidth();
      std::cout << "No initial CP spacing given: using diffeo kernel width of " << m_CPSpacing << std::endl;
    }

    /// Generate a regular lattice.
    const std::shared_ptr<const DeformableMultiObjectType> temp = GetTemplate();
    VectorType Xmin = temp->GetBoundingBox().get_column(0);
    VectorType Xmax = temp->GetBoundingBox().get_column(1);

    std::vector<VectorType> pointList;
    VectorType v(Dimension);
    switch (Dimension) {
      case 2: {
        ScalarType offsetX = 0.5f * (Xmax[0] - Xmin[0] - m_CPSpacing * floor((Xmax[0] - Xmin[0]) / m_CPSpacing));
        ScalarType offsetY = 0.5f * (Xmax[1] - Xmin[1] - m_CPSpacing * floor((Xmax[1] - Xmin[1]) / m_CPSpacing));

        for (ScalarType x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_CPSpacing)
          for (ScalarType y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_CPSpacing) {
            v[0] = x;
            v[1] = y;
            pointList.push_back(v);
          }
        break;
      }
      case 3: {
        ScalarType offsetX = 0.5f * (Xmax[0] - Xmin[0] - m_CPSpacing * floor((Xmax[0] - Xmin[0]) / m_CPSpacing));
        ScalarType offsetY = 0.5f * (Xmax[1] - Xmin[1] - m_CPSpacing * floor((Xmax[1] - Xmin[1]) / m_CPSpacing));
        ScalarType offsetZ = 0.5f * (Xmax[2] - Xmin[2] - m_CPSpacing * floor((Xmax[2] - Xmin[2]) / m_CPSpacing));

        for (ScalarType x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_CPSpacing)
          for (ScalarType y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_CPSpacing)
            for (ScalarType z = Xmin[2] + offsetZ; z <= Xmax[2]; z += m_CPSpacing) {
              v[0] = x;
              v[1] = y;
              v[2] = z;
              pointList.push_back(v);
            }
        break;
      }
      default:throw std::runtime_error("GenerateInitialCP not implemented in Dimensions other than 2 and 3");
        break;
    }

    unsigned int NumCPs = pointList.size();

    controlPoints.set_size(NumCPs, Dimension);
    for (unsigned int i = 0; i < NumCPs; i++)
      controlPoints.set_row(i, pointList[i]);

    SetControlPoints(controlPoints);
    std::cout << "Set of " << NumCPs << " control points defined." << std::endl;
  }

  /// Sets the dimension of tangent spaces.
  m_NumberOfControlPoints = controlPoints.rows();
  m_DimensionOfTangentSpaces = controlPoints.n_elem();

  if (!m_FreezeControlPointsFlag) {

    /// If needed, sets the control points random effect standard deviation.
    ScalarType controlPointsRandomEffectVarianceSqrt = GetControlPointsRandomEffectVarianceSqrt();
    if (controlPointsRandomEffectVarianceSqrt < 0.0) {
      ScalarType kernelWidth = m_Def->GetKernelWidth();
      std::cout << "No standard deviation given for the control points random effect : "
          "defaulting to 1/50th of the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
      controlPointsRandomEffectVarianceSqrt = kernelWidth / 50.0;
      SetControlPointsRandomEffectVarianceSqrt(controlPointsRandomEffectVarianceSqrt);
    }

    /// Sets the control points prior mean as the initial control points.
    SetControlPointsPriorMean(controlPoints);

    /// If needed, sets the control points prior standard deviation.
    ScalarType controlPointsPriorVarianceSqrt = GetControlPointsPriorVarianceSqrt();
    if (controlPointsPriorVarianceSqrt < 0.0) {
      ScalarType kernelWidth = m_Def->GetKernelWidth();
      std::cout << "No standard deviation given for the control points prior : "
          "defaulting to the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
      controlPointsPriorVarianceSqrt = kernelWidth;
      SetControlPointsPriorVarianceSqrt(controlPointsPriorVarianceSqrt);
    }
  } else {

    /// Cleans the unnecessary entries in the case of frozen control points.
    Superclass::m_FixedEffects.erase("ControlPoints");
    this->m_PopulationRandomEffects.erase("ControlPoints");
    this->m_Priors.erase("ControlPoints");
  }

}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeMomentaVariables() {
  /// Sets the momenta fixed effect and random effect mean.
  MatrixType controlPoints = GetControlPoints();
  MatrixType momenta = GetMomenta();

  if ((momenta.rows() == controlPoints.rows()) && (momenta.cols() == controlPoints.cols())) {
    std::cout << "Using a predefined set of momenta." << std::endl;
  } else {
    if (momenta.rows() > 0)
      std::cout << "Warning: initial momenta file has incompatible number of vectors. Initial momenta reset to zero."
                << std::endl;

    momenta.set_size(controlPoints.rows(), Dimension);
    momenta.fill(0.0);

    SetMomenta(momenta);
  }

  /// If needed, sets the momenta random effect standard deviation.
  ScalarType momentaEffectVarianceSqrt = GetMomentaRandomEffectVarianceSqrt();
  if (momentaEffectVarianceSqrt < 0.0) {
    ScalarType kernelWidth = m_Def->GetKernelWidth();
    std::cout << "No standard deviation given for the momenta random effect : "
        "defaulting to 1/50th of the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
    momentaEffectVarianceSqrt = kernelWidth / 50.0;
    SetMomentaRandomEffectVarianceSqrt(momentaEffectVarianceSqrt);
  }

  /// Sets the momenta prior mean as the initial momenta.
  SetMomentaPriorMean(momenta);

  /// If needed, sets the momenta prior standard deviation.
  ScalarType momentaPriorVarianceSqrt = GetMomentaPriorVarianceSqrt();
  if (momentaPriorVarianceSqrt < 0.0) {
    ScalarType kernelWidth = m_Def->GetKernelWidth();
    std::cout << "No standard deviation given for the momenta prior : "
        "defaulting to the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
    momentaPriorVarianceSqrt = kernelWidth;
    SetMomentaPriorVarianceSqrt(momentaPriorVarianceSqrt);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeModulationMatrixVariables() {
  /// If needed, sets the modulation matrix fixed effect and random effect mean.
  MatrixType modulationMatrix = GetModulationMatrix();
  if (modulationMatrix.rows() == 0) {
    std::cout << "Modulation matrix initialized to zero for " << m_NumberOfSources << "-source signal. "
              << std::endl;
    modulationMatrix.set_size(m_DimensionOfTangentSpaces, m_NumberOfSources);
    modulationMatrix.fill(0.0);
    SetModulationMatrix(modulationMatrix);
  } else {
    m_NumberOfSources = modulationMatrix.cols();
    std::cout << "Using the given modulation matrix for " << m_NumberOfSources << "-source signal." << std::endl;
  }

  /// If needed, sets the modulation matrix random effect standard deviation.
  ScalarType modulationMatrixRandomEffectVarianceSqrt = GetModulationMatrixRandomEffectVarianceSqrt();
  if (modulationMatrixRandomEffectVarianceSqrt < 0.0) {
    ScalarType kernelWidth = m_Def->GetKernelWidth();
    std::cout << "No standard deviation given for the modulation matrix random effect : "
        "defaulting to 1/50th of the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
    modulationMatrixRandomEffectVarianceSqrt = kernelWidth / 50.0;
    SetModulationMatrixRandomEffectVarianceSqrt(modulationMatrixRandomEffectVarianceSqrt);
  }

  /// Sets the modulation matrix prior mean as the initial modulation matrix.
  SetModulationMatrixPriorMean(modulationMatrix);

  /// If needed, sets the modulation matrix prior standard deviation.
  ScalarType modulationMatrixPriorVarianceSqrt = GetModulationMatrixPriorVarianceSqrt();
  if (modulationMatrixPriorVarianceSqrt < 0.0) {
    ScalarType kernelWidth = m_Def->GetKernelWidth();
    std::cout << "No standard deviation given for the modulation matrix prior : "
        "defaulting to the deformation kernel width (which is " << kernelWidth << ")." << std::endl;
    modulationMatrixPriorVarianceSqrt = kernelWidth;
    SetModulationMatrixPriorVarianceSqrt(modulationMatrixPriorVarianceSqrt);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeReferenceTimeVariables() {
  /// If needed, sets the reference time fixed effect and random effect mean.
  ScalarType referenceTime = GetReferenceTime();
  if (referenceTime != referenceTime) { // Test if isnan.
    referenceTime = 0.0;
    std::cout << "No given initial reference time. Defaulting to " << referenceTime << " ( >> ARBITRARY << )."
              << std::endl;
    SetReferenceTime(referenceTime);
  }

  /// Sets the reference time prior mean as the initial reference time.
  SetReferenceTimePriorMean(referenceTime);
  /// If needed, sets the reference time prior standard deviation.
  ScalarType referenceTimePriorVarianceSqrt = GetReferenceTimePriorVarianceSqrt();
  if (referenceTimePriorVarianceSqrt < 0.0) {
    referenceTimePriorVarianceSqrt = 1.0;
    std::cout << "No standard deviation given for the reference time prior. Defaulting to "
              << referenceTimePriorVarianceSqrt << " ( >> ARBITRARY << )." << std::endl;
    SetReferenceTimePriorVarianceSqrt(referenceTimePriorVarianceSqrt);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeSources() {
  VectorType sourcesMean(m_NumberOfSources, 0.0);
  std::static_pointer_cast<MultiScalarNormalDistributionType>(
      this->m_IndividualRandomEffects.at("Sources"))->SetMean(sourcesMean);

  ScalarType sourcesVarianceSqrt(1.0);
  std::static_pointer_cast<MultiScalarNormalDistributionType>(
      this->m_IndividualRandomEffects.at("Sources"))->SetVarianceSqrt(sourcesVarianceSqrt);
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeTimeShiftVariables() {
  /// If needed, sets the time shift variance fixed effect and random effect variance.
  ScalarType timeShiftVariance = GetTimeShiftVariance();
  if (timeShiftVariance < 0.0) {
    ScalarType timeShiftVarianceSqrt = 1.0; // 1.0 ?
    std::cout << "No initial standard deviation given for the time shift. Defaulting to "
              << timeShiftVarianceSqrt << " ( >> ARBITRARY << )." << std::endl;
    timeShiftVariance = pow(timeShiftVarianceSqrt, 2);
    SetTimeShiftVariance(timeShiftVariance);
  }

  /// Initializes the time shift random effect mean.
  ScalarType timeShiftRandomEffectMean = GetReferenceTime();
  SetTimeShiftRandomEffectMean(timeShiftRandomEffectMean);

  /// If needed, initializes the time shift variance prior scale factor.
  ScalarType timeShiftPriorScaleFactor = GetTimeShiftVariancePriorScaleFactor();
  if (timeShiftPriorScaleFactor < 0.0) {
    ScalarType timeShiftRandomEffectVarianceSqrt = GetTimeShiftRandomEffectVarianceSqrt();
    timeShiftPriorScaleFactor = pow(timeShiftRandomEffectVarianceSqrt, 2);
    std::cout << "No scale factor sqrt given for the time shift variance prior. "
        "Defaulting to the initial time shift standard deviation, i.e. "
              << timeShiftRandomEffectVarianceSqrt << " (standard)." << std::endl;
    SetTimeShiftVariancePriorScaleFactor(timeShiftPriorScaleFactor);
  }

  /// If needed, initializes the time shift variance prior degree of freedom.
  ScalarType timeShiftPriorDegreeOfFreedom = GetTimeShiftVariancePriorDegreeOfFreedom();
  if (timeShiftPriorDegreeOfFreedom < 0.0) {
    timeShiftPriorDegreeOfFreedom = 1.0;
    std::cout << "No degree of freedom given for the time shift prior. Defaulting to "
              << timeShiftPriorDegreeOfFreedom << " ( >> ARBITRARY << )." << std::endl;
    SetTimeShiftVariancePriorDegreeOfFreedom(timeShiftPriorDegreeOfFreedom);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeLogAccelerationVariables() {
  /// If needed, sets the log acceleration variance fixed effect and random effect variance.
  if (GetLogAccelerationVariance() < 0.0) {
    ScalarType logAccelerationVarianceSqrt = 0.5; // 0.1 ?
    std::cout << "No initial standard deviation given for the log acceleration. Defaulting to "
              << logAccelerationVarianceSqrt << " ( >> ARBITRARY << )." << std::endl;
    SetLogAccelerationVariance(pow(logAccelerationVarianceSqrt, 2));
  }

  /// Initializes the log acceleration random effect mean.
  VectorType logAccelerationRandomEffectMean(1, 0.0);
  SetLogAccelerationRandomEffectMean(logAccelerationRandomEffectMean);

  /// If needed, initializes the log acceleration prior scale factor.
  ScalarType logAccelerationPriorScaleFactor = GetLogAccelerationVariancePriorScaleFactor();
  if (logAccelerationPriorScaleFactor < 0.0) {
    ScalarType logAccelerationRandomEffectVarianceSqrt = GetLogAccelerationRandomEffectVarianceSqrt();
    logAccelerationPriorScaleFactor = pow(logAccelerationRandomEffectVarianceSqrt, 2);
    std::cout << "No scale factor sqrt given for the log acceleration variance prior. "
        "Defaulting to the initial log acceleration random effect standard deviation, i.e. "
              << logAccelerationRandomEffectVarianceSqrt << std::endl;
    SetLogAccelerationVariancePriorScaleFactor(logAccelerationPriorScaleFactor);
  }

  /// If needed, initializes the log acceleration prior degree of freedom.
  ScalarType logAccelerationPriorDegreeOfFreedom = GetLogAccelerationVariancePriorDegreeOfFreedom();
  if (logAccelerationPriorDegreeOfFreedom < 0.0) {
    logAccelerationPriorDegreeOfFreedom = 1.0;
    std::cout << "No degree of freedom given for the log acceleration variance prior. Defaulting to "
              << logAccelerationPriorDegreeOfFreedom << " ( >> ARBITRARY << )." << std::endl;
    SetLogAccelerationVariancePriorDegreeOfFreedom(logAccelerationPriorDegreeOfFreedom);
  }
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalAtlas<ScalarType, Dimension>
::InitializeNoiseVariables() {
  /// If needed, initializes the noise variance fixed effect.
  VectorType noiseVariance = GetNoiseVariance();
  if (noiseVariance[0] < 0.0) {
    noiseVariance.set_size(m_NumberOfObjects);
    noiseVariance.fill(1.0);
    std::cout << "No initial noise variance given. Defaulting to "
              << noiseVariance[0] << " for each object ( >> ARBITRARY << )." << std::endl;
    SetNoiseVariance(noiseVariance);
  }

  /// If needed, initializes the noise prior scale vector.
  VectorType noisePriorScaleVector = GetNoiseVariancePriorScaleVector();
  if (noisePriorScaleVector[0] < 0.0) {
    noisePriorScaleVector.set_size(m_NumberOfObjects);
    noisePriorScaleVector.fill(1.0);
    std::cout << "No scale vector given for the noise prior. Defaulting to "
              << noisePriorScaleVector[0] << " for each object ( >> ARBITRARY << )." << std::endl;
    SetNoiseVariancePriorScaleVector(noisePriorScaleVector);
  }

  /// Initializes the noise prior degrees of freedom.
  VectorType noisePriorDegreesOfFreedoms = GetNoiseVariancePriorDegreesOfFreedom();
  if (noisePriorDegreesOfFreedoms[0] < 0.0) {
    noisePriorDegreesOfFreedoms.set_size(m_NumberOfObjects);
    noisePriorDegreesOfFreedoms.fill(1.0);
    std::cout << "No degrees of freedom given for the noise prior. Defaulting to "
              << noisePriorDegreesOfFreedoms[0] << " for each object ( >> ARBITRARY << )." << std::endl;
    SetNoiseVariancePriorDegreesOfFreedom(noisePriorDegreesOfFreedoms);
  }
}

template
class LongitudinalAtlas<double, 2>;
template
class LongitudinalAtlas<double, 3>;
