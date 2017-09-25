/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "estimate_longitudinal_atlas.h"

#include <src/core/estimators_tools/samplers/SrwMhwgSampler.h>
#include <src/support/utilities/SimpleTimer.h>
#include <src/io/XmlDataSet.hpp>

/// Core files.
#include "LongitudinalAtlas.h"
#include "McmcSaem.h"
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"

/// Support files.
#include "DeformableObjectReader.h"

#if ITK_VERSION_MAJOR >= 4
#include <itkFFTWGlobalConfiguration.h>
#endif

using namespace def::algebra;
using namespace def::proba;

template<unsigned int Dimension>
void estimateLongitudinalAtlas(std::shared_ptr<const def::io::XmlModel> xml_model) {
  const auto paramDiffeos = xml_model->param_diffeos;
  const auto numObjects = xml_model->num_objects();
  const auto numSubjects = xml_model->num_subjects();

  std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
  std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
  typedef float ScalarType;
    std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif

  std::cout << "Longitudinal diffeomorphic atlas estimation\n===" << std::endl;
  std::cout << "ITK version " << itk::Version::GetITKVersion() << std::endl;
  std::cout << "VTK version " << vtkVersion::GetVTKVersion() << std::endl;
  std::cout << "Number of objects: " << numObjects << std::endl;
  std::cout << "Number of subjects: " << numSubjects << std::endl;
  std::cout << "\n===\n" << std::endl;
  std::cout << "Deformation Parameters :" << std::endl;
  paramDiffeos->PrintSelf(std::cout);
  std::cout << "\n===\n" << std::endl;

  for (auto it : xml_model->param_objects) {
    std::cout << "Object id " << it.first << " : " << std::endl;
    DeformableObjectParameters::Pointer deformable = it.second;
    std::cout << "template file : " << deformable->GetFilename() << std::endl;
    deformable->PrintSelf(std::cout);
    std::cout << "\n===\n" << std::endl;
  }

#if ITK_VERSION_MAJOR >= 4
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

  /// Typedefs.
  typedef Diffeos<ScalarType, Dimension> DiffeosType;
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryListType;
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;
  typedef LongitudinalAtlas<ScalarType, Dimension> LongitudinalAtlasType;
  typedef LongitudinalDataSet<ScalarType, Dimension> DataSetType;
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;
  typedef McmcSaem<ScalarType, Dimension> McmcSaemType;
  typedef SrwMhwgSampler<ScalarType, Dimension> SrwMhwgSamplerType;

  /// Update the parameters' class.
  xml_model->param_diffeos->Update();
  for (auto it : xml_model->param_objects)
    it.second->Update();

  /// Check the <model-type> tag.
  assert(itksys::SystemTools::Strucmp(paramDiffeos->GetModelType().c_str(), "LongitudinalAtlas") == 0);

  /// Check the <optimization-method-type> tag.
  if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientascent") == 0) {
    std::cerr << "It is not possible to estimate a longitudinal atlas with a gradient-based algorithm. "
        "Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fastgradientascent")
      == 0) {
    std::cerr << "It is not possible to estimate a longitudinal atlas with a gradient-based algorithm. "
        "Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "amalasaem") == 0) {
    std::cerr << "It is not possible to estimate a longitudinal atlas with a gradient-based algorithm. "
        "Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "srwmhwgsaem") != 0) {
    std::cout << "Unspecified optimization method : defaulting to the SRW-MHwG-SAEM algorithm." << std::endl;
  }

  /// Set up the kernel factory.
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
  kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

  /// Checks at least two visits are available for each subject.
  if (std::any_of(xml_model->subjects.begin(),
                  xml_model->subjects.end(),
                  [](const def::io::XmlSubject &s) { return s.visits.size() < 2; }))
    throw std::runtime_error(std::string(
        "At least two visits for each subject are required to estimate a longitudinal atlas model."));

  /// Creates template and target objects.
  AbstractGeometryListType tempList(numObjects);
  std::vector<std::vector<AbstractGeometryListType>> targetsList(numSubjects);

  for (unsigned int i = 0; i < numSubjects; ++i) {
    targetsList[i].resize(xml_model->subjects[i].visits.size());
    for (unsigned int t = 0; t < xml_model->subjects[i].visits.size(); ++t) {
      AbstractGeometryListType Aux(numObjects);
      targetsList[i][t] = std::move(Aux);
    }
  }

  /// Scans the list of objects.
  std::size_t i_object = 0;
  for (auto it : xml_model->param_objects) {
    std::vector<std::vector<std::string>>
        filenames = xml_model->subjects_filename_by_object_id_several_visits(it.first);
    for (unsigned int i = 0; i < numSubjects; ++i) {
      for (unsigned int t = 0; t < xml_model->subjects[i].visits.size(); ++t) {
        ObjectReaderType reader;
        reader.SetObjectParameters(it.second);
        reader.SetFileName(filenames[i][t]);
        reader.Update();
        targetsList[i][t][i_object] = reader.GetOutput();
      }
    }

    ObjectReaderType reader;
    reader.SetObjectParameters(it.second);
    reader.SetTemplateType();
    auto fn = it.second->GetFilename();
    reader.SetFileName(fn);
    reader.Update();

    tempList[i_object++] = reader.GetOutput();
  }

  /// Scans the list of times.
  xml_model->check_age_consistency();
  std::vector<std::vector<ScalarType>> times(numSubjects);

  unsigned int count = 0;
  ScalarType meanTime = 0.0;
  ScalarType stdTime = 0.0;

  for (unsigned int i = 0; i < numSubjects; ++i) {
    times[i].resize(xml_model->subjects[i].visits.size());
    for (unsigned int t = 0; t < xml_model->subjects[i].visits.size(); ++t) {
      times[i][t] = xml_model->subjects[i].visits[t].age;

      meanTime += times[i][t];
      stdTime += times[i][t] * times[i][t];
      ++count;
    }
  }
  meanTime /= count;
  stdTime = std::sqrt(stdTime / count - meanTime * meanTime);

  /// Creates multi-objects for template and targets.
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> targets(numSubjects);
  for (unsigned int i = 0; i < numSubjects; ++i) {
    targets[i].resize(xml_model->subjects[i].visits.size());
    for (unsigned int t = 0; t < xml_model->subjects[i].visits.size(); ++t) {
      targets[i][t] = std::make_shared<DeformableMultiObjectType>();
      targets[i][t]->SetObjectList(targetsList[i][t]);
      targets[i][t]->Update();
    }
  }

  std::shared_ptr<DeformableMultiObjectType> temp = std::make_shared<DeformableMultiObjectType>();
  temp->SetObjectList(tempList);
  temp->Update();

  /// Compute the noise dimension.
  std::vector<unsigned long> noiseDimension(numObjects);
  std::vector<unsigned long> nd_temp = temp->GetDimensionOfDiscretizedObjects();
  std::vector<unsigned long> nd_targ = targets[0][0]->GetDimensionOfDiscretizedObjects();

  for (unsigned long k = 0; k < numObjects; ++k) {
    if (tempList[k]->IsOfLandmarkKind()) { // Landmark kind => current & varifold metrics.
      noiseDimension[k] = nd_temp[k];
    } else { // Image kind (LinearInterpImage or ParametricImage) => ALL SUBJECTS' DATA ARE ASSUMED SAMPLED EQUALLY
      noiseDimension[k] = nd_targ[k];
    }
  }

  /// Creates the deformation object.
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetT0(paramDiffeos->GetT0());
  def->SetTN(paramDiffeos->GetTN());
  def->SetKernelWidth(paramDiffeos->GetKernelWidth());
  def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
  def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor()); // to define the bounding box
  if (not(paramDiffeos->UseImprovedEuler()))
    def->UseStandardEuler();
  if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0) {
    def->SetKernelType(P3M);
  }
#ifdef USE_CUDA
    else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") == 0) {
      def->SetKernelType(CUDAExact);
    }
#endif
  else {
    if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
      std::cerr << "Unknown kernel type for the deformation : defaulting to exact." << std::endl;
    def->SetKernelType(Exact);
  }

  ///////////////////////////////////
  /// Creates the dataset object. ///
  ///////////////////////////////////
  DataSetType *dataSet = new DataSetType();
  dataSet->SetDeformableMultiObjects(targets);
  dataSet->SetTimes(times);
  dataSet->SetSubjectIds(xml_model->get_subject_ids());
  dataSet->Update();

  //////////////////////////////////////////////
  /// Creates the longitudinal atlas object. ///
  //////////////////////////////////////////////
  std::shared_ptr<LongitudinalAtlasType> model = std::make_shared<LongitudinalAtlasType>();
  model->SetDimensionOfDiscretizedObjects(noiseDimension);

  /// Initializes the number of sources.
  model->SetNumberOfSources(paramDiffeos->GetNumberOfSources());

  /// Initializes the time-shift variables.
  model->SetTimeShiftVariance(stdTime);

  /// Initializes the log-acceleration variables.
  ScalarType logAccelerationRandomEffectSqrt = paramDiffeos->GetLogAccelerationRandomEffectStd();
  if (logAccelerationRandomEffectSqrt > 0.0) {
    model->SetLogAccelerationVariance(pow(logAccelerationRandomEffectSqrt, 2));
  }

  /// Initializes the noise variance prior.
  VectorType noiseVariancePriorScaleVector(numObjects, 0.0);
  VectorType noiseVariancePriorDegreesOfFreedom(numObjects, 0.0);
  std::size_t i = 0;
  for (auto it : xml_model->param_objects) {
    noiseVariancePriorScaleVector(i) = it.second->GetDataSigma_Prior();
    noiseVariancePriorDegreesOfFreedom(i) =
        it.second->GetDataSigma_Normalized_Hyperparameter() * noiseDimension[i] * numSubjects;
    ++i;
  }
  model->SetNoiseVariancePriorScaleVector(noiseVariancePriorScaleVector);
  model->SetNoiseVariancePriorDegreesOfFreedom(noiseVariancePriorDegreesOfFreedom);

  /// Initializes the reference time.
  const ScalarType referenceTimePriorMean = paramDiffeos->GetReferenceTimePriorMean();
  if (referenceTimePriorMean != referenceTimePriorMean) { model->SetReferenceTime(meanTime); } // Test if isnan.
  else { model->SetReferenceTime(referenceTimePriorMean); }
  const ScalarType referenceTimePriorStd = paramDiffeos->GetReferenceTimePriorStd();
  if (referenceTimePriorStd < 0) { model->SetReferenceTimePriorVarianceSqrt(stdTime); }
  else { model->SetReferenceTimePriorVarianceSqrt(referenceTimePriorStd); }

  /// Miscellaneous parameters.
  model->SetConcentrationOfTimePointsForReferenceGeodesic(
      paramDiffeos->GetConcentrationOfTimePointsForReferenceGeodesic());
  model->SetNumberOfTimePointsForExponentiation(paramDiffeos->GetNumberOfTimePointsForExponentiation());
  model->SetMarginOnGeodesicLength(paramDiffeos->GetMarginOnGeodesicLength());

  model->SetTemplate(temp);
  model->SetFreezeControlPointsFlag(paramDiffeos->FreezeCP());
  const std::string cp_fn = paramDiffeos->GetInitialCPPosition_fn();
  model->SetControlPoints(cp_fn);
  const std::string mom_fn = paramDiffeos->GetInitialMomenta_fn();
  model->SetMomenta(mom_fn);
  model->SetCPSpacing(paramDiffeos->GetInitialCPSpacing());
  model->SetDiffeos(def);
  model->Update();

  /// Compute and set the bounding boxes.
  /* There are 3 boxes:
     1. The model one : tighlty enclosing template data and control points
     2. The diffeo one : padded version with padding equal to 0.5*PaddingFactor*DeformationKernelWidth (points should never move outside this box, used in Diffeos->CheckBoundingBox)
     3. The p3m kernel one : padded version with padding equal to PaddingFactor*DeformationKernelWidth (used to define p3MKernel grid to account for periodic boundary conditions)
  */
  MatrixType dataDomain = model->GetBoundingBox(); // include template domain and control point positions

  std::cout << "Working domain (only template) : \torigin = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 0) << "] \tlength = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 1) - dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 1) - dataDomain(Dimension - 1, 0) << "]" << std::endl;

  for (unsigned int i = 0; i < numSubjects; ++i) {
    for (unsigned int t = 0; t < targets[i].size(); ++t) {
      MatrixType BB = targets[i][t]->GetBoundingBox();
      for (unsigned int d = 0; d < Dimension; ++d) {
        dataDomain(d, 0) = (dataDomain(d, 0) < BB(d, 0) ? dataDomain(d, 0) : BB(d, 0));
        dataDomain(d, 1) = (dataDomain(d, 1) > BB(d, 1) ? dataDomain(d, 1) : BB(d, 1));
      }
    }
  }

  std::cout << "Working domain : \torigin = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 0) << "] \tlength = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 1) - dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 1) - dataDomain(Dimension - 1, 0) << "]" << std::endl;

  model->SetBoundingBox(dataDomain);

  /////////////////////////////////////
  /// Creates the estimator object. ///
  /////////////////////////////////////

  /// Instanciation.
  std::shared_ptr<SrwMhwgSamplerType> sampler = std::make_shared<SrwMhwgSamplerType>();
  sampler->SetStatisticalModel(model);
  sampler->SetDataSet(dataSet);

  McmcSaemType *estimator = new McmcSaemType();
  estimator->SetSampler(sampler);

  /// Initialization of the population random effect realizations (RER).
  LinearVariableMapType popRER;
  popRER["TemplateData"] = model->GetTemplateData();
  if (!paramDiffeos->FreezeCP()) { popRER["ControlPoints"] = model->GetControlPoints(); }
  popRER["Momenta"] = model->GetMomenta();
  popRER["ModulationMatrix"] = model->GetModulationMatrix();
  estimator->InitializePopulationRER(popRER);

  /// Initialization of the individual random effect realizations (RER).
  LinearVariablesMapType indRER;
  std::vector<LinearVariableType> sourcesRERs(numSubjects);
  std::vector<LinearVariableType> logAccelerationRERs(numSubjects);
  std::vector<LinearVariableType> timeShiftsRERs(numSubjects);
  for (unsigned int i = 0; i < numSubjects; ++i) {
    sourcesRERs[i] = VectorType(model->GetNumberOfSources(), 0.0);
    logAccelerationRERs[i] = 0.0;
    timeShiftsRERs[i] = model->GetTimeShiftRandomEffectMean();
  }
  indRER["Sources"] = sourcesRERs;
  indRER["LogAcceleration"] = logAccelerationRERs;
  indRER["TimeShift"] = timeShiftsRERs;
  estimator->InitializeIndividualRER(indRER);

  /// Initialization of the adaptive abilities.
  estimator->SetMemoryWindowSize(paramDiffeos->GetPrintAcceptanceRatesWindow());
  estimator->SetAdaptMemoryWindowSize(paramDiffeos->GetAdaptiveAcceptanceRatesWindow());
  const unsigned int arw = paramDiffeos->GetAdaptiveAcceptanceRatesWindow();
  if (arw < 1) { estimator->UnsetAdaptiveMode(); }
  else {
    estimator->SetAdaptiveMode();
    estimator->SetAdaptMemoryWindowSize(arw);
    sampler->SetAcceptanceRatesTarget(paramDiffeos->GetAdaptiveAcceptanceRatesTarget());
  }

  /// Initialization of the optional tempering.
  if (paramDiffeos->UseTempering()) {
    estimator->SetUseTempering();
    estimator->SetInitialTemperature(paramDiffeos->GetInitialTemperature());
    estimator->SetTemperingDurationRatio(paramDiffeos->GetTemperingDurationRatio());
  } else { estimator->UnsetUseTempering(); }

  /// Initialization of the sampler's proposal distributions.
  AbstractNormalDistributionMapType popPD, indPD; // PD : Proposal Distribution.
  const unsigned int proposalBlocksize = paramDiffeos->GetSrwProposalBlocksize();

  // Template data proposal distribution.
  const unsigned int templateData_nbElem = model->GetTemplateData().n_elem();
  unsigned int templateData_blocksize = templateData_nbElem; // Always all the template at once. TODO: possibility to reduce ?
  auto templateDataPD = std::make_shared<DisplacementFieldNormalDistributionType>();
  templateDataPD->SetMean(VectorType(templateData_blocksize, 0.0));
  templateDataPD->SetVarianceSqrt(paramDiffeos->GetSrwTemplateDataProposalStd());
  templateDataPD->SetKernelType(def->GetKernelType());
  ScalarType templateDataPropKernelWidth = paramDiffeos->GetSrwTemplateDataProposalKernelWidth();
  if (templateDataPropKernelWidth > 0) { templateDataPD->SetKernelWidth(templateDataPropKernelWidth); }
  else { templateDataPD->SetKernelWidth(paramDiffeos->GetKernelWidth()); }
  const unsigned int
      subsamplingStepSize = templateData_nbElem / paramDiffeos->GetSrwTemplateDataProposalNumberOfControlPoints();
  templateDataPD->SetSubsamplingStepSize(subsamplingStepSize);
  templateDataPD->SetShapeDimension(Dimension);
  popPD["TemplateData"] = std::static_pointer_cast<AbstractNormalDistributionType>(templateDataPD);

  // Control points proposal distribution.
  if (!paramDiffeos->FreezeCP()) {
    const unsigned int controlPoints_nbElem = model->GetControlPoints().n_elem();
    unsigned int controlPoints_blocksize = std::min(proposalBlocksize, controlPoints_nbElem);
    if (controlPoints_blocksize < 1) { controlPoints_blocksize = controlPoints_nbElem; }
    auto controlPointsPD = std::make_shared<MultiScalarNormalDistributionType>();
    controlPointsPD->SetMean(VectorType(controlPoints_blocksize, 0.0));
    controlPointsPD->SetVarianceSqrt(paramDiffeos->GetSrwControlPointsProposalStd());
    popPD["ControlPoints"] = std::static_pointer_cast<AbstractNormalDistributionType>(controlPointsPD);
  }

  // Momenta proposal distribution.
  const unsigned int momenta_nbElem = model->GetMomenta().n_elem();
  unsigned int momenta_blocksize = std::min(proposalBlocksize, momenta_nbElem);
  if (momenta_blocksize < 1) { momenta_blocksize = momenta_nbElem; }
  auto momentaPD = std::make_shared<MultiScalarNormalDistributionType>();
  momentaPD->SetMean(VectorType(momenta_blocksize, 0.0));
  momentaPD->SetVarianceSqrt(paramDiffeos->GetSrwMomentaProposalStd());
  popPD["Momenta"] = std::static_pointer_cast<AbstractNormalDistributionType>(momentaPD);

  // Modulation matrix proposal distribution.
  const unsigned int sources_nbElem = model->GetNumberOfSources();
  const unsigned int modulationMatrix_nbElem = model->GetModulationMatrix().n_elem();
  const unsigned int independentComponent_nbElem = modulationMatrix_nbElem / sources_nbElem;
  unsigned int modulationMatrix_blocksize = std::min(proposalBlocksize, independentComponent_nbElem);
  if (modulationMatrix_blocksize < 1) { modulationMatrix_blocksize = independentComponent_nbElem; }
  auto modulationMatrixPD = std::make_shared<MultiScalarNormalDistributionType>();
  modulationMatrixPD->SetMean(VectorType(modulationMatrix_blocksize, 0.0));
  modulationMatrixPD->SetVarianceSqrt(paramDiffeos->GetSrwModulationMatrixProposalStd());
  popPD["ModulationMatrix"] = std::static_pointer_cast<AbstractNormalDistributionType>(modulationMatrixPD);

  // Sources proposal distribution.
  unsigned int sources_blocksize = std::min(proposalBlocksize, sources_nbElem);
  if (sources_blocksize < 1) { sources_blocksize = sources_nbElem; }
  auto sourcesPD = std::make_shared<MultiScalarNormalDistributionType>();
  sourcesPD->SetMean(VectorType(sources_blocksize, 0.0));
  sourcesPD->SetVarianceSqrt(paramDiffeos->GetSrwSourcesProposalStd());
  indPD["Sources"] = std::static_pointer_cast<AbstractNormalDistributionType>(sourcesPD);

  // TimeShift proposal distribution.
  auto timeShiftPD = std::make_shared<MultiScalarNormalDistributionType>();
  timeShiftPD->SetMean(0.0);
  timeShiftPD->SetVarianceSqrt(paramDiffeos->GetSrwTimeShiftProposalStd());
  indPD["TimeShift"] = std::static_pointer_cast<AbstractNormalDistributionType>(timeShiftPD);

  // Log acceleration proposal distribution.
  auto logAccelerationPD = std::make_shared<MultiScalarNormalDistributionType>();
  logAccelerationPD->SetMean(0.0);
  logAccelerationPD->SetVarianceSqrt(paramDiffeos->GetSrwLogAccelerationProposalStd());
  indPD["LogAcceleration"] = std::static_pointer_cast<AbstractNormalDistributionType>(logAccelerationPD);

  // Finalization.
  sampler->SetPopulationProposalDistribution(popPD);
  sampler->SetIndividualProposalDistribution(indPD);

  /// Miscellaneous parameters.
  estimator->SetStatisticalModel(model);
  estimator->SetDataSet(dataSet);
  estimator->SetMaxIterations(paramDiffeos->GetMaxIterations());
  estimator->SetPrintEveryNIters(paramDiffeos->GetPrintEveryNIters());
  estimator->SetSaveEveryNIters(paramDiffeos->GetSaveEveryNIters());

  /// Creates output names for saving.
  std::string modelName = paramDiffeos->GetModelName();
  if (modelName.empty()) { modelName = "LongitudinalAtlas"; }

  std::vector<std::string> objectsName(numObjects);
  std::vector<std::string> objectsNameExtension(numObjects);

  auto template_iterator = xml_model->param_objects.begin();
  for (unsigned int i = 0; i < numObjects; ++i) {
    std::string sourcefn;
    sourcefn.assign(template_iterator->second->GetFilename());
    int length = sourcefn.size();
    int index = length - 1;
    while ((sourcefn[index] != '.') && (index > 0)) { index--; }

    if (index == 0) { throw std::runtime_error("template file name has no extension"); }

    int index2 = index - 1;
    while ((sourcefn[index2] != '/' && (index2 > 0))) { index2--; }

    if (index2 > 0) { index2 += 1; }

    objectsName[i] = sourcefn.substr(index2, index - index2);
    objectsNameExtension[i] = sourcefn.substr(index, length - index);
    ++template_iterator;
  }

  model->SetName(modelName);
  model->SetTemplateObjectsName(objectsName);
  model->SetTemplateObjectsNameExtension(objectsNameExtension);

  ///////////////
  /// Launch. ///
  ///////////////

  SimpleTimer timer;
  estimator->Update();
  timer.Stop();

  std::cout << "Longitudinal atlas estimation took "
            << timer.GetElapsedHours() << " hours, "
            << timer.GetElapsedMinutes() << " minutes, "
            << timer.GetElapsedSeconds() << " seconds"
            << std::endl;

  delete dataSet;
  delete estimator;
  KernelFactoryType::Delete();
}

template void estimateLongitudinalAtlas<2>(std::shared_ptr<const def::io::XmlModel> xml_model);
template void estimateLongitudinalAtlas<3>(std::shared_ptr<const def::io::XmlModel> xml_model);
