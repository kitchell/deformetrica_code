/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "estimate_longitudinal_registration.h"

/// Core files.
#include "LongitudinalRegistration.h"
#include "PowellsMethod.h"

/// Support files.
#include "LinearAlgebra.h"

/// IO files.
#include "DeformableObjectReader.h"

/// Miscellaneous files.
#include <src/support/utilities/SimpleTimer.h>
#include <src/io/XmlDataSet.hpp>

#if ITK_VERSION_MAJOR >= 4
#include <itkFFTWGlobalConfiguration.h>
#endif

using namespace def::algebra;
using namespace def::proba;

template<unsigned int Dimension>
void estimateLongitudinalRegistration(std::shared_ptr<const def::io::XmlModel> xml_model) {
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

  std::cout << "Longitudinal registration estimation\n===" << std::endl;
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
  typedef LongitudinalRegistration<ScalarType, Dimension> LongitudinalRegistrationType;
  typedef LongitudinalDataSet<ScalarType, Dimension> DataSetType;
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;
  typedef PowellsMethod<ScalarType, Dimension> PowellsMethodType;

  /// Update the parameters' class.
  xml_model->param_diffeos->Update();
  for (auto it : xml_model->param_objects)
    it.second->Update();

  /// Check the <model-type> tag.
  assert(itksys::SystemTools::Strucmp(paramDiffeos->GetModelType().c_str(), "LongitudinalAtlas") == 0);

  /// Check the <optimization-method-type> tag.
  if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientascent") == 0) {
    std::cerr << "It is currently not possible to estimate a longitudinal registration with another algorithm "
        "than the Powell's method. Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fastgradientascent")
      == 0) {
    std::cerr << "It is currently not possible to estimate a longitudinal registration with another algorithm "
        "than the Powell's method. Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "amalasaem") == 0) {
    std::cerr << "It is currently not possible to estimate a longitudinal registration with another algorithm "
        "than the Powell's method. Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "srwmhwgsaem") == 0) {
    std::cerr << "It is currently not possible to estimate a longitudinal registration with another algorithm "
        "than the Powell's method. Defaulting to the SRWMH-SAEM algorithm." << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "powellsmethod") != 0) {
    std::cout << "Unspecified optimization method : defaulting to the Powell's method." << std::endl;
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
        "At least two visits for each subject are required to estimate a longitudinal registration."));

  /// Check that a single subject is given as input.
  if (numSubjects != 1) {
    std::cerr << "Error : a single subject must be specified to estimate a longitudinal registration." << std::endl;
  }

  /// Creates template and target objects.
  AbstractGeometryListType tempList(numObjects);
  std::vector<std::vector<AbstractGeometryListType>> targetsList(numSubjects);

  targetsList[0].resize(xml_model->subjects[0].visits.size());
  for (unsigned int t = 0; t < xml_model->subjects[0].visits.size(); ++t) {
    AbstractGeometryListType Aux(numObjects);
    targetsList[0][t] = std::move(Aux);
  }

  /// Scans the list of objects.
  VectorType noiseVariance(numObjects);
  std::size_t i_object = 0;
  for (auto it : xml_model->param_objects) {
    std::vector<std::vector<std::string>>
        filenames = xml_model->subjects_filename_by_object_id_several_visits(it.first);
    for (unsigned int t = 0; t < xml_model->subjects[0].visits.size(); ++t) {
      ObjectReaderType reader;
      reader.SetObjectParameters(it.second);
      reader.SetFileName(filenames[0][t]);
      reader.Update();
      targetsList[0][t][i_object] = reader.GetOutput();
    }

    ObjectReaderType reader;
    reader.SetObjectParameters(it.second);
    reader.SetTemplateType();
    auto fn = it.second->GetFilename();
    reader.SetFileName(fn);
    reader.Update();

    tempList[i_object] = reader.GetOutput();
    noiseVariance(i_object++) = it.second->GetDataSigma();
  }

  /// Scans the list of times.
  xml_model->check_age_consistency();
  std::vector<std::vector<ScalarType>> times(numSubjects);

  times[0].resize(xml_model->subjects[0].visits.size());
  for (unsigned int t = 0; t < xml_model->subjects[0].visits.size(); ++t) {
    times[0][t] = xml_model->subjects[0].visits[t].age;
  }

  /// Creates multi-objects for template and targets.
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> targets(numSubjects);
  targets[0].resize(xml_model->subjects[0].visits.size());
  for (unsigned int t = 0; t < xml_model->subjects[0].visits.size(); ++t) {
    targets[0][t] = std::make_shared<DeformableMultiObjectType>();
    targets[0][t]->SetObjectList(targetsList[0][t]);
    targets[0][t]->Update();
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

  /////////////////////////////////////////////////////
  /// Creates the longitudinal registration object. ///
  /////////////////////////////////////////////////////
  std::shared_ptr<LongitudinalRegistrationType> model = std::make_shared<LongitudinalRegistrationType>();
  model->SetDimensionOfDiscretizedObjects(noiseDimension);

  /// Initializes the model hyperparameters.
  model->SetTemplate(temp);
  model->SetControlPoints(xml_model->param_diffeos->GetInitialCPPosition_fn());
  model->SetMomenta(xml_model->param_diffeos->GetInitialMomenta_fn());
  model->SetModulationMatrix(xml_model->param_diffeos->GetModulationMatrix_fn());
  model->SetReferenceTime(xml_model->param_diffeos->GetReferenceTime());
  model->SetTimeShiftVariance(xml_model->param_diffeos->GetTimeShiftVariance());
  model->SetLogAccelerationVariance(xml_model->param_diffeos->GetLogAccelerationVariance());
  model->SetNoiseVariance(noiseVariance);

  /// Initializes the log-acceleration variables.
  ScalarType logAccelerationRandomEffectSqrt = paramDiffeos->GetLogAccelerationRandomEffectStd();
  if (logAccelerationRandomEffectSqrt > 0.0) {
    model->SetLogAccelerationVariance(pow(logAccelerationRandomEffectSqrt, 2));
  }

  /// Miscellaneous parameters.
  model->SetConcentrationOfTimePointsForReferenceGeodesic(
      paramDiffeos->GetConcentrationOfTimePointsForReferenceGeodesic());
  model->SetNumberOfTimePointsForExponentiation(paramDiffeos->GetNumberOfTimePointsForExponentiation());
  model->SetMarginOnGeodesicLength(paramDiffeos->GetMarginOnGeodesicLength());
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

  for (unsigned int t = 0; t < targets[0].size(); ++t) {
    MatrixType BB = targets[0][t]->GetBoundingBox();
    for (unsigned int d = 0; d < Dimension; ++d) {
      dataDomain(d, 0) = (dataDomain(d, 0) < BB(d, 0) ? dataDomain(d, 0) : BB(d, 0));
      dataDomain(d, 1) = (dataDomain(d, 1) > BB(d, 1) ? dataDomain(d, 1) : BB(d, 1));
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
  PowellsMethodType *estimator = new PowellsMethodType();

  /// Initializes the uncertainty on the longitudinal registration fixed effects.
  LinearVariableMapType initialUncertainty;
  initialUncertainty["TimeShift"] = 5 * model->GetTimeShiftVariance();
  initialUncertainty["LogAcceleration"] = 5 * model->GetLogAccelerationVariance();
  initialUncertainty["Sources"] = 5;
  estimator->SetInitialUncertainty(initialUncertainty);

  /// Miscellaneous parameters.
  estimator->SetStatisticalModel(model);
  estimator->SetDataSet(dataSet);
  estimator->SetMaxIterations(paramDiffeos->GetMaxIterations());
  estimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
  estimator->SetPrintEveryNIters(paramDiffeos->GetPrintEveryNIters());
  estimator->SetSaveEveryNIters(paramDiffeos->GetSaveEveryNIters());

  /// Creates output names for saving.
  std::string modelName = paramDiffeos->GetModelName();
  if (modelName.empty()) { modelName = "LongitudinalRegistration"; }

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

  std::cout << "Longitudinal registration estimation took "
            << timer.GetElapsedHours() << " hours, "
            << timer.GetElapsedMinutes() << " minutes, "
            << timer.GetElapsedSeconds() << " seconds"
            << std::endl;

  delete dataSet;
  delete estimator;
  KernelFactoryType::Delete();
}

template void estimateLongitudinalRegistration<2>(std::shared_ptr<const def::io::XmlModel> xml_model);
template void estimateLongitudinalRegistration<3>(std::shared_ptr<const def::io::XmlModel> xml_model);
