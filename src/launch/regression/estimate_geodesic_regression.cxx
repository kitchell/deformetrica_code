/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "estimate_geodesic_regression.h"

/// Core files.
#include "GradientAscent.h"
#include "FastGradientAscent.h"

#include <src/io/XmlDataSet.hpp>
#if ITK_VERSION_MAJOR >= 4
#include <itkFFTWGlobalConfiguration.h>
#endif

/// Core files.
#include <src/support/utilities/SimpleTimer.h>
#include "Regression.h"

/// Non-core files.
#include "DeformableObjectReader.h"
#include "SparseDiffeoParametersXMLFile.h"

#include "ProbabilityDistributions.h"
#include "LinearAlgebra.h"

using namespace def::algebra;
using namespace def::proba;

/**
 estimateGeodesicRegression()
 Manages the estimation of sparse geodesic regression.  Reads and organizes input and handles output of results.

 @param[in] paramDiffeos An object containing diffeo parameters
 @param[in] numObservations The number of observations per object
 @param[in] numObjects The number of unique objects (eg 2 for multiple observations of left/right hemisphere)
 @param[in] paramObjectsList A vector of objects containing parameter values for each object
 @param[in] templatefnList A vector of vectors of string paths to data defining the inital template for each observation and object
 @param[in] observationfnList A vector of vectors of string paths to data defining the observations for each observation and object
 @param[in] observationTimesList A vector of vectors of times for each observation and object
 */

template<unsigned int Dimension>
void estimateGeodesicRegression(std::shared_ptr<const def::io::XmlModel> xml_model)
/*
void estimateGeodesicRegression(SparseDiffeoParameters* paramDiffeos, int numObservations, int numObjects,
								std::vector<DeformableObjectParameters::Pointer>& paramObjectsList, const std::vector<std::string>& templatefnList,
								const std::vector<std::vector<std::string> >& observationfnList, const std::vector<std::vector<double> >& observationTimesList)
*/
{
  const auto paramDiffeos = xml_model->param_diffeos;
  const auto numObjects = xml_model->num_objects();
  const auto numSubjects = xml_model->num_subjects();

  if (numSubjects != 1)
    throw std::runtime_error(std::string("Regression need only one subject"));

  xml_model->check_age_consistency();

  std::vector<DeformableObjectParameters::Pointer> paramObjectsList;
  std::vector<std::string> templatefnList;
  std::vector<std::vector<std::string>> observationfnList;
  std::vector<std::vector<double>> observationTimesList;

  observationfnList.resize(xml_model->subjects[0].visits.size());
  observationTimesList.resize(observationfnList.size());

  for (std::size_t j = 0; j < observationfnList.size(); ++j) {
    observationfnList[j].resize(xml_model->param_objects.size());
    observationTimesList[j].resize(xml_model->param_objects.size());
  }

  auto begin = xml_model->param_objects.begin();
  for (std::size_t i = 0; i < numObjects; ++i, ++begin) {

    paramObjectsList.push_back(begin->second);
    templatefnList.push_back(begin->second->GetFilename());

    auto &&files = xml_model->subjects_filename_by_object_id(begin->first);

    for (std::size_t j = 0; j < files.size(); ++j) {
      observationfnList[j][i] = files[j];
      observationTimesList[j][i] = xml_model->subjects[0].visits[j].age;
    }
  }

  const auto numObservations = xml_model->subjects[0].visits.size();

  std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
  std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
  typedef float ScalarType;
  std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif

  std::cout << "Sparse diffeomorphic regression estimation\n===" << std::endl;
  std::cout << "ITK version " << itk::Version::GetITKVersion() << std::endl;
  std::cout << "VTK version " << vtkVersion::GetVTKVersion() << std::endl;
  std::cout << "Number of objects: " << numObjects << std::endl;
  std::cout << "Number of observations: " << numObservations << std::endl;
  std::cout << "\n===\n" << std::endl;
  std::cout << "Deformation parameters:" << std::endl;
  paramDiffeos->PrintSelf(std::cout);
  std::cout << "\n===\n" << std::endl;
  for (unsigned int i = 0; i < numObjects; i++) {
    std::cout << "Object " << i << ": " << std::endl;
    std::cout << "template file: " << templatefnList[i] << std::endl;
    paramObjectsList[i]->PrintSelf(std::cout);
    std::cout << "\n===\n" << std::endl;
  }

#if ITK_VERSION_MAJOR >= 4
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

  /// Typedefs.
  typedef Diffeos<ScalarType, Dimension> DiffeosType;

  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryList;
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;

  typedef Regression<ScalarType, Dimension> RegressionType;
  typedef TimeSeriesDataSet<ScalarType, Dimension> DataSetType;

  typedef AbstractEstimator <ScalarType, Dimension> AbstractEstimatorType;
  typedef GradientAscent <ScalarType, Dimension> GradientAscentType;
  typedef FastGradientAscent<ScalarType, Dimension> FastGradientAscentType;

  /// Updates the parameter objects.
  paramDiffeos->Update();
  for (unsigned int k = 0; k < numObjects; ++k)
    paramObjectsList[k]->Update();

  bool useGradientAscent(0), useFastGradientAscent(1);
  if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientAscent") == 0) {
    useGradientAscent = 1;
    useFastGradientAscent = 0;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "SrwMhwgSaem") == 0) {
    std::cerr << "SrwMhwgSaem algorithm not compatible with the Regression model. Defaulting to fast gradient ascent."
              << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "AmalaSaem") == 0) {
    std::cerr << "AmalaSaem algorithm not compatible with the Regression model. Defaulting to fast gradient ascent."
              << std::endl;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fastgradientAscent")
      != 0) {
    std::cerr << "Unknown optimization method type: defaulting to fast gradient ascent." << std::endl;
  }

  // Set up the kernel factory - can't use any Update() before setting this up!
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
  kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

  // Create the deformation object
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetT0(paramDiffeos->GetT0());
  def->SetTN(paramDiffeos->GetTN());
  def->SetKernelWidth(paramDiffeos->GetKernelWidth());
  def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
  def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor()); // to define the bounding box
  if (not(paramDiffeos->UseImprovedEuler()))
    def->UseStandardEuler();

  if (paramDiffeos->ComputeTrueInverseFlow() == SparseDiffeoParameters::On) {
    std::cout << "Warning : the compute-true-inverse-flow integration scheme is indeed advised for image regression, "
        "but can be less accurate."
              << std::endl;
    def->SetComputeTrueInverseFlow();
  } else if (paramDiffeos->ComputeTrueInverseFlow() == SparseDiffeoParameters::Off) {
    std::cout << "Warning : for image regression, the compute-true-inverse-flow integration scheme is much faster "
        "and is usually advised."
              << std::endl;
    def->UnsetComputeTrueInverseFlow();
    def->SetRegressionFlag();
  } else {
    def->SetComputeTrueInverseFlow();
  }

  if (paramDiffeos->UseImplicitEuler()) { def->SetUseImplicitEuler(); }
  else { def->UnsetUseImplicitEuler(); }

  if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0)
    def->SetKernelType(P3M);
#ifdef USE_CUDA
    else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") == 0)
      def->SetKernelType(CUDAExact);
#endif
  else {
    if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
      std::cerr << "Unknown kernel type for the deformation : defaulting to exact." << std::endl;
    def->SetKernelType(Exact);
  }

  // Create template and target objects
  AbstractGeometryList templateObjectList(numObjects);
  std::vector<AbstractGeometryList> targetObjectList(numObservations);
  for (int s = 0; s < numObservations; s++) {
    AbstractGeometryList Aux(numObjects);
    targetObjectList[s] = Aux;
  }

  // Scan the list of objects
  typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;

  std::vector<unsigned int> timeIndices;
  timeIndices.resize(numObservations);

  for (int i = 0; i < numObjects; i++) {
    for (int s = 0; s < numObservations; s++) {
      ObjectReaderType *reader = new ObjectReaderType();
      reader->SetObjectParameters(paramObjectsList[i]);
      reader->SetFileName(observationfnList[s][i]);
      reader->Update();

      targetObjectList[s][i] = reader->GetOutput();

      // Get time point information
      ScalarType timept = observationTimesList[s][i];

      // Figure out which time index the time point corresponds to
      unsigned int T = paramDiffeos->GetNumberOfTimePoints();
      ScalarType tn = paramDiffeos->GetTN();
      ScalarType t0 = paramDiffeos->GetT0();
      int timeIndex = int((timept - t0) * ((T - 1) / (tn - t0)) + 0.5f);

      std::cout << "Time = " << timept << "\tTime index = " << timeIndex << std::endl;
      timeIndices[s] = timeIndex;
    }

    ObjectReaderType *reader = new ObjectReaderType();
    reader->SetObjectParameters(paramObjectsList[i]);
    reader->SetTemplateType();
    reader->SetFileName(templatefnList[i]);
    reader->Update();

    templateObjectList[i] = reader->GetOutput();
  }

  // Creating multi-objects for template and targets
  typename std::vector<std::shared_ptr<DeformableMultiObjectType>> target(numObservations);
  for (unsigned int t = 0; t < numObservations; ++t) {
    target[t] = std::make_shared<DeformableMultiObjectType>();
    target[t]->SetObjectList(targetObjectList[t]);
    target[t]->Update();
  }

  std::shared_ptr<DeformableMultiObjectType> templateObjects = std::make_shared<DeformableMultiObjectType>();
  templateObjects->SetObjectList(templateObjectList);
  templateObjects->Update();

  std::shared_ptr<RegressionType> model = std::make_shared<RegressionType>();

  model->SetDiffeos(def);
  model->SetTemplate(templateObjects);
  model->SetFreezeTemplateFlag(paramDiffeos->FreezeTemplate());
  model->SetFreezeControlPointsFlag(paramDiffeos->FreezeCP());
  model->SetCPSpacing(paramDiffeos->GetInitialCPSpacing());
  std::string cp_fn = paramDiffeos->GetInitialCPPosition_fn();
  if (strlen(cp_fn.c_str())) { model->SetControlPoints(cp_fn); }
  else if (paramDiffeos->OptimizeInitialControlPoints()) { model->InitializeControlPoints(true); }
  std::string mom_fn = paramDiffeos->GetInitialMomenta_fn();
  if (strlen(mom_fn.c_str())) { model->SetInitialMomenta(mom_fn); }
  model->SetSmoothingKernelWidth(paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
  model->SetWriteFullTrajectoriesFlag(paramDiffeos->WriteFullTrajectories());

  VectorType DataSigmaSquared(numObjects, 0.0);
  for (unsigned int i = 0; i < numObjects; i++) {
    ScalarType DSS = paramObjectsList[i]->GetDataSigma();
    DataSigmaSquared(i) = DSS * DSS;
  }
  model->SetDataSigmaSquared(DataSigmaSquared);

  std::string regressionName = paramDiffeos->GetModelName();
  if (regressionName.empty())
    regressionName = "Regression";

  std::vector<std::string> objectsName(numObjects);
  std::vector<std::string> objectsNameExtension(numObjects);

  for (int i = 0; i < numObjects; i++) {
    std::string sourcefn;
    sourcefn.assign(templatefnList[i]);
    int length = sourcefn.size();
    int index = length - 1;
    while ((sourcefn[index] != '.') && (index > 0))
      index--;

    if (index == 0)
      throw std::runtime_error("template file name has no extension");

    int index2 = index - 1;
    while ((sourcefn[index2] != '/' && (index2 > 0)))
      index2--;

    if (index2 > 0)
      index2 += 1;

    objectsName[i] = sourcefn.substr(index2, index - index2);
    objectsNameExtension[i] = sourcefn.substr(index, length - index);
  }

  model->SetName(regressionName);
  model->SetTemplateObjectsName(objectsName);
  model->SetTemplateObjectsNameExtension(objectsNameExtension);

  DataSetType *dataSet = new DataSetType();
  dataSet->SetDeformableMultiObjects(target);
  dataSet->SetTimeIndices(timeIndices);
  dataSet->Update();

  AbstractEstimatorType *estimator;
  if (useGradientAscent) { estimator = new GradientAscentType(); }
  else if (useFastGradientAscent) { estimator = new FastGradientAscentType(); }

  /// Final initialization of the regression model.
  model->Update();

  /// Purely informative.
  MatrixType dataDomain = model->GetBoundingBox();

  std::cout << "Working domain (only template) : \torigin = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 0) << "] \tlength = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 1) - dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 1) - dataDomain(Dimension - 1, 0) << "]" << std::endl;

  for (unsigned int t = 0; t < target.size(); ++t) {
    MatrixType BB = target[t]->GetBoundingBox();
    for (int d = 0; d < Dimension; d++) {
      dataDomain(d, 0) = (dataDomain(d, 0) < BB(d, 0) ? dataDomain(d, 0) : BB(d, 0));
      dataDomain(d, 1) = (dataDomain(d, 1) > BB(d, 1) ? dataDomain(d, 1) : BB(d, 1));
    }
  }

  std::cout << "Working domain : \torigin = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 0) << "] \tlength = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 1) - dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 1) - dataDomain(Dimension - 1, 0) << "]" << std::endl;

  estimator->SetStatisticalModel(model);
  estimator->SetDataSet(dataSet);
  estimator->SetMaxIterations(paramDiffeos->GetMaxIterations());
  estimator->SetPrintEveryNIters(paramDiffeos->GetPrintEveryNIters());
  estimator->SetSaveEveryNIters(paramDiffeos->GetSaveEveryNIters());

  if (estimator->IsGradientAscent()) {
    GradientAscentType *gradientAscentEstimator = static_cast<GradientAscentType *>(estimator);
    gradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
    gradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
    gradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
    gradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
    gradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());

  } else if (estimator->IsFastGradientAscent()) {
    FastGradientAscentType *fastGradientAscentEstimator = static_cast<FastGradientAscentType *>(estimator);
    fastGradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
    fastGradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
    fastGradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
    fastGradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
    fastGradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
  }

  /// Creates timer.
  SimpleTimer timer;

  estimator->Update();

  timer.Stop();

  std::cout << "Regression estimation took "
            << timer.GetElapsedHours() << " hours, "
            << timer.GetElapsedMinutes() << " minutes, "
            << timer.GetElapsedSeconds() << " seconds"
            << std::endl;

  delete dataSet;
  delete estimator;
}

template void estimateGeodesicRegression<2>(std::shared_ptr<const def::io::XmlModel> xml_model);
template void estimateGeodesicRegression<3>(std::shared_ptr<const def::io::XmlModel> xml_model);
