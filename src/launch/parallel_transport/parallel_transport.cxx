/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "parallel_transport.h"

#include <src/support/kernels/KernelFactory.h>
#include <src/core/model_tools/deformations/Diffeos.h>
#include <src/io/SparseDiffeoParametersXMLFile.h>
#include <src/io/DeformableObjectParametersXMLFile.h>
#include <src/support/linear_algebra/LinearAlgebra.h>
#include <src/io/DeformableObjectReader.h>

#if ITK_VERSION_MAJOR >= 4
#include "itkFFTWGlobalConfiguration.h"
#endif

#include <src/io/MatrixDLM.h>
#include "MatrixDLM.h"

#include "LinearAlgebra.h"

///Core
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <src/io/XmlDataSet.hpp>
#include <src/support/utilities/GeneralSettings.h>

using namespace def::algebra;

template<unsigned int Dimension>
void parallel_transport(std::shared_ptr<const def::io::XmlModel> xml_model) {
  const auto paramDiffeos = xml_model->param_diffeos;
  const auto numObjects = xml_model->num_objects();
  const auto numSubjects = xml_model->num_subjects();

  if (numSubjects != 1)
    throw std::runtime_error(std::string("Parallel transport needs only one subject"));

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

  std::cout << "ParallelTransport\n===" << std::endl;
  std::cout << "Deformations Parameters:" << std::endl;
  paramDiffeos->PrintSelf(std::cout);
  std::cout << "\n===\n" << std::endl;

#if ITK_VERSION_MAJOR >= 4
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

  ///Typedefs.
  typedef Diffeos<ScalarType, Dimension> DiffeosType;
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryList;
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;


  //Read initial control point positions and momenta;
  auto initialCPFile = paramDiffeos->GetInitialCPPosition_fn();
  std::cout << "Reading initial control points from file : " << initialCPFile << std::endl;
  MatrixType CPRegression = readMatrixDLM<ScalarType>(initialCPFile.c_str());
  auto initialMomentaFile = paramDiffeos->GetInitialMomenta_fn();
  std::cout << "Reading initial momenta from file : " << initialMomentaFile << std::endl;
  MatrixType MomRegression = readMatrixDLM<ScalarType>(initialMomentaFile.c_str());
  //Read initial control point positions and momenta for the deformation to be transported
  auto initialCPToTransportFile = paramDiffeos->GetInitialCPPositionForTransport_fn();
  std::cout << "Reading initial control points to transport from file : " << initialCPToTransportFile << std::endl;
  MatrixType CPMatching = readMatrixDLM<ScalarType>(initialCPToTransportFile.c_str());
  auto initialMomentaToTransportFile = paramDiffeos->GetInitialMomentaForTransport_fn();
  std::cout << "Reading initial momenta to transport from file : " << initialMomentaToTransportFile << std::endl;
  MatrixType MomMatching = readMultipleMatrixDLM<ScalarType>(initialMomentaToTransportFile.c_str())[0];

  ///Sanity check
  assert(CPRegression.rows() == MomRegression.rows());
  assert(CPMatching.rows() == MomMatching.rows());

  std::vector<std::string> objectsName(numObjects);
  std::vector<std::string> objectsNameExtension(numObjects);
  std::vector<ScalarType> targetTimes;
  AbstractGeometryList templateObjectList(numObjects);


  ///Get the initial multiobject, later used to shoot and get the parallel trajectory.
  for (int i = 0; i < numObjects; i++) {
    std::shared_ptr<ObjectReaderType> reader = std::make_shared<ObjectReaderType>();
    reader->SetObjectParameters(paramObjectsList[i]);
    reader->SetTemplateType();
    reader->SetFileName(templatefnList[i]);
    std::cout << "Reading file from : " << templatefnList[i] << std::endl;
    reader->Update();
    templateObjectList[i] = reader->GetOutput();
    ///We also get the names and the extensions of the objects composing the template.
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

  std::shared_ptr<DeformableMultiObjectType> originObject = std::make_shared<DeformableMultiObjectType>();
  originObject->SetObjectList(templateObjectList);
  originObject->Update();

  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<DiffeosType> def;

  int regressionNumberOfTimePoints = paramDiffeos->GetNumberOfTimePoints();
  int matchingNumberOfTimePoints = paramDiffeos->GetMatchingNumberOfTimePoints();
  bool useExpParallelization = paramDiffeos->UseExpParallelization();
  double matchingT0 = paramDiffeos->GetMatchingT0();
  double matchingTN = paramDiffeos->GetMatchingTN();

  def = std::make_shared<DiffeosType>();
  def->SetKernelWidth(paramDiffeos->GetKernelWidth());


  ///Now the computation begins. If the user added a matching number of time points in the model.xml, then we proceed to a Schiratti style of transport : transport the regression
  ///along the matching and shoot at each point (exp-parallel trajectory). Otherwise we transport the regression and shoot once along the whole time.
  if (useExpParallelization) {
    std::cout << "We transport the matching along the regression (exp-parallelization) " << std::endl;

    //Create the deformation and set all the parameters
    def->SetT0(paramDiffeos->GetT0());
    def->SetTN(paramDiffeos->GetTN());
    def->SetNumberOfTimePoints(regressionNumberOfTimePoints);
    targetTimes.resize(regressionNumberOfTimePoints);
    for (unsigned int i = 0; i < regressionNumberOfTimePoints; ++i)
      targetTimes[i] = paramDiffeos->GetT0()
          + i * (paramDiffeos->GetTN() - paramDiffeos->GetT0()) / (regressionNumberOfTimePoints - 1);
    def->SetStartMomentas(MomRegression);
    def->SetStartPositions(CPRegression);
  }

    ///Else we use Tom Fletcher approach : transport the regression along the matching, and shoot.
  else {
    std::cout << "We transport the regression along the matching (geodesic-parallelization)" << std::endl;
    def->SetT0(matchingT0);
    def->SetTN(matchingTN);
    def->SetNumberOfTimePoints(matchingNumberOfTimePoints);
    targetTimes.resize(matchingNumberOfTimePoints);
    for (unsigned int i = 0; i < matchingNumberOfTimePoints; ++i)
      targetTimes[i] = matchingT0 + i * (matchingTN - matchingT0) / (matchingNumberOfTimePoints - 1);
    def->SetStartMomentas(MomMatching);
    def->SetStartPositions(CPMatching);
  }

  def->SetPaddingFactor(paramDiffeos->GetP3MWorkingSpacingRatio()); // we don't care if we shoot out of the box.
  if (not(paramDiffeos->UseImprovedEuler()))
    def->UseStandardEuler();
  if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0) {
    def->SetKernelType(P3M);
  }
#ifdef USE_CUDA
    else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") == 0)
{
def->SetKernelType(CUDAExact);
}
#endif
  else {
    if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
      std::cerr << "Unknown kernel type for the deformation : defaulting to exact" << std::endl;
    def->SetKernelType(Exact);
  }

  def->SetDeformableMultiObject(originObject);
  MatrixType boundingBox = originObject->GetBoundingBox();
  def->SetDataDomain(boundingBox);
  kfac->SetDataDomain(boundingBox);
  kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
  kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

  ///Updating the diffeos to get the trajectory along which we transport.
  def->Update();

  ///Only required when transporting the matching along the regression.
  std::vector<std::shared_ptr<DeformableMultiObjectType>> regressionObjects(regressionNumberOfTimePoints);
  if (useExpParallelization) {
    for (unsigned int t = 0; t < regressionNumberOfTimePoints; t++) {
      regressionObjects[t] = def->GetDeformedObjectAt(t);
    }
  }

  ///When transporting the regression along the matching, we need to shoot from the end point of the matching (it could alternatively be given)
  std::shared_ptr<DeformableMultiObjectType> endOfMatchingObject;
  if (not(useExpParallelization))
    endOfMatchingObject = def->GetDeformedObject();

  ///Velocities : it will be filled by our call to paralleltransport.
  MatrixListType velocities;
  std::cout << "Parallel Transporting" << std::endl;
  def->SetKernelType(Exact);///exact kernel required for inverse convolution.
  MatrixListType transportedMomentas =
      def->ParallelTransport(MomMatching, CPMatching, paramDiffeos->GetT0(), targetTimes, velocities);
  MatrixListType trajectoryPositions = def->GetTrajectoryPositions();
  MatrixType lastMomentas = transportedMomentas.at(transportedMomentas.size() - 1);
  MatrixType lastControlPoints = trajectoryPositions.at(trajectoryPositions.size() - 1);

  std::cout << "Shooting and Writing output Files" << std::endl;
  def = std::make_shared<DiffeosType>();
  def->SetKernelWidth(paramDiffeos->GetKernelWidth());
  def->SetPaddingFactor(10000.); // we don't care if we shoot out of the box. p3m not advised then !
  if (not(paramDiffeos->UseImprovedEuler()))
    def->UseStandardEuler();
  def->SetDataDomain(boundingBox);

  if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0) {
    def->SetKernelType(P3M);
  }
#ifdef USE_CUDA
    else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") == 0)
{
def->SetKernelType(CUDAExact);
}
#endif
  else {
    if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
      std::cerr << "Unknown kernel type for the deformation : defaulting to exact" << std::endl;
    def->SetKernelType(Exact);
  }

  const std::string outputDir = def::utils::settings.output_dir;

  if (useExpParallelization) {
    ///From each point of the regression, we shoot (assuming again the matching has 10 tps).
    for (int t = 0; t < regressionNumberOfTimePoints; t++) {
      def->SetT0(matchingT0);
      def->SetTN(matchingTN);
      def->SetNumberOfTimePoints(matchingNumberOfTimePoints);
      def->SetStartMomentas(transportedMomentas[t]);
      def->SetStartPositions(trajectoryPositions[t]);
      def->SetDeformableMultiObject(regressionObjects[t]);
      def->Update();
      std::shared_ptr<DeformableMultiObjectType> deformed = def->GetDeformedObject();

      std::vector<std::string> outputName(numObjects, "");
      for (int i = 0; i < numObjects; i++) {
        std::ostringstream oss;
        oss << outputDir << "ParallelTransport_" << objectsName[i] << "_" << "t=" << t << objectsNameExtension[i];
        outputName[i] = oss.str();
      }
      deformed->WriteMultiObject(outputName);

      ///Writing the reference geodesic
      for (int i = 0; i < numObjects; i++) {
        std::ostringstream oss;
        oss << outputDir << "Reference_geodesic" << objectsName[i] << "_" << "t=" << t << objectsNameExtension[i];
        outputName[i] = oss.str();
      }
      regressionObjects[t]->WriteMultiObject(outputName);
    }
  } else {
    def->SetStartMomentas(lastMomentas);
    def->SetStartPositions(lastControlPoints);
    def->SetT0(paramDiffeos->GetT0());
    def->SetTN(paramDiffeos->GetTN());
    def->SetNumberOfTimePoints(regressionNumberOfTimePoints);
    def->SetDeformableMultiObject(endOfMatchingObject);
    def->Update();

    ///Now we write the output objects.
    for (unsigned int t = 0; t < regressionNumberOfTimePoints; t++) {
      std::shared_ptr<DeformableMultiObjectType> deformed = def->GetDeformedObjectAt(t);

      std::vector<std::string> outputName(numObjects, "");
      for (int i = 0; i < numObjects; i++) {
        std::ostringstream oss;
        oss << "ParallelTransport_" << objectsName[i] << "_" << "t=" << t << objectsNameExtension[i];
        outputName[i] = oss.str();
      }
      deformed->WriteMultiObject(outputName);
    }

    ///Computing and writing the reference geodesic (it has not been computed before in the geo-parallelization case)
    def->SetDeformableMultiObject(originObject);
    def->SetStartMomentas(MomRegression);
    def->SetStartPositions(CPRegression);
    def->Update();

    for (unsigned int t = 0; t < regressionNumberOfTimePoints; t++) {
      std::shared_ptr<DeformableMultiObjectType> deformed = def->GetDeformedObjectAt(t);

      std::vector<std::string> outputName(numObjects, "");
      for (int i = 0; i < numObjects; i++) {
        std::ostringstream oss;
        oss << "Reference_geodesic_" << objectsName[i] << "_" << "t=" << t << objectsNameExtension[i];
        outputName[i] = oss.str();
      }
      deformed->WriteMultiObject(outputName);
    }
  }

  ///Finally, we write the momentas and controlPoints trajectory.
  for (unsigned int t = 0; t < trajectoryPositions.size(); t++) {
    //We also save the momenta/CP files :
    std::ostringstream oss;
    oss << outputDir << "Transported_ControlPoints_" << t << ".txt" << std::ends;
    writeMatrixDLM<ScalarType>(oss.str().c_str(), trajectoryPositions[t]);

    std::ostringstream oss2;
    oss2 << outputDir << "TransportedMomentas_" << t << ".txt" << std::ends;
    writeMatrixDLM<ScalarType>(oss2.str().c_str(), transportedMomentas[t]);
  }
}


template void parallel_transport<2>(std::shared_ptr<const def::io::XmlModel> xml_model);
template void parallel_transport<3>(std::shared_ptr<const def::io::XmlModel> xml_model);
