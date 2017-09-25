/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

/// Core files.
#include "DeformableMultiObject.h"
#include "AbstractGeometry.h"
#include "Diffeos.h"

/// Input-output files.
#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"
#include "MatrixDLM.h"
#include "DeformableObjectReader.h"

/// Support files.
#include "KernelFactory.h"

/// Librairies files.
#if ITK_VERSION_MAJOR >= 4
#include "itkFFTWGlobalConfiguration.h"
#endif
#include "itksys/SystemTools.hxx"
#include <cstring>
#include <iostream>
#include <sstream>
#include <src/io/XmlDataSet.hpp>
#include <src/support/utilities/GeneralSettings.h>

template<class ScalarType, unsigned int Dimension>
void deform(SparseDiffeoParameters *paramDiffeos,
            bool useInverseFlow,
            const char *CP_fn,
            MatrixType MOM0_i,
            int numObjects,
            std::vector<DeformableObjectParameters::Pointer> paramObjectsList,
            const std::vector<std::string> objectfnList) {
  std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

  std::cout << "Shoot and Flow\n===" << std::endl;
  std::cout << "Number of objects to flow: " << numObjects << std::endl;
  std::cout << "\n===\n" << std::endl;
  std::cout << "Deformations Parameters:" << std::endl;
  paramDiffeos->PrintSelf(std::cout);
  std::cout << "\n===\n" << std::endl;
  for (int i = 0; i < numObjects; i++) {
    std::cout << "Object " << i << ": " << std::endl;
    std::cout << "Source file: " << objectfnList[i] << std::endl;
    paramObjectsList[i]->PrintSelf(std::cout);
    std::cout << "\n===\n" << std::endl;
  }

#if ITK_VERSION_MAJOR >= 4
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

  typedef Diffeos<ScalarType, Dimension> DiffeosType;
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryList;
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;


  // Read initial control point positions.
  MatrixType CP0 = readMatrixDLM<ScalarType>(CP_fn);
  std::cout << "number of CPs: " << CP0.rows() << std::endl;

  if (CP0.rows() != MOM0_i.rows())
    throw std::runtime_error("Number of CPs and Momentas mismatched for the given subject !");

  // Set up the kernel factory
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
  kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());

  // create the deformation object
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetKernelWidth(paramDiffeos->GetKernelWidth());
  def->SetT0(paramDiffeos->GetT0());
  def->SetTN(paramDiffeos->GetTN());
  def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
  def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor()); // to define the bounding box
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
  else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "compact") == 0)
  {
    def->SetKernelType(COMPACT);
  }
  else {
    if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
      std::cerr << "Unknown kernel type for the deformation : defaulting to exact" << std::endl;
    def->SetKernelType(Exact);
  }
  ///Compact kernel usage :
  if (paramDiffeos->UseFastConvolutions()) {def->SetUseFastConvolutions();}
  else {def -> UnsetUseFastConvolutions();}

  // create source object
  AbstractGeometryList objectList(numObjects);
  for (int i = 0; i < numObjects; i++) {
    typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;
    ObjectReaderType *reader = new ObjectReaderType();
    reader->SetObjectParameters(paramObjectsList[i]);

    reader->SetFileName(objectfnList[i]);
    //reader->SetTemplateType();
    reader->Update();

    objectList[i] = reader->GetOutput();
  }

  std::shared_ptr<DeformableMultiObjectType> object = std::make_shared<DeformableMultiObjectType>();
  object->SetObjectList(objectList);
  object->Update();

  // Determine the bounding box
  MatrixType DataDomain = object->GetBoundingBox();
  for (int k = 0; k < CP0.rows(); k++) {
    for (int d = 0; d < Dimension; d++) {
      DataDomain(d, 0) = (CP0(k, d) < DataDomain(d, 0) ? CP0(k, d) : DataDomain(d, 0));
      DataDomain(d, 1) = (CP0(k, d) > DataDomain(d, 1) ? CP0(k, d) : DataDomain(d, 1));
    }
  }
  // Pass DataDomain to kernel factory (for p3m kernel in the default mode) and deformation (for the definition of the "out of box" exception)
  kfac->SetDataDomain(DataDomain);
  def->SetDataDomain(DataDomain);
  std::cout << "Working domain: origin =  [" << DataDomain.get_column(0) << "] length = ["
            << DataDomain.get_column(1) - DataDomain.get_column(0) << "]" << std::endl;

  // shoot and flow
  std::cout << "Shoot and Flow" << std::endl;

  def->SetStartPositions(CP0);
  def->SetStartMomentas(MOM0_i);
  def->SetDeformableMultiObject(object);
  def->Update();

  if (def->OutOfBox()) { throw std::runtime_error("Out of box"); }

  const std::string outputDir = def::utils::settings.output_dir;


  // Write trajectories of control points and momenta
  MatrixListType TrajectoryPositions = def->GetTrajectoryPositions();
  MatrixListType TrajectoryMomentas = def->GetTrajectoryMomentas();
  for (int t = 0; t < def->GetNumberOfTimePoints(); t++) {
    std::ostringstream oss;
    oss << outputDir << "CP_t_" << t << ".txt" << std::ends;
    writeMatrixDLM<ScalarType>(oss.str().c_str(), TrajectoryPositions[t]);

    std::ostringstream oss2;
    oss2 << outputDir << "MOM_t_" << t << ".txt" << std::ends;
    writeMatrixDLM<ScalarType>(oss2.str().c_str(), TrajectoryMomentas[t]);
  }

  // Write deformation of objects
  std::cout << "Write output files" << std::endl;

  std::vector<std::string> outname(numObjects);
  std::vector<std::string> outext(numObjects);
  for (int i = 0; i < numObjects; i++) {
    std::string sourcefn;
    sourcefn.assign(objectfnList[i]);
    int length = sourcefn.size();

    int index_ext = length - 1;
    while ((sourcefn[index_ext] != '.')) { index_ext--; }
    int index_name = index_ext - 1;
    while ((sourcefn[index_name - 1] != '/') && (index_name >= 0)) { index_name--; }

    outname[i] = sourcefn.substr(index_name, index_ext - index_name);
    outext[i] = sourcefn.substr(index_ext, length - index_ext);

    std::ostringstream oss;
    oss << outputDir << outname[i] << "_flow";
    outname[i] = oss.str();
  }

  def->WriteFlow(outname, outext);
}

template<unsigned int Dimension>
void deformation(std::shared_ptr<const def::io::XmlModel> xml_model) {

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif

  const auto paramDiffeos = xml_model->param_diffeos;
  const auto numObjects = xml_model->num_objects();
  const auto numSubjects = xml_model->num_subjects();

  std::vector<DeformableObjectParameters::Pointer> paramObjectsList;
  std::vector<std::vector<std::string>> subjectfn;

  subjectfn.reserve(xml_model->param_objects.size());

  /**
   * Create objects in the old command line style
   */
  for (auto it : xml_model->param_objects) {
    std::vector<std::string> v;
    for (auto subject_filename : xml_model->subjects_filename_by_object_id(it.first)) {
      v.push_back(subject_filename);
    }
    subjectfn.push_back(std::move(v));
    paramObjectsList.push_back(it.second);
  }

  bool useInverseFlow = paramDiffeos->UseForwardShooting();
  auto CP_fn = paramDiffeos->GetInitialCPPosition_fn();
  auto Mom_fn = paramDiffeos->GetInitialMomenta_fn();

  if (CP_fn.empty()) {
    throw std::runtime_error("Missing <initial-cp-position> xml tag");
  }

  if (Mom_fn.empty()) {
    throw std::runtime_error("Missing <initial-mom-values> xml tag");
  }

  if (paramObjectsList.size() < subjectfn.size()) {
    throw std::runtime_error(
        "Mismatching number of objects and subjects : Please, remove the unused objects on <template> tag and check 'object_id' values; note: tag <filename>ignore</filename> must be present");
  }

  std::vector<std::string> objectfnList;
  std::vector<DeformableObjectParameters::Pointer> paramObjectsListOldStyle;

  /**
   * Converting objects into the old command line style
   **/
  {
    size_t i_object = 0;
    for (auto it : subjectfn) {
      for (auto ik : it) {
        objectfnList.push_back(ik);
        paramObjectsListOldStyle.push_back(paramObjectsList[i_object]);
      }
      ++i_object;
    }
  }

  if (numObjects != objectfnList.size())
    throw std::runtime_error(
        "Mismatching number of objects and subjects ; note: tag <filename>ignore</filename> must be present");

  MatrixListType MOM0 = readMultipleMatrixDLM<ScalarType>(Mom_fn.c_str());

  if (MOM0.size() > 0 && MOM0[0].rows() > 0 && MOM0[0].cols() > 0) {
    for (unsigned int subjectId = 0; subjectId < MOM0.size(); ++subjectId)
      deform<ScalarType, Dimension>(
          paramDiffeos, useInverseFlow, CP_fn.c_str(), MOM0[subjectId],
          numObjects, paramObjectsListOldStyle, objectfnList);
  } else {
    deform<ScalarType, Dimension>(
        paramDiffeos, useInverseFlow, CP_fn.c_str(), readMatrixDLM<ScalarType>(Mom_fn.c_str()),
        numObjects, paramObjectsListOldStyle, objectfnList);
  }
}

template void deform<double, 2>(SparseDiffeoParameters *paramDiffeos,
                       bool useInverseFlow,
                       const char *CP_fn,
                       MatrixType MOM0_i,
                       int numObjects,
                       std::vector<DeformableObjectParameters::Pointer> paramObjectsList,
                       const std::vector<std::string> objectfnList);
template void deform<double, 3>(SparseDiffeoParameters *paramDiffeos,
                       bool useInverseFlow,
                       const char *CP_fn,
                       MatrixType MOM0_i,
                       int numObjects,
                       std::vector<DeformableObjectParameters::Pointer> paramObjectsList,
                       const std::vector<std::string> objectfnList);
template void deformation<2>(std::shared_ptr<const def::io::XmlModel> xml_model);
template void deformation<3>(std::shared_ptr<const def::io::XmlModel> xml_model);
