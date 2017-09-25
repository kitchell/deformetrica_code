/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "DeformetricaConfig.h"

/// Core files.
#include "LongitudinalAtlas.h"
#include <NonOrientedPolyLine.h>

/// Support files.
#include <src/support/utilities/GeneralSettings.h>
#include "LinearAlgebra.h"

/// Libraries files.
#include <vtkPolyDataReader.h>
#include <memory>

using namespace def::algebra;

int main(int argc, char **argv) {

  if (argc != 5) {
    std::cerr << "Usage : " << argv[0]
              << " numberbOfSubjects"
              << " numberOfTimePoints"
              << " momentaMultiplier"
              << " modulationMatrixMultiplier"
              << std::endl;
    return -1;
  }

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif

  /// Arguments.
  unsigned int nbOfSubjects = atoi(argv[1]);
  unsigned int nbOfTimePoints = atoi(argv[2]);
  double momMultiplier = atof(argv[3]);
  double modMultiplier = atof(argv[4]);

  /// Typedefs.
  typedef Diffeos<ScalarType, 2> DiffeosType;
  typedef LongitudinalAtlas<ScalarType, 2> LongitudinalAtlasType;
  typedef LongitudinalDataSet<ScalarType, 2> DataSetType;
  typedef AbstractGeometry<ScalarType, 2> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryListType;
  typedef DeformableMultiObject<ScalarType, 2> DeformableMultiObjectType;
  typedef NonOrientedPolyLine<ScalarType, 2> NonOrientedPolyLineType;


  /// Model parameters.
  // Template vtk file.
  std::shared_ptr<NonOrientedPolyLineType> shape = std::make_shared<NonOrientedPolyLineType>();
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName("template_shape.vtk");
  reader->Update();
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
  shape->SetAnatomicalCoordinateSystem("LPS");
  shape->SetPolyData(polyData);
  shape->SetKernelType(Exact);
  shape->SetKernelWidth(1.0);

  AbstractGeometryListType tempList(1);
  tempList[0] = shape;

  std::shared_ptr<DeformableMultiObjectType> temp = std::make_shared<DeformableMultiObjectType>();
  temp->SetObjectList(tempList);
  temp->Update();

  /// Deformation parameters.
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetT0(0.0);
  def->SetTN(1.0);
  def->SetKernelWidth(1.0);
  def->SetNumberOfTimePoints(10);
  def->SetPaddingFactor(3); // to define the bounding box
  def->UseImprovedEuler();
  def->SetKernelType(Exact);
  def->SetDeformableMultiObject(temp);

  // Text files.
  MatrixType controlPoints = readMatrixDLM<ScalarType>("template_cp.txt");
  MatrixType momenta = momMultiplier * readMatrixDLM<ScalarType>("template_mom.txt");
  MatrixType modulationMatrix = modMultiplier * readMatrixDLM<ScalarType>("modmat.txt");

  // Hard-coded parameters.
  ScalarType referenceTime = 0.0;
  ScalarType timeShiftVariance = pow(1, 2);
  ScalarType logAccelerationVariance = pow(0.1, 2);
  ScalarType noiseVariance = pow(0.001, 2);

  /// Creates the dataset object.
  DataSetType *dataSet = new DataSetType();

  const ScalarType startingTime = - std::floor(nbOfTimePoints / 2);
  std::vector<std::vector<ScalarType>> times(nbOfSubjects);
  for (unsigned int i = 0; i < nbOfSubjects; ++i) {
    times[i].resize(nbOfTimePoints);
    for (unsigned int t = 0; t < nbOfTimePoints; ++t) {
      times[i][t] = startingTime + t;
    }
  }
  dataSet->SetTimes(times);

  /// Creates the longitudinal atlas object.
  std::shared_ptr<LongitudinalAtlasType> model = std::make_shared<LongitudinalAtlasType>();

  model->SetTemplate(temp);
  model->SetControlPoints(controlPoints);
  model->SetMomenta(momenta);
  model->SetModulationMatrix(modulationMatrix);
  model->SetReferenceTime(referenceTime);
  model->SetTimeShiftVariance(timeShiftVariance);
  model->SetLogAccelerationVariance(logAccelerationVariance);
  model->SetNoiseVariance(VectorType(1, noiseVariance));

  std::vector<unsigned long> noiseDimension = temp->GetDimensionOfDiscretizedObjects();
  model->SetDimensionOfDiscretizedObjects(noiseDimension);

  std::vector<std::string> objectsName(1);
  std::vector<std::string> objectsNameExtension(1);
  objectsName[0] = "starman";
  objectsNameExtension[0] = ".vtk";
  model->SetTemplateObjectsName(objectsName);
  model->SetTemplateObjectsNameExtension(objectsNameExtension);

  model->SetConcentrationOfTimePointsForReferenceGeodesic(5);
  model->SetNumberOfTimePointsForExponentiation(10);
  model->SetMarginOnGeodesicLength(1.25);
  model->SetDiffeos(def);
  model->SetName("SimulatedData");
  model->Update();

  /// Samples from the longitudinal atlas model and writes the results.
  // Sample.
  LinearVariableMapType popParams;
  LinearVariablesMapType indParams;
  model->Sample(dataSet, popParams, indParams);
  dataSet->Update();

  // Remove pre-existing files.
  def::utils::settings.output_dir = "output/";

  namespace bfs=boost::filesystem;
  const std::string dataDir = "data/";
  bfs::path dataPath(dataDir);
  for (bfs::directory_iterator end_dir_it, it(dataPath); it != end_dir_it; ++it) {
    if (it->path().string()[dataDir.size()] != '.') { bfs::remove_all(it->path()); }
  }

  // Write.
  model->Write(dataSet, popParams, indParams);

  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> samples = dataSet->GetDeformableMultiObjects();
  for (unsigned int i = 0; i < nbOfSubjects; ++i) {
    for (unsigned int t = 0; t < times[i].size(); ++t) {
      std::vector<std::string> filenames(1);
      std::ostringstream filename;
      filename << dataDir << "subject_" << i << "__tp_" << t << ".vtk";
      filenames[0] = filename.str();
      samples[i][t]->WriteMultiObject(filenames);
    }
  }

  return EXIT_SUCCESS;
}


