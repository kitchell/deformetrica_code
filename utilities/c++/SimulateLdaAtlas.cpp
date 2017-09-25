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
#include "LdaAtlas.h"
#include <vtkPolyDataReader.h>
#include "NonOrientedPolyLine.h"

/// Support files.
#include <src/support/utilities/GeneralSettings.h>
#include "LinearAlgebra.h"
#include <memory>
#include "MatrixDLM.h"
#include <vector>


using namespace def::algebra;

int main(int argc, char **argv) {

  if (argc != 3) {
    std::cerr << "Usage : " << argv[0]
              << "NumberOfClasses"
              << " NumberOfSubjectsPerClass"
              << std::endl;
    return -1;
  }

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif

  /// Arguments.
  unsigned int numberOfClasses = atoi(argv[1]);
  unsigned int numberOfSubjectsPerClass = atoi(argv[2]);

  /// Typedefs.
  typedef Diffeos<ScalarType, 2> DiffeosType;
  typedef LdaAtlas<ScalarType, 2> LdaAtlasType;
  typedef CrossSectionalDataSet<ScalarType, 2> DataSetType;
  typedef AbstractGeometry<ScalarType, 2> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryListType;
  typedef DeformableMultiObject<ScalarType, 2> DeformableMultiObjectType;
  typedef NonOrientedPolyLine<ScalarType, 2> NonOrientedPolyLineType;

  /// Deformation parameters.
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetT0(0.0);
  def->SetTN(1.0);
  def->SetKernelWidth(1.5);
  def->SetNumberOfTimePoints(10);
  def->SetPaddingFactor(3.); // to define the bounding box
  def->UseImprovedEuler();
  def->SetKernelType(Exact);

  /// Model parameters.
  // Template vtk file.
  std::shared_ptr<NonOrientedPolyLineType> shape = std::make_shared<NonOrientedPolyLineType>();
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName("template_shape_ground_truth.vtk");
  reader->Update();
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
  shape->SetAnatomicalCoordinateSystem("LPS");
  shape->SetPolyData(polyData);
  shape->SetKernelType(Exact);
  shape->SetKernelWidth(1.5);

  AbstractGeometryListType tempList(1);
  tempList[0] = shape;

  std::shared_ptr<DeformableMultiObjectType> temp = std::make_shared<DeformableMultiObjectType>();
  temp->SetObjectList(tempList);
  temp->Update();

  // Text files.
  MatrixType controlPoints = readMatrixDLM<ScalarType>("cp_ground_truth.txt");
  MatrixType F = readMatrixDLM<ScalarType>("F_ground_truth.txt");
  MatrixType G = readMatrixDLM<ScalarType>("G_ground_truth.txt");
  MatrixListType alpha = readMultipleMatrixDLM<ScalarType>("alpha_ground_truth.txt");

  // Hard-coded parameters.
  unsigned int intraClassPCADimension = G.cols();
  std::cout << "I detected " << intraClassPCADimension << " PGA dimensions" << std::endl;

  /// Creates the dataset object.
  DataSetType *dataSet = new DataSetType();


  /// Creates the longitudinal atlas object.
  std::shared_ptr<LdaAtlasType> model = std::make_shared<LdaAtlasType>();
  LinearVariableMapType fixedEffects;
  if (numberOfClasses>1) {
    fixedEffects["alpha"] = alpha;
    fixedEffects["F"] = F;
  }
  fixedEffects["G"] = G;
  fixedEffects["TemplateData"] = temp->GetLandmarkPoints();
  fixedEffects["ControlPoints"] = controlPoints;
  model->SetTemplate(temp);
  model->SetFixedEffects(fixedEffects);
  model->SetControlPoints(controlPoints);
  model->SetNumberOfClasses(numberOfClasses);
  model->SetNumberOfSubjectsPerClass(numberOfSubjectsPerClass);
  model->SetIntraClassPCADimension(intraClassPCADimension);
  model->SetDiffeos(def);
  model->SetName("SimulatedData");

  ///Name of objects
  std::vector<std::string> objectsName(1);
  std::vector<std::string> objectsNameExtension(1);
  objectsName[0] = "starman";
  objectsNameExtension[0] = ".vtk";
  model->SetTemplateObjectsName(objectsName);
  model->SetTemplateObjectsNameExtension(objectsNameExtension);

  ///Updating the model
  model->Update();

  // Sample.
  LinearVariableMapType popRER;
  LinearVariablesMapType indRER;
  model->Sample(dataSet, popRER, indRER);
  dataSet->Update();

  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> components = model->GetComponents();

  // Remove pre-existing files.
  def::utils::settings.output_dir = "output/";

  namespace bfs=boost::filesystem;
  const std::string dataDir = "data/";
  bfs::path dataPath(dataDir);
  for (bfs::directory_iterator end_dir_it, it(dataPath); it != end_dir_it; ++it) {
    if (it->path().string()[dataDir.size()] != '.') { bfs::remove_all(it->path()); }
  }


  // Write the betas for the generated data.
  std::vector<std::shared_ptr<DeformableMultiObjectType>> samples = dataSet->GetDeformableMultiObjects();
  std::vector<MatrixType> betas = recast<MatrixType>(indRER["beta"]);
  for (unsigned int c = 0; c < numberOfClasses; ++c) {
    for (unsigned int i=0; i<numberOfSubjectsPerClass;++i) {
      std::vector<std::string> filenames(1);
      std::ostringstream filename;
      filename << dataDir << "class_" << c << "_subject_" << i  << ".vtk";
      filenames[0] = filename.str();
      samples[numberOfSubjectsPerClass * c + i]->WriteMultiObject(filenames);
      std::ostringstream betaName;
      betaName << dataDir << "beta_class_" << c << "_subject_" << i << ".txt";
      std::string betaNameStr = betaName.str();
      writeMatrixDLM<ScalarType>(betaNameStr, betas[numberOfSubjectsPerClass * c + i]);
    }
  }

  // Write the names of the objects in the order (to correspond with betas and momenta)
  ofstream myfile (dataDir+"names.txt");
  if (myfile.is_open())
  {
    for (unsigned int c = 0; c < numberOfClasses; ++c)
      for (unsigned int i=0; i<numberOfSubjectsPerClass;++i) {
        myfile << "class_" << c << "_subject_" << i << ".vtk" << std::endl;
      }
    myfile.close();
    }
  else cout << "Unable to open file";

  ///Writing the betas.
  MatrixListType beta(betas);
  writeMultipleMatrixDLM<ScalarType>("Betas_ground_truth.txt", beta);


  std::vector<unsigned int> sizeParams(2);
  sizeParams[0] = controlPoints.rows();
  sizeParams[1] = controlPoints.cols();

  ///Writing the momenta
  MatrixListType momentas(beta.size());
  for (unsigned int c=0;c<numberOfClasses;++c)
    for (unsigned int s=0;s<numberOfSubjectsPerClass;++s){
      VectorType aux = G * betas[c*numberOfSubjectsPerClass+s].vectorize();
      momentas[c*numberOfSubjectsPerClass + s] = recast<MatrixType>(unvectorize(aux, sizeParams));
      if (numberOfClasses>1)
        momentas[c*numberOfSubjectsPerClass + s] += recast<MatrixType>(unvectorize(F*vectorize(LinearVariableType(alpha[c])), sizeParams));
    }
  writeMultipleMatrixDLM<ScalarType>("Momentas_ground_truth.txt", momentas);

  ///Writing the labels
  MatrixType labels(numberOfClasses*numberOfSubjectsPerClass,1.,0.);
  for (unsigned int c=0;c<numberOfClasses;++c)
    for (unsigned int s=0;s<numberOfSubjectsPerClass;++s)
      labels(c*numberOfSubjectsPerClass+s,0.) = c;
  writeMatrixDLM<ScalarType>("Labels.txt", labels);

  for (unsigned int i=0;i<intraClassPCADimension;++i) {
    for (unsigned int t = 0; t < components[i].size(); ++t) {
      std::vector<std::string> names(1);
      std::ostringstream oss;
      oss << dataDir << "principal_component_" << i << "_time_" << t << ".vtk" << std::ends;
      names[0] = oss.str();
      components[i][t]->WriteMultiObject(names);
    }
  }


  return EXIT_SUCCESS;
}


