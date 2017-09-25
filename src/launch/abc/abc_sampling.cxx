#pragma once

#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"
#include "DeformetricaConfig.h"
#include "DeformableObjectReader.h"
#include "KernelFactory.h"
#include "Diffeos.h"
#include "MatrixDLM.h"
#include <lib/ThreadPool/ThreadPool.h>
#include <stdlib.h>
#include <src/support/utilities/GeneralSettings.h>

using namespace def::algebra;
using namespace def::proba;


///Saves a vector of strings into a txt file.
template<class elementsClass>
void writeVectorOfStrings(std::vector<elementsClass> v, std::string saveName){

  const std::string outputDir = def::utils::settings.output_dir;

  std::ofstream file(outputDir+saveName);
  for (unsigned int j = 0; j < v.size(); ++j)
    file << v[j] << std::endl;
  file.close();

}

///Returns mat where mat[i] is a vector containing the distances between objects[i] and the targets.
template<unsigned int Dimension>
MatrixType computeDistancesTo(std::vector<std::shared_ptr<DeformableMultiObject<ScalarType, Dimension>>> objects, std::vector<std::shared_ptr<DeformableMultiObject<ScalarType, Dimension>>> targets, const unsigned int numThread){

  std::cout << "Computing the distance for " << objects.size() << " objects to " << targets.size() << " targets." << std::endl;

  MatrixType distances(objects.size(), targets.size(), 0.);

  ///TODO : multithread this.

  {
    ThreadPool pool(numThread);
    for (unsigned int i = 0; i < objects.size(); ++i) {
      pool.enqueue([&, i]() {
        for (unsigned int j=0;j<targets.size();++j) {
          distances(i, j) = objects[i]->ComputeMatch(targets[j])[0];
        }
      });
    }
  }





//
//
//  ///Measure the distances
//  for (unsigned int i=0;i < objects.size();++i)
//  {
//    for (unsigned int j = 0; j < targets.size(); ++j)
//    {
//      distances(i,j) = objects[i]->ComputeMatch(targets[j])[0];
//    }
//  }

  return distances;
}



template<unsigned int Dimension>
void abc_sampling(std::shared_ptr<const def::io::XmlModel> xml_model)
{
  //Type definitions :
  typedef Diffeos<ScalarType, Dimension> DiffeosType;
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryList;
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;

  ///Some parameters
  const auto paramDiffeos = xml_model->param_diffeos;
  const auto numObjects = xml_model->num_objects();
  const auto numSubjects = xml_model->num_subjects();
  const auto numThread = paramDiffeos->GetNumberOfThreads();
  const auto saveEveryNIter = paramDiffeos->GetSaveEveryNIters();
  const auto numberOfSamples = paramDiffeos->GetNumberOfSamples();
  const auto numberOfControlPoints = paramDiffeos->GetNumberOfCps();
  const auto trainingSetSize = paramDiffeos->GetTrainingSetSize();
  const auto kernelWidth = paramDiffeos->GetKernelWidth();
  const std::string outputDir = def::utils::settings.output_dir;

  std::vector<std::string> templateObjectFileNames(numObjects, "");
  std::vector<std::string> templateObjectExtension(numObjects, "");
  std::vector<AbstractGeometryList> targetObjectList(numSubjects);
  std::vector<std::string> targetsFileNames(numSubjects);

  ///initializing the objects lists.
  for (unsigned int i = 0; i < numSubjects; ++i)
    targetObjectList[i] = AbstractGeometryList(numObjects);

  std::cout << "Starting abc procedure" << std::endl;
  std::cout << "Number of samples : " << numberOfSamples << std::endl;
  std::cout << "Number of control points : " << numberOfControlPoints << std::endl;
  std::cout << "Number of subjects : " << numSubjects << std::endl;
  std::cout << "Number of objects : " << numObjects << std::endl;
  std::cout << "Training set size : " << trainingSetSize << std::endl;
  std::cout << "Number of threads :" << numThread << std::endl;
  std::vector<DeformableObjectParameters::Pointer> paramObjectsList;
  std::vector<std::string> templatefnList;
  std::vector<std::vector<std::string>> observationfnList;
  std::vector<std::vector<double>> observationTimesList;
  std::vector<std::string> fileNames(numSubjects);
  std::vector<std::string> subjectClasses(numSubjects);
  std::vector<std::string> classesFound;
  std::vector<bool> belongsToTestSet(numSubjects);
  srand(time(NULL));

  std::cout << "Parsing the xml and building objects..." << std::endl;

  ///Parse the xml object to get the filenames and the properties of the objects and template.
  std::size_t i_object = 0;
  std::size_t i_subject = 0;

  for (auto it : xml_model->param_objects)
  {


    std::cout << "Reading files for object " << i_object << std::endl;

    ///Reading filenames and corresponding objects
    for (auto subject_filename : xml_model->subjects_filename_by_object_id(it.first)) {
      ObjectReaderType reader;
      reader.SetObjectParameters(it.second);
      reader.SetFileName(subject_filename);
      targetsFileNames[i_subject] = subject_filename;
      reader.Update();
      fileNames[i_subject] = subject_filename;
      targetObjectList[i_subject++][i_object] = reader.GetOutput();
    }

    i_subject = 0;

    ///Reading test or train info for each subject.
    for (auto subject_belongs_to_test_set : xml_model->subjects_test_or_train_label_by_object_id(it.first)) {
      belongsToTestSet[++i_subject] = subject_belongs_to_test_set;
    }


    auto fn = it.second->GetFilename();

    ///Find extension
    int index = fn.size() - 1;
    while ((fn[index] != '.') && (index > 0))
      index--;
    if (index == 0)
      throw std::runtime_error("template file name has no extension");
    ///Find name.
    int index2 = index - 1;
    while ((fn[index2] != '/' && (index2 > 0)))
      index2--;
    templateObjectFileNames[i_object] = fn.substr(index2 + 1, index - index2);
    templateObjectExtension[i_object++] = fn.substr(index + 1, fn.size() - index);

  }

  i_subject = 0;

  ///Reading class labels for each subject.
  for (auto subject_class_label : xml_model->subjects_class_label()) {
    subjectClasses[i_subject++] = subject_class_label;
    if (std::find(classesFound.begin(), classesFound.end(), subject_class_label) == classesFound.end())
      classesFound.push_back(subject_class_label);
  }

  std::cout << classesFound.size() << " classes identified in the data set." << std::endl;


  ///fill target with a deformablemultiobject for each subject.
  typename std::vector<std::shared_ptr<DeformableMultiObjectType>> target(numSubjects);
  for (unsigned int s = 0; s < numSubjects; s++)
  {
    target[s] = std::make_shared<DeformableMultiObjectType>();
    target[s]->SetObjectList(targetObjectList[s]);
    target[s]->Update();
  }


  std::cout << "Initializing diffeo object..." << std::endl;
  std::cout << ">>>>>>>>> Watch out : abc method not ready to take care of multi objects " << std::endl;
  ///Initialize the diffeo object.
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetKernelWidth(kernelWidth);
  def->SetT0(paramDiffeos->GetT0());
  def->SetTN(paramDiffeos->GetTN());
  def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
  def->SetPaddingFactor(paramDiffeos->GetP3MWorkingSpacingRatio()); // we don't care if we shoot out of the box.

  ///Setting kernel type
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


  def->SetDeformableMultiObject(target[0]);
  ///Setting the working domain.
  MatrixType boundingBox = target[0]->GetBoundingBox();
  def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());
  def->SetDataDomain(boundingBox);
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetDataDomain(boundingBox);
  kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
  kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());
  MatrixType BoundingBox = target[0]->GetBoundingBox();

  ///Probability distributions for the sampling procedure :
  UniformDistributionType uniformSampler;
  uniformSampler = UniformDistributionType(BoundingBox.transpose().get_row(0), BoundingBox.transpose().get_row(1));
  NormalDistributionType normalSampler;
  ScalarType momPropVariance = paramDiffeos->GetMomentaPropositionVariance();
  std::cout << "Momenta proposition variance" << momPropVariance << std::endl;
  MatrixType id(2, 2, 0.);
  id.set_identity();
  VectorType mean(2, 0.);
  normalSampler = NormalDistributionType();
  normalSampler.SetMean(mean);
  normalSampler.SetCovariance(id);

  unsigned int randomIntControlPoint;
  long int randomIntTemplate(-1);
  VectorType aux;
  MatrixType mom(numberOfControlPoints, Dimension, 0.);
  MatrixType cp(numberOfControlPoints, Dimension, 0.);
  std::shared_ptr<DeformableMultiObjectType> randomTemplate;
  std::shared_ptr<DeformableMultiObjectType> randomTemplateMI;
  std::vector<ScalarType> kernelWidths(numberOfSamples);
  MatrixListType momentas(numberOfSamples);
  MatrixListType controlPoints(numberOfSamples);
  std::vector<std::string> chosenTemplates(numberOfSamples,"");
  typename std::vector<std::shared_ptr<DeformableMultiObjectType>> sampledObjects(numberOfSamples);
  std::vector<std::shared_ptr<DeformableMultiObject<ScalarType, Dimension>>> trainingSet;
  std::vector<std::string> trainingSetFileNames;

  std::cout << "Starting simulation !" << std::endl;
  std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();




  ///We want n samples.
  ///We have a training set of size t.
  ///for each class
  ///for each of the first t elements e (templates)
  ///deform e n/(t*nbclass) times.
  ///This recipe only worls if n % t*nclass == 0
  assert(numberOfSamples % classesFound.size()*trainingSetSize == 0);

  unsigned int currSample(0);

  for (unsigned int i=0; i < classesFound.size();++i)
  {
    for (unsigned int t=0;t<trainingSetSize;++t)
    {

      if (currSample%10==0)
        std::cout << "Iteration " << i << std::endl;

      ++randomIntTemplate;
      ///We want to get the t elements of the class
      while (subjectClasses[randomIntTemplate]!= classesFound[i])
      {
        ++randomIntTemplate;
        if (randomIntTemplate > numSubjects)
          std::__throw_runtime_error("Training set larger that a class size ?");
      }
      randomTemplate = target[randomIntTemplate];
      trainingSet.push_back(randomTemplate);
      trainingSetFileNames.push_back(targetsFileNames[randomIntTemplate]);
      std::cout << "Chosen template :" << fileNames[randomIntTemplate] << std::endl;

      ///Deforming this template.
      for (unsigned int k=0;k< (numberOfSamples/(trainingSetSize*classesFound.size())); ++k) {

        ///Sampling parameters
        for (unsigned j = 0; j < numberOfControlPoints; ++j) {
          unsigned int nbPoints = randomTemplate->GetNumberOfPoints()[0];
          randomIntControlPoint = rand() % nbPoints;
          if (randomTemplate->IsOfLandmarkKind()[0])
            aux = randomTemplate->GetLandmarkPoints().get_row(randomIntControlPoint);
          else
            aux = uniformSampler.Sample();
          cp.set_row(j, aux);
          ///Sampling a momenta :
          aux = std::sqrt(momPropVariance) * normalSampler.Sample();
          mom.set_row(j, aux);
        }
        def->SetKernelWidth(kernelWidth);
        def->SetStartPositions(cp);
        def->SetStartMomentas(mom);
        def->SetDeformableMultiObject(randomTemplate);
        def->SetModified();
        def->Update();

        std::shared_ptr<DeformableMultiObjectType> deformedTemplate;
        deformedTemplate = def->GetDeformedObject();
        deformedTemplate->Update();
        sampledObjects[currSample] = deformedTemplate;

        std::vector<std::string> outputName(numObjects, "");
        for (int j = 0; j < numObjects; ++j) {
          std::ostringstream oss;
          oss << outputDir << templateObjectFileNames[j] << "__deformed_" << currSample << "."
              << templateObjectExtension[j];
          outputName[j] = oss.str();
        }
        deformedTemplate->WriteMultiObject(outputName);

        chosenTemplates[currSample] = fileNames[randomIntTemplate];
        momentas[currSample] = mom;
        controlPoints[currSample] = cp;
        kernelWidths[currSample] = def->GetKernelWidth();

        currSample++;
      }
    }
  }


  std::vector<std::shared_ptr<DeformableMultiObject<ScalarType, Dimension>>> testSet;
  std::vector<std::string> testSetFileNames;


  ///NOTE that the training set form now on consists only of the elements used to sample + the samples.
  for (unsigned int i=0;i<target.size();++i){
    if (belongsToTestSet[i]){
      testSet.push_back(target[i]);
      testSetFileNames.push_back(targetsFileNames[i]);
    }
  }


  ///Extending the training set with the sampled objects
  for (unsigned int i=0;i<sampledObjects.size();++i) {
    trainingSet.push_back(sampledObjects[i]);
    trainingSetFileNames.push_back("Generated_from_"+chosenTemplates[i]);
  }

  ///Computing the distances between both training set and test set to test set.
  MatrixType distancesTrainingSet = computeDistancesTo(trainingSet, trainingSet, numThread);
  MatrixType distancesTestSet = computeDistancesTo(testSet, trainingSet, numThread);

  std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
  std::cout << "It took " << duration << " seconds to get the " << numberOfSamples << " samples." << std::endl;

  //Saves the output distances for the training set:
  std::ostringstream ossDistancesTrainingSet;
  ossDistancesTrainingSet << outputDir << "DistancesTrainingSet.txt" << std::ends;
  writeMatrixDLM<ScalarType>(ossDistancesTrainingSet.str().c_str(), distancesTrainingSet);

  //Saves the output distances for the test set:
  std::ostringstream ossDistancesTestSet;
  ossDistancesTestSet << outputDir << "DistancesTestSet.txt" << std::ends;
  writeMatrixDLM<ScalarType>(ossDistancesTestSet.str().c_str(), distancesTestSet);


  ///Writing some info.
  writeVectorOfStrings<std::string>(chosenTemplates, "chosenTemplate.txt");
  writeVectorOfStrings<ScalarType>(kernelWidths, "kernelWidths_Sampled.txt");
  writeVectorOfStrings<std::string>(testSetFileNames, "testSetFileNames.txt");
  writeVectorOfStrings<std::string>(trainingSetFileNames, "trainingSetFileNames.txt");


  std::ostringstream oss2;
  oss2 << outputDir << "Momentas_Sampled.txt" << std::ends;
  writeMultipleMatrixDLM<ScalarType>(oss2.str().c_str(), momentas);
  std::ostringstream oss;
  oss << outputDir << "ControlPoints_Sampled.txt" << std::ends;
  writeMultipleMatrixDLM<ScalarType>(oss.str().c_str(), controlPoints);

}




