/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/
#include "estimate_atlas.h"
using namespace def::algebra;
using namespace def::proba;
template<unsigned int Dimension>
void estimateAtlas(std::shared_ptr<const def::io::XmlModel> xml_model) {
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
  std::cout << "Sparse diffeomorphic atlas estimation\n===" << std::endl;
  std::cout << "ITK version " << itk::Version::GetITKVersion() << std::endl;
  std::cout << "VTK version " << vtkVersion::GetVTKVersion() << std::endl;
  std::cout << "Number of objects: " << numObjects << std::endl;
  std::cout << "Number of subjects: " << numSubjects << std::endl;
  std::cout << "\n===\n" << std::endl;
  std::cout << "Deformation Parameters:" << std::endl;
  paramDiffeos->PrintSelf(std::cout);
  std::cout << "\n===\n" << std::endl;
  for (auto it : xml_model->param_objects) {
    std::cout << "Object id " << it.first << ": " << std::endl;
    DeformableObjectParameters::Pointer deformable = it.second;
    std::cout << "template file: " << deformable->GetFilename() << std::endl;
    deformable->PrintSelf(std::cout);
    std::cout << "\n===\n" << std::endl;
  }
#if ITK_VERSION_MAJOR >= 4
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

  typedef Diffeos<ScalarType, Dimension> DiffeosType;
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryListType;
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;
  typedef LinearInterpImage<ScalarType, Dimension> LIImageType;
  typedef ParametricImage<ScalarType, Dimension> ParametricImageType;
  typedef CrossSectionalDataSet<ScalarType, Dimension> DataSetType;
  typedef AbstractStatisticalModel<ScalarType, Dimension> AbstractStatisticalModelType;
  typedef AbstractAtlas<ScalarType, Dimension> AbstractAtlasType;
  typedef DeterministicAtlas<ScalarType, Dimension> DeterministicAtlasType;
  typedef BayesianAtlas<ScalarType, Dimension> BayesianAtlasType;
  typedef LdaAtlas<ScalarType, Dimension> LdaAtlasType;
  typedef BayesianAtlasMixture<ScalarType, Dimension> BayesianAtlasMixtureType;
  typedef AbstractEstimator<ScalarType, Dimension> AbstractEstimatorType;
  typedef GradientAscent<ScalarType, Dimension> GradientAscentType;
  typedef FastGradientAscent<ScalarType, Dimension> FastGradientAscentType;
  typedef McmcSaem<ScalarType, Dimension> McmcSaemType;
  typedef SrwMhwgSampler<ScalarType, Dimension> SrwMhwgSamplerType;
  typedef MalaSampler<ScalarType, Dimension> MalaSamplerType;
  typedef AmalaSampler<ScalarType, Dimension> AmalaSamplerType;
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef NormalDistribution<ScalarType> NormalDistributionType;
  xml_model->param_diffeos->Update();
  for (auto it : xml_model->param_objects)
    it.second->Update();
  /// Get the atlas type (either deterministic or bayesian).
  bool useDeterministicAtlas(0), useBayesianAtlas(1), useBayesianAtlasMixture(0), useLdaAtlas(0);
  if (numSubjects == 1) // Special case of the matching model.
  {
    std::cout << "The matching problem will be seen as a deterministic atlas estimation with fixed source object."
              << std::endl;
    useDeterministicAtlas = 1;
    useBayesianAtlas = 0;
    useBayesianAtlasMixture = 0;
    paramDiffeos->SetFreezeTemplate();
  } else {
    if (itksys::SystemTools::Strucmp(paramDiffeos->GetModelType().c_str(), "deterministicatlas") == 0) {
      useDeterministicAtlas = 1;
      useBayesianAtlas = 0;
    } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetModelType().c_str(), "mixture") == 0) {
      useBayesianAtlasMixture = 1;
      useBayesianAtlas = 0;
    } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetModelType().c_str(), "ldaatlas") == 0) {
      useLdaAtlas = 1;
      useBayesianAtlas = 0;
    } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetModelType().c_str(), "bayesianatlas") != 0) {
      std::cerr << "Unknown atlas model type: defaulting to bayesian." << std::endl;
    }
  }
  bool useGradientAscent(0), useFastGradientAscent(1), useSrwMhwgSaem(0), useMalaSaem(0), useAmalaSaem(0);
  if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientascent") == 0) {
    useGradientAscent = 1;
    useFastGradientAscent = 0;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "srwmhwgsaem") == 0) {
    if (useDeterministicAtlas)
      std::cerr << "It is not possible to estimate a deterministic atlas with an MCMC-SAEM algorithm."
          "Defaulting to fast gradient ascent." << std::endl;
    else if (useLdaAtlas)
      std::cerr << "It is not possible to estimate a lda atlas with an MCMC-SAEM algorithm."
          "Defaulting to fast gradient ascent." << std::endl;
    else
      useFastGradientAscent = 0;
    useSrwMhwgSaem = 1;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "malasaem") == 0) {
    if (useDeterministicAtlas)
      std::cerr << "It is not possible to estimate a deterministic atlas with an MCMC-SAEM algorithm."
          "Defaulting to fast gradient ascent." << std::endl;
    else if (useLdaAtlas){
      std::cerr << "It is not possible to estimate a bayesian atlas mixture with a MALA-SAEM algorithm."
          "Defaulting to fast gradient ascent" << std::endl;
    } else if (useBayesianAtlasMixture) {
      std::cerr << "It is not possible to estimate a bayesian atlas mixture with a MALA-SAEM algorithm."
          "Defaulting to SRW-MHWG-SAEM." << std::endl;
      useSrwMhwgSaem = 1;
      useFastGradientAscent = 0;
    } else
      useFastGradientAscent = 0;
    useMalaSaem = 1;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "amalasaem") == 0) {
    if (useDeterministicAtlas)
      std::cerr << "It is not possible to estimate a deterministic atlas with an MCMC-SAEM algorithm."
          "Defaulting to fast gradient ascent." << std::endl;
    else if (useBayesianAtlasMixture) {
      std::cerr << "It is not possible to estimate a bayesian atlas mixture with an AMALA-SAEM algorithm."
          "Defaulting to SRW-MHWG-SAEM." << std::endl;
      useSrwMhwgSaem = 1;
      useFastGradientAscent = 0;
    } else
      useFastGradientAscent = 0;
    useAmalaSaem = 1;
  } else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fastgradientascent")
      != 0) {
    std::cout << "Unspecified optimization method : defaulting to fast gradient ascent." << std::endl;
  } else if (useBayesianAtlasMixture) {
    std::cerr << "It is not possible to estimate a bayesian atlas mixture with a fast gradient ascent algorithm."
        "Defaulting to SRW-MHWG-SAEM." << std::endl;
    useSrwMhwgSaem = 1;
    useFastGradientAscent = 0;
  }
  /// Sets up the kernel factory - can't use any Update() before setting this up!
  KernelFactoryType *kfac = KernelFactoryType::Instantiate();
  kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
  kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());
  /**
   * Only one visit per subject is admited
   */
  if (std::any_of(xml_model->subjects.begin(),
                  xml_model->subjects.end(),
                  [](const def::io::XmlSubject &s) { return s.visits.size() > 1; }))
    throw std::runtime_error(std::string("Atlas require only one visit per subject"));
  /// Creates template and target objects.
  AbstractGeometryListType templateObjectList(numObjects);
  std::vector<AbstractGeometryListType> targetObjectList(numSubjects);
  for (int s = 0; s < numSubjects; s++) {
    AbstractGeometryListType Aux(numObjects);
    targetObjectList[s] = std::move(Aux);
  }
  /// Scans the list of objects.
  typedef DeformableObjectReader<ScalarType, Dimension> ObjectReaderType;
  std::size_t i_object = 0;
  for (auto it : xml_model->param_objects) {
    std::size_t i_subject = 0;
    for (auto subject_filename : xml_model->subjects_filename_by_object_id(it.first)) {
      ObjectReaderType reader;
      reader.SetObjectParameters(it.second);
      reader.SetFileName(subject_filename);
      reader.Update();
      targetObjectList[i_subject++][i_object] = reader.GetOutput();
    }
    ObjectReaderType reader;
    reader.SetObjectParameters(it.second);
    reader.SetTemplateType();
    auto fn = it.second->GetFilename();
    reader.SetFileName(fn);
    reader.Update();
    templateObjectList[i_object++] = reader.GetOutput();
  } // End scanning objects.
  /// Creates multi-objects for template aÒ/nd targets.
  typename std::vector<std::shared_ptr<DeformableMultiObjectType>> target(numSubjects);
  for (unsigned int s = 0; s < numSubjects; s++) {
    target[s] = std::make_shared<DeformableMultiObjectType>();
    target[s]->SetObjectList(targetObjectList[s]);
    target[s]->Update();
  }
  std::shared_ptr<DeformableMultiObjectType> templateObjects = std::make_shared<DeformableMultiObjectType>();
  templateObjects->SetObjectList(templateObjectList);
  templateObjects->Update();

  ///Get the class labels for all the objects, for the lda atlas modeL.
  ///What we actually do is to take the labels as a string and number them for later convenience.
  std::vector<unsigned int> classes(numSubjects);
  std::vector<std::string> classesFound;
  unsigned int nbClasses;

  if (useLdaAtlas) {
    size_t i_subject(0);
    for (auto subject_class_label : xml_model->subjects_class_label()){
      if (std::find(classesFound.begin(), classesFound.end(), subject_class_label) == classesFound.end())
        classesFound.push_back(subject_class_label);
      auto a = std::find(classesFound.begin(), classesFound.end(), subject_class_label) - classesFound.begin();
      classes[i_subject++] = a;
    }
    nbClasses = classesFound.size();
    std::cout << nbClasses << " classes identified in the data set." << std::endl;
    for (unsigned int i=0;i<classesFound.size();++i)
      std::cout << i << " corresponds to class " << classesFound[i] << std::endl;
  }

  /// Warnings.
  bool isOfBayesianType = useBayesianAtlas;
  /// Computes the noise dimension, for the bayesian atlas noise model.
  std::vector<unsigned long> noiseDimension(numObjects);
  if (isOfBayesianType or useLdaAtlas) {
    std::vector<unsigned long> nd_temp = templateObjects->GetDimensionOfDiscretizedObjects();
    std::vector<unsigned long> nd_targ = target[0]->GetDimensionOfDiscretizedObjects();
    for (unsigned long i = 0; i < numObjects; ++i) {
      if (templateObjectList[i]->IsOfLandmarkKind()) // Landmark kind => current & varifold metrics.
      {
        noiseDimension[i] = nd_temp[i];
      } else // Image kind (LinearInterpImage or ParametricImage) => ALL SUBJECTS' DATA ARE ASSUMED SAMPLED EQUALLY
      {
        noiseDimension[i] = nd_targ[i];
      }
    }
  }
  if ((numSubjects == 1) && (!paramDiffeos->FreezeTemplate())) {
    std::cout << "Warning: for what is obviously a matching problem, source object will be froozen" << std::endl;
    paramDiffeos->SetFreezeTemplate();
  } else if ((paramDiffeos->FreezeTemplate()) && isOfBayesianType) {
    std::cout << "Warning: template freeze is not implemented for Bayesian framework: template will be updated!"
              << std::endl;
    paramDiffeos->UnsetFreezeTemplate();
  }
  /// Creates the deformation object.
  std::shared_ptr<DiffeosType> def = std::make_shared<DiffeosType>();
  def->SetKernelWidth(paramDiffeos->GetKernelWidth());
  def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
  def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor()); // to define the bounding box
  ///Compact kernel usage :
  if (paramDiffeos->UseFastConvolutions()) {def->SetUseFastConvolutions();}
  else {def -> UnsetUseFastConvolutions();}
  if (not(paramDiffeos->UseImprovedEuler()))
    def->UseStandardEuler();
  if (paramDiffeos->ComputeTrueInverseFlow() == SparseDiffeoParameters::On) {
    std::cout
        << "Warning : an active compute-true-inverse-flow flag is usually not advised when used for the atlas model."
        << std::endl;
    def->SetComputeTrueInverseFlow();
  } else if (paramDiffeos->ComputeTrueInverseFlow() == SparseDiffeoParameters::Off) {
    def->UnsetComputeTrueInverseFlow();
  } else {
    def->UnsetComputeTrueInverseFlow();
  }
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
  /// Creates the dataSet object.
  DataSetType *dataSet = new DataSetType();
  dataSet->SetDeformableMultiObjects(target);
  if (useLdaAtlas)
    dataSet->SetClasses(classes);
  dataSet->Update();
  /// Creates the atlas and the estimator objects.
  std::shared_ptr<AbstractStatisticalModelType> model;
  AbstractEstimatorType *estimator;
  if (useBayesianAtlas) {
    std::shared_ptr<BayesianAtlasType> bayesianAtlasModel = std::make_shared<BayesianAtlasType>();
    bayesianAtlasModel->SetCovarianceMomenta_HyperParameter(
        numSubjects * paramDiffeos->GetCovarianceMomenta_Normalized_Hyperparameter());
    bayesianAtlasModel->SetCovarianceMomenta_Prior_Inverse(paramDiffeos->GetCovarianceMomenta_Prior_Inverse_fn());
    bayesianAtlasModel->SetNoiseDimension(noiseDimension);
    VectorType DataSigmaSquared_Hyperparameter(numObjects, 0.0);
    VectorType DataSigmaSquared_Prior(numObjects, 0.0);
    std::size_t i = 0;
    for (auto it : xml_model->param_objects) {
      DataSigmaSquared_Hyperparameter(i)
          = it.second->GetDataSigma_Normalized_Hyperparameter() * noiseDimension[i] * numSubjects;
      DataSigmaSquared_Prior(i) = it.second->GetDataSigma_Prior();
      ++i;
    }
    bayesianAtlasModel->SetDataSigmaSquared_HyperParameter(DataSigmaSquared_Hyperparameter);
    bayesianAtlasModel->SetDataSigmaSquared_Prior(DataSigmaSquared_Prior);
    if (useGradientAscent) {
      GradientAscentType *gradientAscentEstimator = new GradientAscentType();
      gradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
      gradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
      gradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
      gradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
      gradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
      estimator = static_cast<AbstractEstimatorType *>(gradientAscentEstimator);
    } else if (useFastGradientAscent) {
      FastGradientAscentType *fastGradientAscentEstimator = new FastGradientAscentType();
      fastGradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
      fastGradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
      fastGradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
      fastGradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
      fastGradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
      estimator = static_cast<AbstractEstimatorType *>(fastGradientAscentEstimator);
    } else if (useSrwMhwgSaem) {
      McmcSaemType *mcmcSaemEstimator = new McmcSaemType();
      std::shared_ptr<SrwMhwgSamplerType> sampler = std::make_shared<SrwMhwgSamplerType>();
      bayesianAtlasModel->SetUsePriorOnControlPoints();
      bayesianAtlasModel->SetUseRandomControlPoints();
      sampler->SetStatisticalModel(bayesianAtlasModel);
      sampler->SetDataSet(dataSet);
      mcmcSaemEstimator->SetSampler(sampler);
      estimator = static_cast<AbstractEstimatorType *>(mcmcSaemEstimator);
    } else if (useMalaSaem) {
      McmcSaemType *mcmcSaemEstimator = new McmcSaemType();
      std::shared_ptr<MalaSamplerType> sampler = std::make_shared<MalaSamplerType>();
      bayesianAtlasModel->SetUsePriorOnControlPoints();
      bayesianAtlasModel->SetUseRandomControlPoints();
      sampler->SetStatisticalModel(bayesianAtlasModel);
      sampler->SetDataSet(dataSet);
      mcmcSaemEstimator->SetSampler(sampler);
      estimator = static_cast<AbstractEstimatorType *>(mcmcSaemEstimator);
    } else if (useAmalaSaem) {
      McmcSaemType *mcmcSaemEstimator = new McmcSaemType();
      std::shared_ptr<AmalaSamplerType> sampler = std::make_shared<AmalaSamplerType>();
      bayesianAtlasModel->SetUsePriorOnControlPoints();
      bayesianAtlasModel->SetUseRandomControlPoints();
      sampler->SetStatisticalModel(bayesianAtlasModel);
      sampler->SetDataSet(dataSet);
      mcmcSaemEstimator->SetSampler(sampler);
      estimator = static_cast<AbstractEstimatorType *>(mcmcSaemEstimator);
    }
    model = std::static_pointer_cast<AbstractStatisticalModelType>(bayesianAtlasModel);
  } else if (useDeterministicAtlas) {
    std::shared_ptr<DeterministicAtlasType> deterministicAtlasModel = std::make_shared<DeterministicAtlasType>();
    std::string CMI_fn = paramDiffeos->GetCovarianceMomentaInverse_fn();
    if (CMI_fn.empty()) { deterministicAtlasModel->SetRKHSNormForRegularization(); }
    else { deterministicAtlasModel->SetCovarianceMomentaInverse(CMI_fn); }
    VectorType DataSigmaSquared(numObjects, 0.0);
    std::size_t i = 0;
    for (auto it : xml_model->param_objects) {
      ScalarType DSS = it.second->GetDataSigma();
      DataSigmaSquared(i++) = DSS * DSS;
    }
    deterministicAtlasModel->SetDataSigmaSquared(DataSigmaSquared);
    if (useGradientAscent) {
      GradientAscentType *gradientAscentEstimator = new GradientAscentType();
      gradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
      gradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
      gradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
      gradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
      gradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
      estimator = static_cast<AbstractEstimatorType *>(gradientAscentEstimator);
    } else if (useFastGradientAscent) {
      FastGradientAscentType *fastGradientAscentEstimator = new FastGradientAscentType();
      fastGradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
      fastGradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
      fastGradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
      fastGradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
      fastGradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
      estimator = static_cast<AbstractEstimatorType *>(fastGradientAscentEstimator);
    } else {
      std::cerr << "Error in the estimate_atlas.cxx launcher : "
          "a deterministic atlas can only be estimated by a gradient ascent algorithm." << std::endl;
    }
    model = std::static_pointer_cast<AbstractStatisticalModelType>(deterministicAtlasModel);
  } else if (useBayesianAtlasMixture) {
    VectorType DataSigmaSquared_Hyperparameter(numObjects, 0.0);
    VectorType DataSigmaSquared_Prior(numObjects, 0.0);
    std::size_t i = 0;
    for (auto it : xml_model->param_objects) {
      DataSigmaSquared_Hyperparameter(i)
          = it.second->GetDataSigma_Normalized_Hyperparameter() * noiseDimension[i] * numSubjects;
      DataSigmaSquared_Prior(i) = it.second->GetDataSigma_Prior();
      ++i;
    }
    std::shared_ptr<BayesianAtlasMixtureType> bayesianAtlasMixtureModel = std::make_shared<BayesianAtlasMixtureType>();
    const unsigned int nbClusters = 2;
    std::vector<std::shared_ptr<BayesianAtlasType>> atlases(nbClusters);
    for (unsigned int k = 0; k < nbClusters; ++k) {
      atlases[k] = std::make_shared<BayesianAtlasType>();
      atlases[k]->SetCovarianceMomenta_HyperParameter(
          numSubjects * paramDiffeos->GetCovarianceMomenta_Normalized_Hyperparameter());
      atlases[k]->SetCovarianceMomenta_Prior_Inverse(paramDiffeos->GetCovarianceMomenta_Prior_Inverse_fn());
      atlases[k]->SetNoiseDimension(noiseDimension);
      atlases[k]->SetDataSigmaSquared_HyperParameter(DataSigmaSquared_Hyperparameter);
      atlases[k]->SetDataSigmaSquared_Prior(DataSigmaSquared_Prior);
    }
    bayesianAtlasMixtureModel->SetAtlases(atlases);
    if (useSrwMhwgSaem) {
      McmcSaemType *mcmcSaemEstimator = new McmcSaemType();
      std::shared_ptr<SrwMhwgSamplerType> sampler = std::make_shared<SrwMhwgSamplerType>();
      for (unsigned int k = 0; k < nbClusters; ++k) {
        atlases[k]->SetUsePriorOnControlPoints();
        atlases[k]->SetUseRandomControlPoints();
      }
      sampler->SetStatisticalModel(model);
      sampler->SetDataSet(dataSet);
      mcmcSaemEstimator->SetSampler(sampler);
      estimator = static_cast<AbstractEstimatorType *>(mcmcSaemEstimator);
    } else {
      std::cerr << "Error in the estimate_atlas.cxx launcher : "
          "a bayesian atlas mixture can only be estimated by a SrwMhwgSaem algorithm (for now at least)." << std::endl;
    }
    model = std::static_pointer_cast<AbstractStatisticalModelType>(bayesianAtlasMixtureModel);
  } else if (useLdaAtlas){
    std::shared_ptr<LdaAtlasType> LdaAtlasModel = std::make_shared<LdaAtlasType>();
    LdaAtlasModel->SetNoiseDimension(noiseDimension);
    LdaAtlasModel->SetNumberOfClasses(nbClasses);
    LdaAtlasModel->SetIntraClassPCADimension(paramDiffeos->GetIntraClassPCADimension());

    if (useGradientAscent) {
      GradientAscentType *gradientAscentEstimator = new GradientAscentType();
      gradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
      gradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
      gradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
      gradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
      gradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
      estimator = static_cast<AbstractEstimatorType *>(gradientAscentEstimator);
    } else if (useFastGradientAscent) {
      FastGradientAscentType *fastGradientAscentEstimator = new FastGradientAscentType();
      fastGradientAscentEstimator->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
      fastGradientAscentEstimator->SetInitialStepSize(paramDiffeos->GetInitialStepSize());
      fastGradientAscentEstimator->SetAdaptiveShrink(paramDiffeos->GetStepShrink());
      fastGradientAscentEstimator->SetAdaptiveExpand(paramDiffeos->GetStepExpand());
      fastGradientAscentEstimator->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
      estimator = static_cast<AbstractEstimatorType *>(fastGradientAscentEstimator);
    }
    model = std::static_pointer_cast<AbstractStatisticalModelType>(LdaAtlasModel);
  }
  /////////////////////////
  /// Detailed options. ///
  /////////////////////////
  const unsigned int nbSubjects = dataSet->GetNumberOfSubjects();
  unsigned int nbControlPoints = 0, nbClusters = 0;
  MatrixType dataDomain;
  std::vector<std::shared_ptr<BayesianAtlasType>> atlases;
  if (model->IsDeterministicAtlas() || model->IsBayesianAtlas() || model->IsLdaAtlas()) {
    std::shared_ptr<AbstractAtlasType> aux = std::dynamic_pointer_cast<AbstractAtlasType>(model);
    aux->SetDiffeos(def);
    aux->SetTemplate(templateObjects);
    aux->SetFreezeTemplateFlag(paramDiffeos->FreezeTemplate());
    aux->SetFreezeControlPointsFlag(paramDiffeos->FreezeCP());
    aux->SetCPSpacing(paramDiffeos->GetInitialCPSpacing());
    std::string CP_fn = paramDiffeos->GetInitialCPPosition_fn();
    if (strlen(CP_fn.c_str())) { aux->SetControlPoints(CP_fn); }
    else if (paramDiffeos->OptimizeInitialControlPoints()) { aux->InitializeControlPoints(true); }
    aux->SetSmoothingKernelWidth(paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
    aux->SetNumberOfThreads(paramDiffeos->GetNumberOfThreads());
    model->Update();
    nbControlPoints = aux->GetControlPoints().rows();
    dataDomain = aux->GetBoundingBox();
  } else if (model->IsBayesianAtlasMixture()) {
    std::shared_ptr<BayesianAtlasMixtureType> aux = std::dynamic_pointer_cast<BayesianAtlasMixtureType>(model);
    atlases = aux->GetAtlases();
    nbClusters = atlases.size();
    // Quick and dirty differentiated initialization of the templates for a specific example (US postal, 2 atlases).
    std::shared_ptr<ParametricImageType> I0
        = std::static_pointer_cast<ParametricImageType>(templateObjects->GetObjectList()[0]->Clone());
    std::shared_ptr<ParametricImageType> I1
        = std::static_pointer_cast<ParametricImageType>(templateObjects->GetObjectList()[0]->Clone());
    I0->Initialize(std::static_pointer_cast<LIImageType>(dataSet->GetDataForSubject(0)->GetObjectList()[0]));
    I1->Initialize(std::static_pointer_cast<LIImageType>(dataSet->GetDataForSubject(2)->GetObjectList()[0]));
    std::shared_ptr<DeformableMultiObjectType> temp0 = std::make_shared<DeformableMultiObjectType>();
    std::shared_ptr<DeformableMultiObjectType> temp1 = std::make_shared<DeformableMultiObjectType>();
    std::vector<std::shared_ptr<AbstractGeometryType>> agl0(1);
    agl0[0] = I0;
    std::vector<std::shared_ptr<AbstractGeometryType>> agl1(1);
    agl1[0] = I1;
    temp0->SetObjectList(agl0);
    temp0->Update();
    temp1->SetObjectList(agl1);
    temp1->Update();
    atlases[0]->SetTemplate(temp0);
    atlases[1]->SetTemplate(temp1);
    for (unsigned int k = 0; k < nbClusters; ++k) {
//            atlases[k]->SetTemplate(templateObjects->Clone());
      atlases[k]->SetFreezeTemplateFlag(paramDiffeos->FreezeTemplate());
//            std::string CP_fn = paramDiffeos->GetInitialCPPosition_fn();
//            atlases[k]->SetControlPoints(CP_fn);
      atlases[k]->SetCPSpacing(paramDiffeos->GetInitialCPSpacing());
      atlases[k]->SetDiffeos(def->Clone());
      atlases[k]->SetSmoothingKernelWidth(
          paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
      atlases[k]->SetNumberOfThreads(paramDiffeos->GetNumberOfThreads());
    }
    model->Update();
    nbControlPoints = atlases[0]->GetControlPoints().rows();
    dataDomain = atlases[0]->GetBoundingBox();
  }
  /// Purely informative.
  std::cout << "Working domain (only template) : \torigin = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 0) << "] \tlength = [";
  for (unsigned int d = 0; d < Dimension - 1; ++d) { std::cout << dataDomain(d, 1) - dataDomain(d, 0) << " ; "; }
  std::cout << dataDomain(Dimension - 1, 1) - dataDomain(Dimension - 1, 0) << "]" << std::endl;
  for (unsigned int s = 0; s < numSubjects; s++) {
    MatrixType BB = target[s]->GetBoundingBox();
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
  /// Initialization of the momentas. Either a fixed or a random effect depending on the model model.
  const std::string fn = paramDiffeos->GetInitialMomenta_fn();
  std::vector<MatrixType> momentas;
  LinearVariablesMapType indRER;
  if (strlen(fn.c_str())) {
    momentas = readMultipleMatrixDLM<ScalarType>(fn.c_str());
    std::cout << "Using an initial momentas matrix list of size " << momentas.size()
              << " from file: " << fn << std::endl;
  } else {
    momentas.resize(nbSubjects);
    for (int s = 0; s < nbSubjects; s++) {
      momentas[s].set_size(nbControlPoints, Dimension);
      momentas[s].fill(0.0);
    }
  }
  if (not(useLdaAtlas))
    indRER["Momenta"] = momentas;
  estimator->InitializeIndividualRER(indRER);
  /// Optional initialization of the random control points realization, if a MCMC-SAEM algorithm is used.
  if (model->IsBayesianAtlas()) {
    std::shared_ptr<BayesianAtlasType> batlas = std::dynamic_pointer_cast<BayesianAtlasType>(model);
    if (batlas->UseRandomControlPoints()) {
      LinearVariableMapType popRER;
      popRER["ControlPoints"] = batlas->GetControlPoints();
      estimator->InitializePopulationRER(popRER);
    }
  } else if (model->IsBayesianAtlasMixture()) {
    if (atlases[0]->UseRandomControlPoints()) {
      LinearVariableMapType popRER;
      for (unsigned int k = 0; k < nbClusters; ++k)
        popRER["ControlPoints_Atlas" + std::to_string(k)] = atlases[k]->GetControlPoints();
      estimator->InitializePopulationRER(popRER);
    }
  }
  /// Optional initialization of the mixture parameters, in case of an atlas mixture model.
  if (model->IsBayesianAtlasMixture()) {
    std::vector<VectorType> mixtureCoefficients(nbSubjects);
    for (unsigned int i = 0; i < nbSubjects; ++i) {
      mixtureCoefficients[i].set_size(nbClusters);
      mixtureCoefficients[i].fill(1.0 / nbClusters);
    }
    indRER["MixtureCoefficients"] = mixtureCoefficients;
    estimator->InitializeIndividualRER(indRER);
  }
  ///Optional initialization of the variables of the lda atlas:
  if (model->IsLdaAtlas()){
    std::shared_ptr<LdaAtlasType> ldaAtlas = std::dynamic_pointer_cast<LdaAtlasType>(model);
    ldaAtlas->SetFreezeF(paramDiffeos->FreezeF());
    ldaAtlas->SetFreezeG(paramDiffeos->FreezeG());
    std::string G_fn = paramDiffeos->GetInitialG_fn();
    if (G_fn.length() > 0){
      MatrixType G = readMatrixDLM<ScalarType>(G_fn.c_str());
      ldaAtlas->SetG(G);
    }
    std::string F_fn = paramDiffeos->GetInitialF_fn();
    if (F_fn.length() > 0){
      MatrixType F = readMatrixDLM<ScalarType>(F_fn.c_str());
      ldaAtlas->SetF(F);
    }
    std::string alpha_fn = paramDiffeos->GetInitialAlpha_fn();
    if (alpha_fn.length() > 0){
      MatrixListType alpha = readMultipleMatrixDLM<ScalarType>(alpha_fn.c_str());
      ldaAtlas->SetAlpha(alpha);
    }
    std::string Beta_fn = paramDiffeos->GetInitialBeta_fn();
    if (Beta_fn.length() > 0){
      MatrixListType betaAux = readMultipleMatrixDLM<ScalarType>(Beta_fn.c_str());
      std::vector<MatrixType> betas(betaAux.size());
      for (unsigned int i=0; i<betas.size();++i)
        betas[i] = betaAux.at(i);
      indRER["beta"] = betas;
    }
    else{
      std::vector<MatrixType> betas(nbSubjects);
      std::shared_ptr<NormalDistributionType> betaDistribution = std::make_shared<NormalDistributionType>();
      MatrixType mean(paramDiffeos->GetIntraClassPCADimension(),1,0.);
      betaDistribution->SetMean(mean);
      MatrixType cov = diagonal_matrix(paramDiffeos->GetIntraClassPCADimension(), 0.001);
      betaDistribution->SetCovariance(cov);
      for (unsigned int i=0;i<nbSubjects;++i) {
        betas[i] = MatrixType(paramDiffeos->GetIntraClassPCADimension(), 1., 0.);
        betas[i].set_column(0, betaDistribution->Sample());
      }
      indRER["beta"] = betas;
    }
    estimator->InitializeIndividualRER(indRER);
  }
  /// Optional initialization of the sampler, in case of a Mcmc-Saem estimation algorithm.
  if (estimator->IsMcmcSaem()) {
    const unsigned int nmax = nbControlPoints * Dimension;
    unsigned int n = std::min(paramDiffeos->GetSrwProposalBlocksize(), nmax);
    if (n < 1) { n = nmax; }
    auto saem = dynamic_cast<McmcSaemType *>(estimator);
    auto sampler = saem->GetSampler();
    saem->SetMemoryWindowSize(paramDiffeos->GetPrintAcceptanceRatesWindow());
    saem->SetAdaptMemoryWindowSize(paramDiffeos->GetAdaptiveAcceptanceRatesWindow());
    const unsigned int arw = paramDiffeos->GetAdaptiveAcceptanceRatesWindow();
    if (arw < 1) { saem->UnsetAdaptiveMode(); }
    else {
      saem->SetAdaptiveMode();
      saem->SetAdaptMemoryWindowSize(arw);
      sampler->SetAcceptanceRatesTarget(paramDiffeos->GetAdaptiveAcceptanceRatesTarget());
    }
    if (sampler->IsSrwMhwg()) {
      std::shared_ptr<SrwMhwgSamplerType> sampler
          = std::dynamic_pointer_cast<SrwMhwgSamplerType>(saem->GetSampler());
      MultiScalarNormalDistributionMapType popPD, indPD; // PD : "Proposal Distribution"
      // Std normal distributions. Block size = 1.
      auto cp = std::make_shared<MultiScalarNormalDistributionType>();
      auto mom = std::make_shared<MultiScalarNormalDistributionType>();
      // Block normal distribution, blockSize = random variable size.
      cp->SetMean(VectorType(n, 0.0));
      cp->SetVarianceSqrt(paramDiffeos->GetSrwControlPointsProposalStd());
      mom->SetMean(VectorType(n, 0.0));
      mom->SetVarianceSqrt(paramDiffeos->GetSrwMomentaProposalStd());
      if (model->IsBayesianAtlas()) {
        popPD["ControlPoints"] = cp;
        indPD["Momenta"] = mom;
      } else if (model->IsBayesianAtlasMixture()) {
        for (unsigned int k = 0; k < nbClusters; ++k) {
          popPD["ControlPoints_Atlas" + std::to_string(k)] = cp;
          indPD["Momenta_Atlas" + std::to_string(k)] = mom;
        }
        /// TODO the proposal distributions system has to be adapted here.
//        auto mix = std::make_shared<ConditionedNormalDistributionType>();
//        mix->SetMean(VectorType(nbClusters, 0.0));
//        mix->SetCovarianceSqrt(diagonal_matrix<ScalarType>(nbClusters,
//                                                           paramDiffeos->GetSrwMixtureCoefficientsProposalStd()));
//        indPD["MixtureCoefficients"] = mix;
      }
      sampler->SetPopulationProposalDistribution(popPD);
      sampler->SetIndividualProposalDistribution(indPD);
    } else if (sampler->IsMala()) {
      std::shared_ptr<MalaSamplerType> sampler = std::dynamic_pointer_cast<MalaSamplerType>(saem->GetSampler());
      std::map<std::string, ScalarType> scales;
      if (model->IsBayesianAtlas()) {
        scales["ControlPoints"] = paramDiffeos->GetMalaControlPointsScale();
        scales["Momenta"] = paramDiffeos->GetMalaMomentaScale();
      } else if (model->IsBayesianAtlasMixture()) {
        for (unsigned int k = 0; k < nbClusters; ++k) {
          scales["ControlPoints_Atlas" + std::to_string(k)] = paramDiffeos->GetMalaControlPointsScale();
          scales["Momenta_Atlas" + std::to_string(k)] = paramDiffeos->GetMalaMomentaScale();
        }
      }
      sampler->Sethreshold(1000);
      sampler->SetScales(scales);
    } else if (sampler->IsAmala()) {
      std::shared_ptr<AmalaSamplerType> sampler
          = std::dynamic_pointer_cast<AmalaSamplerType>(saem->GetSampler());
      std::map<std::string, ScalarType> meanScales, covarianceScales, regularizations;
      if (model->IsBayesianAtlas()) {
        meanScales["ControlPoints"] = paramDiffeos->GetAmalaControlPointsMeanScale();
        meanScales["Momenta"] = paramDiffeos->GetAmalaMomentaMeanScale();
        covarianceScales["ControlPoints"] = paramDiffeos->GetAmalaControlPointsCovarianceScale();
        covarianceScales["Momenta"] = paramDiffeos->GetAmalaMomentaCovarianceScale();
        regularizations["ControlPoints"] = paramDiffeos->GetAmalaControlPointsRegularization();
        regularizations["Momenta"] = paramDiffeos->GetAmalaMomentaRegularization();
      } else if (model->IsBayesianAtlasMixture()) {
        for (unsigned int k = 0; k < nbClusters; ++k) {
          meanScales["ControlPoints_Atlas" + std::to_string(k)] = paramDiffeos->GetAmalaControlPointsMeanScale();
          meanScales["Momenta_Atlas" + std::to_string(k)] = paramDiffeos->GetAmalaMomentaMeanScale();
          covarianceScales["ControlPoints_Atlas" + std::to_string(k)] =
              paramDiffeos->GetAmalaControlPointsCovarianceScale();
          covarianceScales["Momenta_Atlas" + std::to_string(k)] = paramDiffeos->GetAmalaMomentaCovarianceScale();
          regularizations["ControlPoints_Atlas" + std::to_string(k)] =
              paramDiffeos->GetAmalaControlPointsRegularization();
          regularizations["Momenta_Atlas" + std::to_string(k)] = paramDiffeos->GetAmalaMomentaRegularization();
        }
      }
      sampler->Sethreshold(1000);
      sampler->SetMeanScales(meanScales);
      sampler->SetCovarianceScales(covarianceScales);
      sampler->SetRegularizations(regularizations);
    }
  }
  estimator->SetStatisticalModel(model);
  estimator->SetDataSet(dataSet);
  estimator->SetMaxIterations(paramDiffeos->GetMaxIterations());
  estimator->SetPrintEveryNIters(paramDiffeos->GetPrintEveryNIters());
  estimator->SetSaveEveryNIters(paramDiffeos->GetSaveEveryNIters());
  // estimator->SetNumberOfThreads(paramDiffeos->GetNumberOfThreads()); // TODO
  /// Creates output names for saving.
  // All saving name will start with AtlasName, so this is the place to add the path to a particular folder!
  std::vector<std::string> objectsName(numObjects);
  std::vector<std::string> objectsNameExtension(numObjects);
  auto template_iterator = xml_model->param_objects.begin();
  for (unsigned int i = 0; i < numObjects; ++i, ++template_iterator) {
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
  }
  std::string modelName = paramDiffeos->GetModelName();
  if (model->IsDeterministicAtlas()) {
    if (modelName.empty()) {
      if (numSubjects == 1) { model->SetName("Registration"); }
      else { model->SetName("Atlas"); }
    }
    std::shared_ptr<AbstractAtlasType> aux = std::dynamic_pointer_cast<AbstractAtlasType>(model);
    aux->SetTemplateObjectsName(objectsName);
    aux->SetTemplateObjectsNameExtension(objectsNameExtension);
  } else if (model->IsBayesianAtlas()) {
    if (modelName.empty()) { model->SetName("Atlas"); }
    std::shared_ptr<AbstractAtlasType> aux = std::dynamic_pointer_cast<AbstractAtlasType>(model);
    aux->SetTemplateObjectsName(objectsName);
    aux->SetTemplateObjectsNameExtension(objectsNameExtension);
  } else if (model->IsBayesianAtlasMixture()) {
    if (modelName.empty()) { model->SetName("BayesianAtlasMixture"); }
    for (unsigned int k = 0; k < nbClusters; ++k) {
      atlases[k]->SetName("Atlas" + std::to_string(k));
      atlases[k]->SetTemplateObjectsName(objectsName);
      atlases[k]->SetTemplateObjectsNameExtension(objectsNameExtension);
    }
  } else if (model->IsLdaAtlas()) {
    if (modelName.empty())
      model->SetName("LdaAtlas");
    std::shared_ptr<AbstractAtlasType> aux = std::dynamic_pointer_cast<AbstractAtlasType>(model);
    aux->SetTemplateObjectsName(objectsName);
    aux->SetTemplateObjectsNameExtension(objectsNameExtension);
  }
  // Creates timer.
  SimpleTimer timer;


  if (paramDiffeos->OnlyWrite()) {
    std::cout << "Will not estimate anything, but only write the atlas output." << std::endl;
    model->Update();
    LinearVariableMapType popRER;
    model->Write(dataSet, popRER, indRER);
  }
  else {
    estimator->Update();
    timer.Stop();
    std::cout << "Atlas estimation took "
              << timer.GetElapsedHours() << " hours, "
              << timer.GetElapsedMinutes() << " minutes, "
              << timer.GetElapsedSeconds() << " seconds"
              << std::endl;
    delete dataSet;
    delete estimator;
    KernelFactoryType::Delete();
  }
}
template void estimateAtlas<2>(std::shared_ptr<const def::io::XmlModel> xml_model);
template void estimateAtlas<3>(std::shared_ptr<const def::io::XmlModel> xml_model);