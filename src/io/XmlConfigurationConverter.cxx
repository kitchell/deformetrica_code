/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "XmlConfigurationConverter.h"

namespace def {
namespace io {

void XmlConfigurationConverter::load_xml_file(itk::DOMNodeXMLReader::Pointer& reader, const std::string& filename)  throw(std::runtime_error) {
  try {
    reader = itk::DOMNodeXMLReader::New();
    reader->SetFileName(filename);
    reader->Update();
  } catch (...) {
    throw std::runtime_error(std::string("Error in loading the " + filename + " file. "
        "Please check the given path is correct."));
  }
}

void XmlConfigurationConverter::load_xml_model(const std::string filename) throw(std::runtime_error) {
  load_xml_file(itk_xml_model_, filename);
}

void XmlConfigurationConverter::load_xml_data_set(const std::string& filename) throw(std::runtime_error) {
  load_xml_file(itk_xml_data_set_, filename);
}

void XmlConfigurationConverter::load_xml_optimization_parameters(const std::string& filename) throw(std::runtime_error) {
  load_xml_file(itk_xml_optimization_parameters_, filename);
  has_optimication_parameters_ = true;
}

void XmlConfigurationConverter::validate_configuration() throw(std::runtime_error) {
  XmlDictionary model;
  param_diffeos_ = SparseDiffeoParameters::New();

  programming_dynamic_xml_model_templates(model);

  programming_xml_model(model, param_diffeos_);
  XmlDictionary::xml_reader(itk_xml_model_->GetOutput(), model, "model");

  if (!has_optimication_parameters_) return;

  XmlDictionary optimization_parameters;
  programming_xml_optimization_parameters(optimization_parameters, param_diffeos_);
  XmlDictionary::xml_reader(itk_xml_optimization_parameters_->GetOutput(), optimization_parameters, "optimization-parameters");
}

std::shared_ptr<XmlModel> XmlConfigurationConverter::generate_xml_model() throw(std::runtime_error) {
  XmlDictionary model;

  programming_dynamic_xml_model_object(model);

  programming_xml_model(model, xml_model_->param_diffeos);
  XmlDictionary::xml_reader(itk_xml_model_->GetOutput(), model, "model");

  if (has_optimication_parameters_) {
    XmlDictionary optimization_parameters;
    programming_xml_optimization_parameters(optimization_parameters, xml_model_->param_diffeos);
    XmlDictionary::xml_reader(itk_xml_optimization_parameters_->GetOutput(), optimization_parameters, "optimization-parameters");
  }

  return xml_model_;
}

void XmlConfigurationConverter::programming_xml_model(XmlDictionary& xml, SparseDiffeoParameters::Pointer sp) {
  xml["model-type"]
      .required(true)
      .one_of<std::string>("Shooting", "Registration", "DeterministicAtlas", "BayesianAtlas", "LdaAtlas", "LongitudinalAtlas", "Regression", "ParallelTransport")
      .assign_to<std::string>(sp, &SparseDiffeoParameters::SetModelType);
  xml["model-name"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetModelName);

  xml["initial-cp-spacing"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetInitialCPSpacing);

  xml["initial-cp-position"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialCPPosition_fn);
  xml["initial-mom-values"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialMomenta_fn);
  xml["modmat"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetModulationMatrix_fn);

  xml["reftime"].assign_to<double>(sp, &SparseDiffeoParameters::SetReferenceTime);
  xml["timeshift-var"].assign_to<double>(sp, &SparseDiffeoParameters::SetTimeShiftVariance);
  xml["logacc-var"].assign_to<double>(sp, &SparseDiffeoParameters::SetLogAccelerationVariance);

  xml["initial-cp-to-transport"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialCPPositionForTransport_fn);
  xml["initial-mom-values-to-transport"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialMomentaForTransport_fn);

  xml["matching-number-of-timepoints"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetMatchingNumberOfTimePoints);

  xml["matching-t0"].assign_to<double>(sp, &SparseDiffeoParameters::SetMatchingT0);
  xml["matching-tn"].assign_to<double>(sp, &SparseDiffeoParameters::SetMatchingTN);
  xml["use-exp-parallelization"].filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetUseExpParallelization() : sp->UnsetUseExpParallelization();
      });

  xml["nb-sources"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetNumberOfSources);

  xml["initial-logacc-std"].assign_to<double>(sp, &SparseDiffeoParameters::SetLogAccelerationRandomEffectStd);

  xml["prior-reftime-mean"].assign_to<double>(sp, &SparseDiffeoParameters::SetReferenceTimePriorMean);
  xml["prior-reftime-std"].assign_to<double>(sp, &SparseDiffeoParameters::SetReferenceTimePriorStd);

  xml["deformation-parameters"]["kernel-width"].assign_to<double>(sp, &SparseDiffeoParameters::SetKernelWidth);

  xml["deformation-parameters"]["kernel-type"]
      .one_of<std::string>("EXACT","CUDAEXACT","P3M","COMPACT")
      .assign_to<std::string>(sp, &SparseDiffeoParameters::SetKernelType);

  xml["deformation-parameters"]["t0"]
      .assign_to<double>(sp, &SparseDiffeoParameters::SetT0);
  xml["deformation-parameters"]["tn"]
      .assign_to<double>(sp, &SparseDiffeoParameters::SetTN);

  xml["deformation-parameters"]["number-of-timepoints"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetNumberOfTimePoints);

  ///For abc sampling
  xml["deformation-parameters"]["number-of-cps"]
      .assign_to<unsigned int>(sp,&SparseDiffeoParameters::SetNumberOfCps);

  xml["number-of-samples"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetNumberOfSamples);

  xml["momenta-proposition-variance"]
      .assign_to<double>(sp, &SparseDiffeoParameters::SetMomentaPropositionVariance);

  xml["training-set-size"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetTrainingSetSize);

  xml["deformation-parameters"]["concentration-of-timepoints"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetConcentrationOfTimePointsForReferenceGeodesic);
  xml["deformation-parameters"]["number-of-timepoints-exp"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetNumberOfTimePointsForExponentiation);
  xml["deformation-parameters"]["margin-on-geodesic"]
      .assign_to<double>(sp, &SparseDiffeoParameters::SetMarginOnGeodesicLength);

  ///For LDA model:
  xml["intra-class-pca-dimension"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetIntraClassPCADimension);

  xml["initial-G"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialG_fn);

  xml["initial-beta"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialBeta_fn);

  xml["initial-alpha"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialAlpha_fn);

  xml["initial-F"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetInitialF_fn);

  xml["only-write"].filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetOnlyWrite() : sp->UnsetOnlyWrite();
      });


}

void XmlConfigurationConverter::programming_xml_optimization_parameters(XmlDictionary& xml, SparseDiffeoParameters::Pointer sp) {
  xml["optimization-method-type"]
      .one_of<std::string>("GRADIENTASCENT", "FASTGRADIENTASCENT", "SRWMHWGSAEM", "AMALASAEM", "POWELLSMETHOD")
      .assign_to<std::string>(sp, &SparseDiffeoParameters::SetOptimizationMethodType);

  xml["initial-step-size"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetInitialStepSize);

  xml["p3m-padding-factor"].assign_to<double>(sp, &SparseDiffeoParameters::SetP3MPaddingFactor);
  xml["p3m-working-spacing-ratio"].assign_to<double>(sp, &SparseDiffeoParameters::SetP3MWorkingSpacingRatio);
  xml["max-iterations"]
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetMaxIterations);

  xml["max-line-search-iterations"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetMaxLineSearchIterations);
  xml["step-expand"].assign_to<double>(sp, &SparseDiffeoParameters::SetStepExpand);
  xml["step-shrink"].assign_to<double>(sp, &SparseDiffeoParameters::SetStepShrink);
  xml["adaptive-tolerance"].assign_to<double>(sp, &SparseDiffeoParameters::SetAdaptiveTolerance);
  xml["print-every-n-iters"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetPrintEveryNIters);
  xml["save-every-n-iters"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetSaveEveryNIters);

  xml["print-ar-window"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetPrintAcceptanceRatesWindow);
  xml["adaptive-ar-window"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetAdaptiveAcceptanceRatesWindow);
  xml["adaptive-ar-target"].assign_to<double>(sp, &SparseDiffeoParameters::SetAdaptiveAcceptanceRatesTarget);

  xml["srw-proposal-blocksize"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetSrwProposalBlocksize);

  xml["srw-proposal-std-template"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwTemplateDataProposalStd);
  xml["srw-proposal-kernel-template"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwTemplateDataProposalKernelWidth);
  xml["srw-proposal-ncp-template"]
      .range<unsigned int>(def::io::range::value::positive_exclude_zero)
      .assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetSrwTemplateDataProposalNumberOfControlPoints);


  xml["srw-proposal-std-cp"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwControlPointsProposalStd);
  xml["srw-proposal-std-mom"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwMomentaProposalStd);
  xml["srw-proposal-std-modmat"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwModulationMatrixProposalStd);
  xml["srw-proposal-std-sources"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwSourcesProposalStd);
  xml["srw-proposal-std-logacc"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwLogAccelerationProposalStd);
  xml["srw-proposal-std-tshift"]
      .range<double>(def::io::range::value::positive_exclude_zero)
      .assign_to<double>(sp, &SparseDiffeoParameters::SetSrwTimeShiftProposalStd);

// /  xml["amala-scale"].assign_to<double>(sp, &SparseDiffeoParameters::SetAmalaScale);
  xml["number-of-threads"].assign_to<unsigned int>(sp, &SparseDiffeoParameters::SetNumberOfThreads);
  xml["covariance-momenta-normalized-hyperparameter"].assign_to<double>(sp, &SparseDiffeoParameters::SetCovarianceMomenta_Normalized_Hyperparameter);
  xml["covariance-momenta-prior-inverse-fn"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetCovarianceMomenta_Prior_Inverse_fn);
  xml["covariance-momenta-inverse-fn"].assign_to<std::string>(sp, &SparseDiffeoParameters::SetCovarianceMomentaInverse_fn);
  xml["maximum-number-of-class"].assign_to<int>(sp, &SparseDiffeoParameters::SetMaximumNumberOfClasses);
  xml["smoothing-kernel-with-ratio"].assign_to<double>(sp, &SparseDiffeoParameters::SetSmoothingKernelWidthRatio);

  xml["use-tempering"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetUseTempering() : sp->UnsetUseTempering();
      });
  xml["initial-temperature"].assign_to<double>(sp, &SparseDiffeoParameters::SetInitialTemperature);
  xml["tempering-duration-ratio"].assign_to<double>(sp, &SparseDiffeoParameters::SetTemperingDurationRatio);

  xml["freeze-cp"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetFreezeCP() : sp->UnsetFreezeCP();
      });

  xml["freeze-template"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetFreezeTemplate() : sp->UnsetFreezeTemplate();
      });

  xml["compute-true-inverse-flow"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetComputeTrueInverseFlow() : sp->UnsetComputeTrueInverseFlow();
      });

  xml["use-implicit-euler"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetUseImplicitEuler() : sp->UnsetUseImplicitEuler();
      });

  xml["use-improved-euler"]
        .filter_with(def::io::filters::lower_case())
        .add_reader([sp](const std::string &v) {
            v == "on" ? sp->SetUseImprovedEuler() : sp->UnsetUseImprovedEuler();
        });

  xml["optimize-initial-cp"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
          v == "on" ? sp->SetOptimizeInitialControlPoints() : sp->UnsetOptimizeInitialControlPoints();
      });

  xml["use-fast-convolutions"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
          v == "on" ? sp->SetUseFastConvolutions() : sp->UnsetUseFastConvolutions();
      });

  xml["multivariate-line-search"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
          v == "on" ? sp->SetMultivariateLineSearch() : sp->UnsetMultivariateLineSearch();
      });

  xml["write-full-trajectories"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
          v == "on" ? sp->SetWriteFullTrajectories() : sp->UnsetWriteFullTrajectories();
      });

    xml["use-forward-shooting"]
        .filter_with(def::io::filters::lower_case())
        .add_reader([sp](const std::string &v) {
      v == "on" ? sp->SetUseForwardShooting() : sp->UnsetUseForwardShooting();
    });

  /// For Lda atlas.
  xml["freeze-F"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetFreezeF() : sp->UnsetFreezeF();
      });

  xml["freeze-G"]
      .filter_with(def::io::filters::lower_case())
      .add_reader([sp](const std::string &v) {
        v == "on" ? sp->SetFreezeG() : sp->UnsetFreezeG();
      });

}

XmlSubjects XmlConfigurationConverter::data_set_subjects() {
  XmlSubjects xml_subjects;

  itk::DOMNode::Pointer dom = itk_xml_data_set_->GetOutput();

  if (itksys::SystemTools::LowerCase(dom->GetName()) != "data-set")
    throw std::runtime_error(std::string("Missing <data-set> initial xml tag instead of <" + dom->GetName() + ">"));

  itk::DOMNode::ChildrenListType subjects;
  dom->GetChildren("subject", subjects);

  if (subjects.empty())
    throw std::runtime_error("Missing <subject> xml tag");

  for (auto it_subject : subjects) {
    XmlSubject xml_subject;

    /* get subject id */
    {
      auto subject_attribute = it_subject->GetAttribute("id");
      if (subject_attribute.empty())
        throw std::runtime_error("Missing 'id' xml attribute for subject");

      xml_subject.id = std::move(subject_attribute);
    }


    /* get subject  Class label*/
    {
      auto class_label = it_subject->GetAttribute("label");
      xml_subject.class_label = std::move(class_label);
    }
    /* get if subject in test or train set*/
    {
      ///test="1" for test, anything else for train.
      auto test = it_subject->GetAttribute("is-in-test-set");
      xml_subject.test = std::move(test);
    }


    itk::DOMNode::ChildrenListType visits;
    it_subject->GetChildren("visit", visits);

    if (visits.empty())
      throw std::runtime_error("Missing <visit> xml tag");

    XmlVisits xml_visits;
    for (auto it_visit : visits) {
      XmlVisit xml_visit;


      /* ID */
      {
        auto visit_id = it_visit->GetAttribute("id");
        if (visit_id.empty())
          throw std::runtime_error("Missing 'id' xml attribute for visit");
        xml_visit.id = std::move(visit_id);
      }


      /* AGE */
      {
        auto age = it_visit->GetChild("age");
        if ((xml_visit.age_assigned = (age && age->GetTextChild()))) {

          std::stringstream ss;
          std::string text = age->GetTextChild()->GetText();
          if (text.empty())
            throw std::runtime_error("Missing value for <age> xml tag");

          ss << text;
          ss >> xml_visit.age;
        }
      }

      /* FILENAME */
      {
        itk::DOMNode::ChildrenListType filenames;
        it_visit->GetChildren("filename", filenames);

        if (filenames.empty())
          throw std::runtime_error("Missing <filename> xml tag");

        for (auto it_filename : filenames) {
          XmlObject xml_object;
          auto id = it_filename->GetAttribute("object_id");
          if (id.empty())
            throw std::runtime_error("Missing 'object_id' xml tag on configuration file");

          xml_object.id = std::move(id);

          if(!it_filename->GetTextChild())
            throw std::runtime_error("Missing data for <filename> xml tag on configuration file");

          auto text = it_filename->GetTextChild()->GetText();
          if (text.empty())
            throw std::runtime_error("Missing value for <filename> xml tag on configuration file");

          xml_object.subject_filename = std::move(text);
          xml_visit.objects.push_back(std::move(xml_object));
        }
      }

      xml_visits.push_back(std::move(xml_visit));
    }

    xml_subject.visits = std::move(xml_visits);
    xml_subjects.push_back(std::move(xml_subject));
  }

  return xml_subjects;
}

void XmlConfigurationConverter::programming_dynamic_xml_model_object(XmlDictionary& templates) {

  xml_model_->subjects = data_set_subjects();
  if (!xml_model_->check_object_id_consistency())
    throw std::runtime_error(std::string("The xml structure of data_set.xml must be the same for each visit and subject : equal number of object_id"));

  for (auto tag : xml_model_->extract_object_id()) {
    tag = def::support::utilities::strtolower(tag);
    auto& object = templates["template"]["object"].add_tag("id", tag);

    DeformableObjectParameters::Pointer def = DeformableObjectParameters::New();
    xml_model_->param_objects[tag] = def;

    object["deformable-object-type"].assign_to<std::string>(def,&DeformableObjectParameters::SetDeformableObjectType);
    object["data-sigma"].assign_to<double>(def,&DeformableObjectParameters::SetDataSigma);
    object["kernel-type"].assign_to<std::string>(def,&DeformableObjectParameters::SetKernelType);
    object["photometric-cp-spacing"].assign_to<double>(def,&DeformableObjectParameters::SetPhotometricCPSpacing);
    object["kernel-width"].assign_to<double>(def,&DeformableObjectParameters::SetKernelWidth);
    object["image-grid-downsampling"].assign_to<double>(def,&DeformableObjectParameters::SetImageGridDownsampling);
    object["data-sigma-normalized-hyperparameter"].assign_to<double>(def,&DeformableObjectParameters::SetDataSigma_Normalized_Hyperparameter);
    object["data-sigma-prior"].assign_to<double>(def,&DeformableObjectParameters::SetDataSigma_Prior);
    object["anatomical-coordinate-system"].assign_to<std::string>(def,&DeformableObjectParameters::SetAnatomicalCoordinateSystem);
    object["reorient-normals"]
        .filter_with(def::io::filters::lower_case())
        .add_reader([def](const std::string &v) {
          v == "on" ? def->SetReOrient() : def->UnsetReOrient();
        });
    object["filename"]
        .required(true)
        .add_reader([this, tag](const std::string &v) {
          xml_model_->param_objects[tag]->SetFilename(v);
        });
  }
}

void XmlConfigurationConverter::programming_dynamic_xml_model_templates(XmlDictionary& templates) {
  std::map<std::string, std::vector<std::string> > files;
  std::map<std::string, std::vector<double> > ages;

  itk::DOMNode::Pointer dom = itk_xml_data_set_->GetOutput();

  if (itksys::SystemTools::LowerCase(dom->GetName()) != "data-set")
    throw std::runtime_error(std::string("Missing <data-set> initial xml tag instead of <" + dom->GetName() + ">"));

  itk::DOMNode::ChildrenListType subjects;
  dom->GetChildren("subject", subjects);

  if (subjects.empty())
    throw std::runtime_error("Missing <subject> xml tag");

//    ASSERT_EQ(subjects.size(), 2);
  for (auto it_subject : subjects) {
    itk::DOMNode::ChildrenListType visits;
    it_subject->GetChildren("visit", visits);

    if (visits.empty())
      throw std::runtime_error("Missing <visit> xml tag");

//      ASSERT_LE(visits.size(), 2);

    for (auto it_visit : visits) {
      /* AGE */
      auto age = it_visit->GetChild("age");
      double age_value = 0;
      if (age && age->GetTextChild()) {
        std::stringstream ss;
        std::string text = age->GetTextChild()->GetText();
        if (text.empty())
          throw std::runtime_error("Missing value for <age> xml tag");
        ss << text;
        ss >> age_value;
      }

      itk::DOMNode::ChildrenListType filenames;
      it_visit->GetChildren("filename", filenames);

      if (filenames.empty())
        throw std::runtime_error("Missing <filename> xml tag");

//        ASSERT_EQ(filenames.size(), 2);

      for (auto it_filename : filenames) {
        auto id = it_filename->GetAttribute("object_id");
        if (id.empty())
          throw std::runtime_error("Missing 'object_id' xml tag on configuration file");

        if(!it_filename->GetTextChild())
          throw std::runtime_error("Missing data for <filename> xml tag on configuration file");

        auto text = it_filename->GetTextChild()->GetText();
        if (text.empty())
          throw std::runtime_error("Missing value for <filename> xml tag on configuration file");

        files[id].push_back(text);
        ages[id].push_back(age_value);
      }
      //TODO: adding more timepoint visit for the future..
      //TODO: check xml structure consistencies
      break;
    }
  }

  /* BUILDING DATA */
  num_subjects_ = subjects.size();
  num_objects_ = files.size();

//    ASSERT_EQ(num_subjects, 2);
//    ASSERT_EQ(num_objects, 2);

  int observation_order = 0;
  objects_order_ = 0;
  templatefn_.resize(num_objects_);
  subjectfn_.resize(num_objects_);
  observations_.resize(num_objects_);

  for (auto it : ages) {
    for (auto ik : it.second)
      observations_[observation_order].push_back(ik);
    ++observation_order;
  }

  for (auto it : files) {
    auto tag = it.first;

    auto file_list = files[tag];
    auto& object = templates["template"]["object"].add_tag("id", tag);

    DeformableObjectParameters::Pointer def = DeformableObjectParameters::New();

    object["deformable-object-type"].assign_to<std::string>(def,&DeformableObjectParameters::SetDeformableObjectType);
    object["data-sigma"].assign_to<double>(def,&DeformableObjectParameters::SetDataSigma);
    object["kernel-type"].assign_to<std::string>(def,&DeformableObjectParameters::SetKernelType);
    object["photometric-cp-spacing"].assign_to<double>(def,&DeformableObjectParameters::SetPhotometricCPSpacing);
    object["kernel-width"].assign_to<double>(def,&DeformableObjectParameters::SetKernelWidth);
    object["image-grid-downsampling"].assign_to<double>(def,&DeformableObjectParameters::SetImageGridDownsampling);
    object["data-sigma-normalized-hyperparameter"].assign_to<double>(def,&DeformableObjectParameters::SetDataSigma_Normalized_Hyperparameter);
    object["data-sigma-prior"].assign_to<double>(def,&DeformableObjectParameters::SetDataSigma_Prior);
    object["anatomical-coordinate-system"].assign_to<std::string>(def,&DeformableObjectParameters::SetAnatomicalCoordinateSystem);
    object["reorient-normals"]
        .filter_with(def::io::filters::lower_case())
        .add_reader([def](const std::string &v) {
          v == "on" ? def->SetReOrient() : def->UnsetReOrient();
        });

    param_objects_.push_back(std::move(def));

    object["filename"]
        .required(true)
        .add_reader([this, file_list](const std::string &v) {
          this->templatefn_[this->objects_order_] = v;
          for (auto file : file_list)
            this->subjectfn_[this->objects_order_].push_back(file);
          ++this->objects_order_;
        });

  }

}

}
}
