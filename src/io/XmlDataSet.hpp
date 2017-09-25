/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

#include <src/io/SparseDiffeoParameters.h>
#include <src/io/DeformableObjectParameters.h>
#include <src/support/utilities/Utils.hpp>
#include <memory>
#include <map>

namespace def {
namespace io {

/**
 * The following structs represent in a simple way the xml of test/data/regression_multiple_objects_and_subjects/data_set.xml file
 */

struct XmlModel;
struct XmlSubject;
struct XmlVisit;
struct XmlObject;

typedef std::vector<XmlSubject> XmlSubjects;
typedef std::vector<XmlVisit> XmlVisits;
typedef std::vector<XmlObject> XmlObjects;

struct XmlSubject {
  std::string id;
  XmlVisits visits;
  std::string class_label;
  std::string test;
};

struct XmlVisit {
  std::string id;
  ///Optional attribute, for supervised learning algos.
  std::string class_label;
  double age;
  bool age_assigned;
  XmlObjects objects;
};

struct XmlObject {
  std::string id;
  std::string subject_filename;
};

/**
 * The inputs files 'model.xml', 'data_set.xml' and 'optimization_parameters.xml'
 * are here described in the related object representation
 */
struct XmlModel {

  XmlModel() : param_diffeos(SparseDiffeoParameters::New()) {}

  /* The object description of model.xml file */
  SparseDiffeoParameters::Pointer param_diffeos;

  /* Mapping each <object id="x"> tag of model.xml file with the corresponding object representation */
  std::map<std::string, DeformableObjectParameters::Pointer> param_objects;

  /* The object representation of the data_set.xml file */
  XmlSubjects subjects;

  /* Looking for all subject with specific object id */
  std::vector<std::string> subjects_filename_by_object_id(const std::string &id) const {
    std::vector<std::string> filenames;

    for (auto &subject : subjects) {
      for (auto &visit : subject.visits) {
        for (auto &object : visit.objects) {
          if (def::support::utilities::strucmp(object.id, id))
            filenames.push_back(object.subject_filename);
        }
      }
    }

    return filenames;
  }

  std::vector<std::string> subjects_class_label() const {
    std::vector<std::string> class_labels;

    for (auto &subject : subjects) {
      class_labels.push_back(subject.class_label);
    }

    return class_labels;
  }

  std::vector<bool> subjects_test_or_train_label_by_object_id(const std::string &id) const {
    std::vector<bool> belongsToTestSet;

    for (auto &subject : subjects) {
      for (auto &visit : subject.visits) {
        for (auto &object : visit.objects) {
          if (def::support::utilities::strucmp(object.id, id))
            belongsToTestSet.push_back((subject.test=="1"));
        }
      }
    }

    return belongsToTestSet;
  }

  /* Looking for all subject with specific object id, visits gathered together */
  std::vector<std::vector<std::string>> subjects_filename_by_object_id_several_visits(const std::string &id) const {
    std::vector<std::vector<std::string>> filenames;

    for (auto &subject : subjects) {
      std::vector<std::string> filenames_subject;
      for (auto &visit : subject.visits) {
        for (auto &object : visit.objects) {
          if (def::support::utilities::strucmp(object.id, id))
            filenames_subject.push_back(object.subject_filename);
        }
      }
      filenames.push_back(filenames_subject);
    }

    return filenames;
  }

  std::set<std::string> extract_object_id() const {
    std::set<std::string> ids;

    for (auto &subject : subjects) {
      for (auto &visit : subject.visits) {
        for (auto &object : visit.objects) {
          ids.emplace(object.id);
        }
      }
    }

    return ids;
  }

  std::vector<std::string> get_subject_ids() const {
    std::vector<std::string> ids;
    for (auto &subject : subjects) {
      ids.push_back(subject.id);
    }
    return ids;
  }

  /* The xml structure of data_set.xml must be the same for each visit and subject : equal number of object_id */
  const bool check_object_id_consistency() const {
    std::set<std::string> objects;

    for (auto &subject : subjects) {
      for (auto &visit : subject.visits) {
        std::set<std::string> current_objects;
        for (auto &object : visit.objects) {
          current_objects.emplace(def::support::utilities::strtolower(object.id));
        }
        if (objects.empty()) {
          objects = current_objects;
        } else {
          if (objects != current_objects) return false;
        }
      }
    }

    return true;
  }

  /* The <age> tag could be required for some computation */
  const bool check_age_consistency() const {

    for (auto &subject : subjects) {
      for (auto &visit : subject.visits) {
        if (!visit.age_assigned) return false;
      }
    }

    return true;
  }

  size_t num_subjects() const { return subjects.size(); }
  size_t num_objects() const { return param_objects.size(); }
};

}
}
