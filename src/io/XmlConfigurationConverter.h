/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef DEFORMETRICA_XMLCONFIGURATIONCONVERTER_H
#define DEFORMETRICA_XMLCONFIGURATIONCONVERTER_H

#include "itkDOMNodeXMLReader.h"
#include <src/io/SparseDiffeoParameters.h>
#include <src/io/SparseDiffeoParametersXMLFile.h>
#include <src/io/DeformableObjectParameters.h>
#include <src/io/DeformableObjectParametersXMLFile.h>
#include <src/io/XmlDictionary.hpp>
#include <src/io/XmlDataSet.hpp>

namespace def {
namespace io {


class XmlConfigurationRequirements {
 public:
  ~XmlConfigurationRequirements() {}

  virtual SparseDiffeoParameters::Pointer param_diffeos() = 0;
  virtual int num_subjects() = 0;
  virtual int num_objects() = 0;
  virtual std::vector<DeformableObjectParameters::Pointer> param_objects() = 0;
  virtual std::vector<std::string> templates() = 0;
  virtual std::vector< std::vector<std::string> > subjects() = 0;
  virtual std::vector< std::vector<double> > observations() = 0;
  virtual XmlSubjects data_set_subjects() = 0;
  virtual std::shared_ptr<XmlModel> generate_xml_model() = 0;
};

class XmlConfigurationConverter : public XmlConfigurationRequirements {
 public:
  XmlConfigurationConverter() : num_subjects_(0), num_objects_(0), has_optimication_parameters_(false)
  {
    xml_model_ = std::make_shared<XmlModel>();
  }

  XmlConfigurationConverter(const std::string& model,
                            const std::string& data_set)
      : num_subjects_(0), num_objects_(0), has_optimication_parameters_(false)
  {
    xml_model_ = std::make_shared<XmlModel>();
    load_xml_model(model);
    load_xml_data_set(data_set);
    validate_configuration();
  }

  XmlConfigurationConverter(const std::string& model,
                            const std::string& data_set,
                            const std::string& optimization_parameters)
      : num_subjects_(0), num_objects_(0), has_optimication_parameters_(false)
  {
    xml_model_ = std::make_shared<XmlModel>();
    load_xml_model(model);
    load_xml_data_set(data_set);
    load_xml_optimization_parameters(optimization_parameters);
    validate_configuration();
  }

  void load_xml_model(const std::string filename) throw(std::runtime_error);
  void load_xml_data_set(const std::string& filename) throw(std::runtime_error);
  void load_xml_optimization_parameters(const std::string& filename) throw(std::runtime_error);
  void validate_configuration() throw(std::runtime_error);
  virtual std::shared_ptr<XmlModel> generate_xml_model() throw(std::runtime_error);

  virtual SparseDiffeoParameters::Pointer param_diffeos() { return param_diffeos_; }
  virtual int num_subjects() { return num_subjects_; }
  virtual int num_objects() { return num_objects_; }
  virtual std::vector<DeformableObjectParameters::Pointer> param_objects() { return param_objects_; }
  virtual std::vector<std::string> templates() { return templatefn_; }
  virtual std::vector< std::vector<std::string> > subjects() { return subjectfn_; }
  virtual std::vector< std::vector<double> > observations() { return observations_; }
  virtual XmlSubjects data_set_subjects();

  void programming_xml_model(XmlDictionary& xml, SparseDiffeoParameters::Pointer sp);
  void programming_xml_optimization_parameters(XmlDictionary& xml, SparseDiffeoParameters::Pointer sp);

 protected:
  void programming_dynamic_xml_model_object(XmlDictionary& templates);
  void programming_dynamic_xml_model_templates(XmlDictionary& templates);

 private:
  void load_xml_file(itk::DOMNodeXMLReader::Pointer& reader, const std::string& filename) throw(std::runtime_error);

 private:
  itk::DOMNodeXMLReader::Pointer itk_xml_model_;
  itk::DOMNodeXMLReader::Pointer itk_xml_data_set_;
  itk::DOMNodeXMLReader::Pointer itk_xml_optimization_parameters_;

  int num_subjects_;
  int num_objects_;
  size_t objects_order_;
  bool has_optimication_parameters_;
  std::vector<std::string> templatefn_;
  std::vector< std::vector<std::string> > subjectfn_;
  std::vector< std::vector<double> > observations_;
  std::vector<DeformableObjectParameters::Pointer> param_objects_;
  SparseDiffeoParameters::Pointer param_diffeos_;
  std::shared_ptr<XmlModel> xml_model_;
};

}
}


#endif //DEFORMETRICA_XMLCONFIGURATIONCONVERTER_H
