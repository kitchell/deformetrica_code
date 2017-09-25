/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestReadConfiguration.h"
#include <src/io/XmlDictionary.hpp>
#include "itkObject.h"
#include "itkDOMReader.h"
#include "itkDOMNodeXMLReader.h"
#include "itkMacro.h"
#include "XmlConfigurationConverter.h"
#include <iostream>
#include <src/io/SparseDiffeoParameters.h>
#include <src/io/SparseDiffeoParametersXMLFile.h>
#include <src/io/DeformableObjectParameters.h>
#include <src/io/DeformableObjectParametersXMLFile.h>
#include <memory>
#include <src/io/XmlDataSet.hpp>
//Note:itk::DOM examples http://yanivresearch.info/writtenMaterial/gong2012.pdf
using namespace def::io;

namespace def {
namespace test {


struct xml_parameters {
  SparseDiffeoParameters::Pointer paramDiffeos;
  int numSubjects;
  int numObjects;
  std::vector<DeformableObjectParameters::Pointer> paramObjects;
  std::vector<char*> templatefn;
  std::vector< std::vector< char* > > subjectfn;
};


auto main_sparse_atlas = [](int argc, char *argv[]) {

  if (argc < 6 )
  {
    std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml NumberOfObjects paramsObject1.xml InitialTemplate1 Subject1 Subject2 Subject3 ... paramsObject2.xml InitialTemplate2 Subject1 Subject2 Subject3 ... " << std::endl;
    assert(false);
  }

  int numObjects = atoi(argv[2]);
//  std::cout << "Number of Objects = " << numObjects << std::endl;

  if ( (argc - 3) % numObjects != 0)
  {
    std::cerr << "Number of files mismatched with the number of objects" << std::endl;
    assert(false);
  }

  int numSubjects = (argc - 3) / numObjects - 2;
//  std::cout << "Number of Subjects = " << numSubjects << std::endl;

  if (numSubjects < 2)
  {
    std::cerr << "Atlas building requires at least 2 subjects" << std::endl;
    assert(false);
  }


  // read general parameters for diffeomorphic matching
  SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
  if (paramDiffeos.IsNull())
  {
    std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
    assert(false);
  }

  // read the list of objects: params, source and target
  std::vector<char*> templatefn;
  std::vector< std::vector< char* > > subjectfn;
  std::vector<DeformableObjectParameters::Pointer> paramObjects;

  templatefn.resize(numObjects);
  paramObjects.resize(numObjects);
  subjectfn.resize(numObjects);
  for (int i = 0; i < numObjects; i++)
    subjectfn[i].resize(numSubjects);

  unsigned int indx = 3;
  for (unsigned int i = 0; i < numObjects; i++)
  {
    paramObjects[i] = readDeformableObjectParametersXML(argv[indx++]);
    if (paramObjects[i].IsNull())
    {
      std::cerr << "Failed creating XML object, bad input file " << argv[indx] << "?" << std::endl;
      assert(false);
    }
    templatefn[i] = argv[indx++];
    for (unsigned s = 0; s < numSubjects; s++)
      subjectfn[i][s] = argv[indx++];
  }

  xml_parameters parameters;
  parameters.paramDiffeos = std::move(paramDiffeos);
  parameters.numSubjects = numSubjects;
  parameters.numObjects = numObjects;
  parameters.paramObjects = std::move(paramObjects);
  parameters.templatefn = std::move(templatefn);
  parameters.subjectfn = std::move(subjectfn);

  return std::move(parameters);
};

void TestReadConfiguration::SetUp() {
  Test::SetUp();
}


TEST_F(TestReadConfiguration, regression_multiple_objects_and_subjects) {
  std::unique_ptr<def::io::XmlConfigurationRequirements> xml_converter(
      new def::io::XmlConfigurationConverter(UNIT_TESTS_DIR"/io/data/generic_multiple_objects_and_subjects/model.xml",
                                             UNIT_TESTS_DIR"/io/data/generic_multiple_objects_and_subjects/data_set.xml")
  );

  return;
  auto xml_model = xml_converter->generate_xml_model();
  XmlSubjects& subjects = xml_model->subjects;

  ASSERT_EQ(subjects.size(), 3);
  ASSERT_EQ(subjects[0].id, "subj1");
  ASSERT_EQ(subjects[1].id, "subj2");
  ASSERT_EQ(subjects[2].id, "subj3");

  ASSERT_EQ(subjects[0].visits.size(), 2);
  ASSERT_EQ(subjects[1].visits.size(), 4);
  ASSERT_EQ(subjects[2].visits.size(), 1);

  ASSERT_EQ(subjects[0].visits[0].age, 65.2);
  ASSERT_EQ(subjects[0].visits[0].objects[0].id, "bundle");
  ASSERT_EQ(subjects[0].visits[0].objects[1].id, "cortico");
  ASSERT_EQ(subjects[0].visits[0].objects[0].subject_filename, "x1");
  ASSERT_EQ(subjects[0].visits[0].objects[1].subject_filename, "x2");

  ASSERT_EQ(subjects[1].visits[2].age, 62.8);
  ASSERT_EQ(subjects[1].visits[2].objects[0].id, "bundle");
  ASSERT_EQ(subjects[1].visits[2].objects[0].subject_filename, "x9");

  ASSERT_EQ(xml_model->param_objects.size(), 2);
  ASSERT_EQ(xml_model->param_objects["bundle"]->GetFilename(), "Bundle_Left_Cortico_Putamen/Template_Bundle_cortico_putamen.vtk");
  ASSERT_EQ(xml_model->param_objects["bundle"]->GetDataSigma(), 0.5);
  ASSERT_EQ(xml_model->param_objects["cortico"]->GetFilename(), "Left_Putamen/template_surface_put.vtk");
  ASSERT_EQ(xml_model->param_objects["cortico"]->GetDataSigma(), 0.6);
}


TEST_F(TestReadConfiguration, sparse_atlas) {

  std::unique_ptr<def::io::XmlConfigurationRequirements> xml_converter(
      new def::io::XmlConfigurationConverter(UNIT_TESTS_DIR"/io/data/sparse_atlas/model.xml",
                                             UNIT_TESTS_DIR"/io/data/sparse_atlas/data_set.xml",
                                             UNIT_TESTS_DIR"/io/data/sparse_atlas/optimization_parameters.xml")
  );

  int num_subjects = xml_converter->num_subjects();
  int num_objects = xml_converter->num_objects();

  ASSERT_EQ(num_subjects, 2);
  ASSERT_EQ(num_objects, 2);

  auto&& param_diffeos = xml_converter->param_diffeos();
  auto&& templatefn = xml_converter->templates();
  auto&& subjectfn = xml_converter->subjects();
  auto&& observations = xml_converter->observations();
  auto&& param_objects = xml_converter->param_objects();

  ASSERT_EQ(param_diffeos->GetInitialCPSpacing(), 1);
  ASSERT_EQ(param_diffeos->GetInitialCPPosition_fn(), "filename");
  ASSERT_EQ(param_diffeos->GetKernelType(), "exact");
  ASSERT_EQ(param_diffeos->GetOptimizationMethodType(), "FastGradientAscent");
  ASSERT_EQ(param_diffeos->GetInitialStepSize(), 0.01);

  ASSERT_EQ(templatefn.size(), 2);
  ASSERT_EQ(templatefn[0], "Bundle_Left_Cortico_Putamen/Template_Bundle_cortico_putamen.vtk");
  ASSERT_EQ(templatefn[1], "Left_Putamen/template_surface_put.vtk");

  ASSERT_EQ(subjectfn.size(), 2);
  ASSERT_EQ(subjectfn[0].size(), 2);
  ASSERT_EQ(subjectfn[0][0], "Data/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj1_MP_std2_perc15.vtk");
  ASSERT_EQ(subjectfn[0][1], "Data/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj2_MP_std2_perc15.vtk");
  ASSERT_EQ(subjectfn[1][0], "Data/Left_Putamen/ngc_1002_CT_DS_RT_subj_1.vtk");
  ASSERT_EQ(subjectfn[1][1], "Data/Left_Putamen/ngc_1002_CT_DS_RT_subj_2.vtk");

/*
  ASSERT_EQ(observations.size(),2);
  ASSERT_EQ(observations[0].size(),2);
  ASSERT_EQ(observations[0][0], 65.2);
  ASSERT_EQ(observations[0][1], 72.3);
  ASSERT_EQ(observations[1].size(),2);
  ASSERT_EQ(observations[1][0], 65.2);
  ASSERT_EQ(observations[1][1], 72.3);
*/

  ASSERT_EQ(param_objects.size(), 2);
  ASSERT_EQ(param_objects[0]->GetDeformableObjectType(), "OrientedSurfaceMesh");
  ASSERT_EQ(param_objects[0]->GetDataSigma(), 0.5);

  ASSERT_EQ(param_objects[1]->GetDeformableObjectType(), "OrientedCurveMesh");
  ASSERT_EQ(param_objects[1]->GetDataSigma(), 0.6);
}

TEST_F(TestReadConfiguration, data_set) {

  itk::DOMNode::Pointer dom;
  itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
  reader->SetFileName(UNIT_TESTS_DIR"/io/data/sparse_atlas/data_set.xml");
  reader->Update();
  dom = reader->GetOutput();

  XmlDictionary xml;

  itk::DOMNode::ChildrenListType subjects;
  dom->GetChildren("subject", subjects);

  std::map<std::string, std::vector<std::string> > files;


  ASSERT_EQ(subjects.size(), 2);
  for (auto it_subject : subjects) {
    itk::DOMNode::ChildrenListType visits;
    it_subject->GetChildren("visit", visits);
    ASSERT_LE(visits.size(), 2);

    for (auto it_visit : visits) {
      /* AGE */
      auto age = it_visit->GetChild("age");
      if (age) {
        //TODO: store age
      }

      itk::DOMNode::ChildrenListType filenames;
      it_visit->GetChildren("filename", filenames);
      ASSERT_EQ(filenames.size(), 2);

      for (auto it_filename : filenames) {
        auto id = it_filename->GetAttribute("object_id");
        if (id.empty())
          throw std::runtime_error("Missing 'object_id' xml tag on configuration file");

        files[id].push_back(it_filename->GetTextChild()->GetText());
      }
      //TODO: adding more timepoint visit for the future..
      //TODO: check xml structure consistencies
      break;
    }
  }

  ASSERT_EQ(files["bundle"][0], "Data/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj1_MP_std2_perc15.vtk");
  ASSERT_EQ(files["bundle"][1], "Data/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj2_MP_std2_perc15.vtk");
  ASSERT_EQ(files["cortico"][0], "Data/Left_Putamen/ngc_1002_CT_DS_RT_subj_1.vtk");
  ASSERT_EQ(files["cortico"][1], "Data/Left_Putamen/ngc_1002_CT_DS_RT_subj_2.vtk");

  /* BUILDING DATA */
  int num_subjects = subjects.size();
  int num_objects = files.size();

  ASSERT_EQ(num_subjects, 2);
  ASSERT_EQ(num_objects, 2);

  int counter = 0;
  std::vector<std::string> templatefn(num_objects);
  std::vector< std::vector<std::string> > subjectfn(num_objects);
  std::vector<DeformableObjectParameters::Pointer> paramObjects;

  XmlDictionary templates;

  for (auto it : files) {
    auto tag = it.first;
    auto& object = templates["template"]["object"].add_tag("id", tag);

    object["filename"].add_reader([&, tag](const std::string &v) {
      templatefn[counter] = v;
      for (auto file : files[tag])
        subjectfn[counter].push_back(file);
      ++counter;
    });
  }

  /* LOADING MODEL TEMPLATES */
  {
    itk::DOMNode::Pointer dom;
    itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
    reader->SetFileName(UNIT_TESTS_DIR"/io/data/sparse_atlas/model.xml");
    reader->Update();
    dom = reader->GetOutput();

    itk::DOMNode::ChildrenListType children;
    dom->GetAllChildren(children);
    XmlDictionary::xml_reader(dom, templates);
  }

  ASSERT_EQ(templatefn[0], "Bundle_Left_Cortico_Putamen/Template_Bundle_cortico_putamen.vtk");
  ASSERT_EQ(templatefn[1], "Left_Putamen/template_surface_put.vtk");

  ASSERT_EQ(subjectfn[0].size(), 2);
  ASSERT_EQ(subjectfn[0][0], "Data/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj1_MP_std2_perc15.vtk");
  ASSERT_EQ(subjectfn[0][1], "Data/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj2_MP_std2_perc15.vtk");

  ASSERT_EQ(subjectfn[1].size(), 2);
  ASSERT_EQ(subjectfn[1][0], "Data/Left_Putamen/ngc_1002_CT_DS_RT_subj_1.vtk");
  ASSERT_EQ(subjectfn[1][1], "Data/Left_Putamen/ngc_1002_CT_DS_RT_subj_2.vtk");
}

TEST_F(TestReadConfiguration, load_template) {

  itk::DOMNode::Pointer dom;
  itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
  reader->SetFileName(UNIT_TESTS_DIR"/io/data/sparse_atlas/model.xml");
  reader->Update();
  dom = reader->GetOutput();

  itk::DOMNode::ChildrenListType children;
  dom->GetAllChildren(children);

  std::vector<DeformableObjectParameters::Pointer> paramObjects;

  std::map<std::string, std::string> file_map;

  XmlDictionary xml;

  for (auto tag : {"bundle","cortico"}) {
    DeformableObjectParameters::Pointer param = DeformableObjectParameters::New();

    auto& object = xml["template"]["object"].add_tag("id", tag);
    object["deformable-object-type"].assign_to<std::string>(param,&DeformableObjectParameters::SetDeformableObjectType);
    object["data-sigma"].assign_to<double>(param,&DeformableObjectParameters::SetDataSigma);
    object["kernel-width"].assign_to<double>(param,&DeformableObjectParameters::SetKernelWidth);
    object["kernel-type"].assign_to<std::string>(param,&DeformableObjectParameters::SetKernelType);
    object["reorient-normals"]
        .filter_with(def::io::filters::lower_case())
        .add_reader([param](const std::string &v) {
          v == "on" ? assert(true) : assert(false);
        });
    object["filename"].add_reader([&file_map, tag](const std::string &v) {
      file_map[tag] = v;
    });

    paramObjects.push_back(std::move(param));
  }

  XmlDictionary bundle, cortico;
  ASSERT_TRUE(xml["template"]["object"].load_tag("id","bundle", bundle));
  ASSERT_TRUE(xml["template"]["object"].load_tag("id","cortico", cortico));

  xml["template"]["object"]["test"]["x"].add_reader([](const std::string&v){
    std::cout  << "__TESST__";
  });
//  ASSERT_EQ(bundle.size(), 6);

  XmlDictionary::xml_reader(dom, xml);

  ASSERT_EQ(paramObjects[0]->GetDeformableObjectType(), "OrientedSurfaceMesh");
  ASSERT_EQ(paramObjects[0]->GetDataSigma(), 0.5);
  ASSERT_EQ(file_map["bundle"], "Bundle_Left_Cortico_Putamen/Template_Bundle_cortico_putamen.vtk");
}

TEST_F(TestReadConfiguration, dictionary) {
  itk::DOMNode::Pointer dom;
  itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
  reader->SetFileName(UNIT_TESTS_DIR"/io/data/sparse_atlas/model.xml");
  reader->Update();
  dom = reader->GetOutput();

  itk::DOMNode::ChildrenListType children;
  dom->GetAllChildren(children);

  XmlDictionary xml;
  SparseDiffeoParameters::Pointer sp = SparseDiffeoParameters::New();

  XmlConfigurationConverter converter;
  converter.programming_xml_model(xml, sp);

//  ASSERT_EQ(xml.size(), 6); // 4 ?
//  ASSERT_EQ(xml["deformation-parameters"].size(), 5);

  XmlDictionary::xml_reader(dom, xml);

  ASSERT_EQ(sp->GetInitialCPSpacing(), 1);
  ASSERT_EQ(sp->GetInitialCPPosition_fn(), "filename");
  ASSERT_EQ(sp->GetKernelWidth(), 1);
  ASSERT_EQ(sp->GetKernelType(), "exact");
  ASSERT_EQ(sp->GetT0(), 0.0);
  ASSERT_EQ(sp->GetTN(), 1.0);
  ASSERT_EQ(sp->GetNumberOfTimePoints(), 10);
}


TEST_F(TestReadConfiguration, compare_oldatlas_cfiguratio) {

  // read a DOM object from an XML file
  try
  {
    char const *cmdline[] = {
        "deformetrica-run", UNIT_TESTS_DIR"/io/data/sparse_atlas/paramDiffeos.xml", "2",
        UNIT_TESTS_DIR"/io/data/sparse_atlas/Bundle_Left_Cortico_Putamen/param_cortico_putamen_1002.xml", UNIT_TESTS_DIR"/data/sparse_atlas/Bundle_Left_Cortico_Putamen/Template_Bundle_cortico_putamen.vtk", UNIT_TESTS_DIR"/data/sparse_atlas/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj1_MP_std2_perc15.vtk", UNIT_TESTS_DIR"/data/sparse_atlas/Bundle_Left_Cortico_Putamen/lpu_to_pial_subj2_MP_std2_perc15.vtk",
        UNIT_TESTS_DIR"/io/data/sparse_atlas/Left_Putamen/param_putamen.xml", UNIT_TESTS_DIR"/data/sparse_atlas/Left_Putamen/template_surface_put.vtk", UNIT_TESTS_DIR"/data/sparse_atlas/Left_Putamen/ngc_1002_CT_DS_RT_subj_1.vtk", UNIT_TESTS_DIR"/data/sparse_atlas/Left_Putamen/ngc_1002_CT_DS_RT_subj_2.vtk"
    };
    char **argv = const_cast<char**>(cmdline);
    int argc = sizeof(cmdline) / sizeof(cmdline[0]);

    xml_parameters&& params = main_sparse_atlas(argc, argv);

    def::io::XmlConfigurationConverter xml_converter;
    xml_converter.load_xml_model(UNIT_TESTS_DIR"/io/data/sparse_atlas/model.xml");
    xml_converter.load_xml_data_set(UNIT_TESTS_DIR"/io/data/sparse_atlas/data_set.xml");
    xml_converter.load_xml_optimization_parameters(UNIT_TESTS_DIR"/io/data/sparse_atlas/optimization_parameters.xml");

    auto&&t = xml_converter.param_diffeos();

  }
  catch (const itk::ExceptionObject& e)
  {
    ASSERT_TRUE(false);
  }
  catch (...)
  {
    ASSERT_TRUE(false);
  }

  ASSERT_TRUE(true);
}

TEST_F(TestReadConfiguration, atlas_3d) {

  def::io::XmlConfigurationConverter xml_converter;
  xml_converter.load_xml_model(UNIT_TESTS_DIR"/io/data/atlas_3d/model.xml");
  xml_converter.load_xml_data_set(UNIT_TESTS_DIR"/io/data/atlas_3d/data_set.xml");
  xml_converter.load_xml_optimization_parameters(UNIT_TESTS_DIR"/io/data/atlas_3d/optimization_parameters.xml");
  xml_converter.validate_configuration();

  ASSERT_EQ(xml_converter.num_subjects(), 3);
  ASSERT_EQ(xml_converter.num_objects(), 1);

}

TEST_F(TestReadConfiguration, xml_model) {

  // read a DOM object from an XML file
  try
  {

    itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
    reader->SetFileName(UNIT_TESTS_DIR"/io/data/conf/data_set.xml");
    reader->Update();
    itk::DOMNode::Pointer dom = reader->GetModifiableOutput();

    ASSERT_EQ(dom->GetName(), "data-set");

    itk::DOMNode::ChildrenListType subjects;
    dom->GetChildren("subject", subjects);
    ASSERT_EQ(subjects.size(), 2);

    for (auto & subject : subjects) {
      itk::DOMNode::ChildrenListType visit;
      subject->GetChildren("visit", visit);

      ASSERT_GE(visit.size(), 2);
      ASSERT_LE(visit.size(), 3);
    }

    auto & subject_1 = subjects[0];
    itk::DOMNode::AttributesListType x;
    itk::DOMNode::ConstChildrenListType data;

    auto gender = subject_1->GetChild("gender");
    ASSERT_EQ(gender->GetTextChild()->GetText(), "male");

    itk::DOMNode::ChildrenListType sub_visits;
    subject_1->GetChildren("visit", sub_visits);

    auto & visit_1 = sub_visits[0];
    ASSERT_EQ(visit_1->GetAttribute("id"), "baseline");

    int age_int;
    itk::FancyString age_text = visit_1->GetChild("age")->GetTextChild()->GetText();
    age_text >> age_int;
    ASSERT_EQ(age_int, 63);

  }
  catch ( const itk::ExceptionObject& espar )
  {
    ASSERT_TRUE(false);
  }
  catch ( ... )
  {
    ASSERT_TRUE(false);
  }

  ASSERT_TRUE(true);
}

}
}
