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

#include <src/io/XmlConfigurationConverter.h>
#include <src/launch/abc/abc_sampling.cxx>
#include <src/launch/atlas/estimate_atlas.h>
#include <src/launch/longitudinal_atlas/estimate_longitudinal_atlas.h>
#include <src/launch/longitudinal_registration/estimate_longitudinal_registration.h>
#include <src/launch/shooting/deform.h>
#include <src/launch/regression/estimate_geodesic_regression.h>
#include <src/launch/parallel_transport/parallel_transport.h>
#include <src/support/utilities/Utils.hpp>
#include <src/support/utilities/GeneralSettings.h>
#include <boost/exception/all.hpp>

enum type {
  _NO_TYPE_,
  _REGISTRATION_,
  _ATLAS_,
  _LONGITUDINAL_,
  _LONGITUDINAL_REGISTRATION_,
  _REGRESSION_,
  _SHOOTING_,
  _PARALLEL_TRANSPORT_,
  _ABC_
};

enum dimension {
  _NO_DIMENSION_,
  _2D_ = 2,
  _3D_
};

enum cmdline_position {
  _PROGRAM_NAME_,
  _TYPE_,
  _DIMENSION_,
  _MODEL_,
  _DATA_SET_,
  _OPTIMIZATION_PARAMETERS_
};

enum cmdline_options {
  _INPUT_STATE_FILE_,
  _OUTPUT_STATE_FILE_,
  _OUTPUT_DIR_
};

void deformetrica(int argc, char **argv) {

  auto cmdline_assert = [&argv](bool assert, std::string msg = std::string()) {
    if (assert) return;

    if (msg.size())
      std::cerr << msg;

    std::cerr << std::endl
              << "Usage: "
              << argv[_PROGRAM_NAME_]
              << " {registration, atlas, regression, longitudinal, longitudinal-registration, parallel-transport} "
              "{2D, 3D} <model.xml> <data_set.xml> <optimization_parameters.xml> "
              "[--input-state-file=<filename.bin>] [--output-state-file=<filename.bin>] [--save-period=<integer>] "
              "[--output-dir=<path>]"
              << std::endl;

    exit(-1);
  };

  /* Looking for optional arguments */
  std::map<unsigned int, std::string> s_argv;
  std::map<unsigned int, std::string> s_opt;

  std::map<std::string, int> index;
  index["--input-state-file="] = _INPUT_STATE_FILE_;
  index["--output-state-file="] = _OUTPUT_STATE_FILE_;
  index["--output-dir="] = _OUTPUT_DIR_;

  std::for_each(argv, argv + argc, [&](char *v) {
    std::string s(v);
    int i = s_argv.size();
    int j = s_opt.size();
    for (std::string op : {"--input-state-file=", "--output-state-file=", "--output-dir="}) {
      if (std::string::npos != s.find(op)) {
        s_opt[index[op]] = s.erase(0, op.size());
        return;
      }
    }

    s_argv[i] = s;
  });

  cmdline_assert(s_argv.size() >= 5 && s_argv.size() <= 9, "Error: Wrong arguments numbers");

  def::utils::settings.save_state = true;
  def::utils::settings.load_state = false;
  def::utils::settings.output_state_filename = "deformetrica-state.bin";
  def::utils::settings.output_dir = "./";

  if (s_opt.size()) {
    if (s_opt.find(_INPUT_STATE_FILE_) != s_opt.end()) {
      def::utils::settings.input_state_filename = s_opt[_INPUT_STATE_FILE_];
      def::utils::settings.load_state = true;
    }

    if (s_opt.find(_OUTPUT_STATE_FILE_) != s_opt.end()) {
      def::utils::settings.output_state_filename = s_opt[_OUTPUT_STATE_FILE_];
    }

    if (s_opt.find(_OUTPUT_DIR_) != s_opt.end()) {
      def::utils::settings.output_dir = s_opt[_OUTPUT_DIR_] + "/";
      if (s_opt.find(_OUTPUT_STATE_FILE_) == s_opt.end()) {
        def::utils::settings.output_state_filename =
            def::utils::settings.output_dir + def::utils::settings.output_state_filename;
      }
    }
  }

  auto &cmp = def::support::utilities::strucmp;
  std::map<decltype(_NO_TYPE_), std::string> type_algo;
  type_algo[_REGISTRATION_] = "registration";
  type_algo[_ATLAS_] = "atlas";
  type_algo[_LONGITUDINAL_] = "longitudinal";
  type_algo[_LONGITUDINAL_REGISTRATION_] = "longitudinal-registration";
  type_algo[_REGRESSION_] = "regression";
  type_algo[_SHOOTING_] = "shooting";
  type_algo[_PARALLEL_TRANSPORT_] = "parallel-transport";
  type_algo[_ABC_] = "abc";

  type algo = _NO_TYPE_;
  for (auto it : type_algo) {
    if (cmp(it.second, s_argv[_TYPE_])) {
      algo = it.first;
      break;
    }
  }

  cmdline_assert(algo != _NO_TYPE_, "Error: available type are 'registration', 'atlas', 'longitudinal', "
      "'longitudinal-registration', 'regression', 'parallel-transport' or 'shooting'");

  std::map<decltype(_NO_DIMENSION_), std::string> type_dim;
  type_dim[_2D_] = "2D";
  type_dim[_3D_] = "3D";

  dimension algo_dimension = _NO_DIMENSION_;
  for (auto it : type_dim) {
    if (cmp(it.second, s_argv[_DIMENSION_])) {
      algo_dimension = it.first;
      break;
    }
  }

  cmdline_assert(algo_dimension != _NO_DIMENSION_);

  const bool argv_has_optimization = (s_argv.size() == 6);

  def::io::XmlConfigurationConverter xml_converter;
  xml_converter.load_xml_model(s_argv[_MODEL_]);
  xml_converter.load_xml_data_set(s_argv[_DATA_SET_]);
  if (argv_has_optimization)
    xml_converter.load_xml_optimization_parameters(s_argv[_OPTIMIZATION_PARAMETERS_]);
  xml_converter.validate_configuration();

  auto xml_model = xml_converter.generate_xml_model();
  int num_subjects = xml_converter.num_subjects();
  cmdline_assert(num_subjects > 0, "Error: missing subjects");

  /// Sets some additional general settings.
  def::utils::settings.number_of_threads = xml_model->param_diffeos->GetNumberOfThreads();

  std::map<type, std::function<void()> > run;

  run[_REGISTRATION_] = [&]() {
    cmdline_assert(num_subjects == 1, "Registration requires only one subject");
    auto atlas_run = algo_dimension == _2D_ ? estimateAtlas<2> : estimateAtlas<3>;
    atlas_run(xml_model);
  };

  run[_ATLAS_] = [&]() {
    cmdline_assert(num_subjects > 1, "Atlas building requires at least two subjects");
    auto atlas_run = algo_dimension == _2D_ ? estimateAtlas<2> : estimateAtlas<3>;
    atlas_run(xml_model);
  };

  run[_LONGITUDINAL_] = [&]() {
    cmdline_assert(num_subjects > 1, "Longitudinal atlas building requires at least two subjects");
    auto longitudinal_run = algo_dimension == _2D_ ? estimateLongitudinalAtlas<2> : estimateLongitudinalAtlas<3>;
    longitudinal_run(xml_model);
  };

  run[_LONGITUDINAL_REGISTRATION_] = [&]() {
    cmdline_assert(num_subjects == 1, "Longitudinal registration requires only one subject");
    auto longitudinal_registration_run =
        algo_dimension == _2D_ ? estimateLongitudinalRegistration<2> : estimateLongitudinalRegistration<3>;
    longitudinal_registration_run(xml_model);
  };

  run[_REGRESSION_] = [&]() {
    cmdline_assert(num_subjects == 1, "Regression requires only one subject");
    auto regression_run = algo_dimension == _2D_ ? estimateGeodesicRegression<2> : estimateGeodesicRegression<3>;
    regression_run(xml_model);
  };

  run[_PARALLEL_TRANSPORT_] = [&]() {
    cmdline_assert(num_subjects >= 1, "Parallel transport requires at least one subject");
    auto parallel_transport_run = algo_dimension == _2D_ ? parallel_transport<2> : parallel_transport<3>;
    parallel_transport_run(xml_model);
  };

  run[_SHOOTING_] = [&]() {
    auto shooting_run = algo_dimension == _2D_ ? deformation<2> : deformation<3>;
    shooting_run(xml_model);
  };

  run[_ABC_] = [&]() {
    auto abc_run = algo_dimension == _2D_ ? abc_sampling<2> : abc_sampling<3>;
    abc_run(xml_model);
  };

  run[algo]();
}
