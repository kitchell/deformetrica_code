/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _GeneralSettings_h_
#define _GeneralSettings_h_

#include <string>

namespace def {
namespace utils {

struct GeneralSettings {
  bool save_state;
  bool load_state;
  std::string input_state_filename;
  std::string output_state_filename;
  std::string output_dir;
  unsigned int number_of_threads;
};

class SingletonGeneralSettings {
 public:
  SingletonGeneralSettings(const SingletonGeneralSettings&) = delete;
  SingletonGeneralSettings(SingletonGeneralSettings&&) = delete;

  SingletonGeneralSettings& operator=(const SingletonGeneralSettings&) = delete;
  SingletonGeneralSettings& operator=(SingletonGeneralSettings&&) = delete;

  static SingletonGeneralSettings* instance();
  GeneralSettings& settings();

 protected:
  SingletonGeneralSettings() {}
  ~SingletonGeneralSettings() {}

  GeneralSettings settings_;
};

extern GeneralSettings settings;
}
}

#endif