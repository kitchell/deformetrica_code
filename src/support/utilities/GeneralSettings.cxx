/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "GeneralSettings.h"

namespace def {
namespace utils {

  GeneralSettings settings = SingletonGeneralSettings::instance()->settings();

  SingletonGeneralSettings* SingletonGeneralSettings::instance() {
    static SingletonGeneralSettings* instance = nullptr;
    if (instance) return instance;
    return (instance = new SingletonGeneralSettings());
  }

  GeneralSettings& SingletonGeneralSettings::settings() {
    return settings_;
  }


}
}

