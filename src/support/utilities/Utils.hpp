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

#include "itksys/SystemTools.hxx"

namespace def {
namespace support {
namespace utilities {

auto strucmp = [](const std::string &v1, const std::string &v2) {
  return (itksys::SystemTools::LowerCase(v1) == itksys::SystemTools::LowerCase(v2));
};

auto strtolower = [](const std::string &v) {
  return itksys::SystemTools::LowerCase(v);
};

auto strtoupper = [](const std::string &v) {
  return itksys::SystemTools::UpperCase(v);
};

}
}
}

