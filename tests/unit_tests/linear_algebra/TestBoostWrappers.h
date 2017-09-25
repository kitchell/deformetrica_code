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

#include "gtest/gtest.h"

/// Support files.
#include "LinearAlgebra.h"

/// Input-output files.
#include "MatrixDLM.h"

/// Librairies files.
#include <typeinfo>
#include <map>
#include <iostream>
#include <cstring>

namespace def {
namespace test {

class TestBoostWrappers : public ::testing::Test {
 public:

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif

 protected:
  virtual void SetUp();
};
}
}