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

#include "LinearAlgebra.h"
#include "MatrixDLM.h"
#include <typeinfo>

namespace def {
    namespace test {

        class TestFloat : public ::testing::Test {
        public:

#ifdef USE_DOUBLE_PRECISION
            typedef double ScalarType;
#else
            typedef float ScalarType;
#endif
            typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
            typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;

        protected:
            virtual void SetUp();
        };
    }
}