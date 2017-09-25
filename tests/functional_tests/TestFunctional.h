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

extern std::vector<std::string> args;

namespace def {
namespace test {

    class TestFunctional : public ::testing::Test {
    protected:
        virtual void SetUp();
    };
}
}