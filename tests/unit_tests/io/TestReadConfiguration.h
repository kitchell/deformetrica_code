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

namespace def {
namespace test {

    class TestReadConfiguration : public ::testing::Test {
    protected:
        virtual void SetUp();
    };
}
}