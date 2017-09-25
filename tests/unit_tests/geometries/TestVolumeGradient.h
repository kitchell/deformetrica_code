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
#include "DeformetricaConfig.h"
#include "OrientedVolumeMesh.h"
#include <vtkUnstructuredGridReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>

using namespace def::algebra;

namespace def {
namespace test {

    class TestVolumeGradient : public ::testing::Test {
    public:

#ifdef USE_DOUBLE_PRECISION
        typedef double ScalarType;
#else
        typedef float ScalarType;
#endif

        template<unsigned int Dimension> std::shared_ptr<OrientedVolumeMesh<ScalarType, Dimension>> ReadOrientedVolumeMesh(const char *filePath);
        MatrixType GetNumericalGradient3(MatrixType &sourcePoints, std::shared_ptr<OrientedVolumeMesh<ScalarType, 3>> source, const std::shared_ptr<OrientedVolumeMesh<ScalarType, 3>> target);
        MatrixType GetNumericalGradient2(MatrixType &sourcePoints, std::shared_ptr<OrientedVolumeMesh<ScalarType, 2>> source, const std::shared_ptr<OrientedVolumeMesh<ScalarType, 2>> target);

    protected:
        virtual void SetUp();
    };
}
}
