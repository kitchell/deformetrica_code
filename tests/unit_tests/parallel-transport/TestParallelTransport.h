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
#include "itkImageFileReader.h"

#include <src/core/observations/deformable_objects/geometries/AbstractGeometry.h>
#include <src/support/linear_algebra/LinearAlgebra.h>


namespace def {
namespace test {

    class TestParallelTransport : public ::testing::Test {
    public:

#ifdef USE_DOUBLE_PRECISION
        typedef double ScalarType;
#else
        typedef float ScalarType;
#endif
        typedef def::algebra::MatrixType MatrixType;
        typedef def::algebra::MatrixListType MatrixListType;

        //We read an image in the example.
        typedef itk::Image<double, 2> ImageType;
        typedef itk::ImageFileReader<ImageType> ImageReaderType;
        typedef typename ImageType::Pointer ImageTypePointer;
        typedef itk::ImageFileReader<ImageType> ReaderType;
        typedef AbstractGeometry<double, 2> AbstractGeometryType;
        typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryList;

 protected:
        virtual void SetUp();
    };
}
}
