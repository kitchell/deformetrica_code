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

#include "AbstractTestKernelPrecision.h"
#include <src/support/kernels/ExactKernel.h>
#include "src/support/kernels/CUDAExactKernel.h"
#include <random>

namespace def {
namespace test {

class TestKernelPrecisionCUDA : public AbstractTestKernelPrecision {
 public:
  // Constructor : initialize data and kernel used in the tests
  TestKernelPrecisionCUDA() {

    exactKernel2D.SetSources(Y2D);
    exactKernel2D.SetKernelWidth(kernel_width);
    exactKernel3D.SetSources(Y3D);
    exactKernel3D.SetKernelWidth(kernel_width);

    cudaExactKernel2D.SetSources(Y2D);
    cudaExactKernel2D.SetKernelWidth(kernel_width);
    cudaExactKernel3D.SetSources(Y3D);
    cudaExactKernel3D.SetKernelWidth(kernel_width);
  }

 protected:

  ExactKernel<ScalarType, 2> exactKernel2D;
  ExactKernel<ScalarType, 3> exactKernel3D;

  CUDAExactKernel<ScalarType, 2> cudaExactKernel2D;
  CUDAExactKernel<ScalarType, 3> cudaExactKernel3D;

};

}
}
