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
#include "src/support/kernels/P3MKernel.h"
#include <random>

namespace def {
namespace test {

class TestKernelPrecisionP3M : public AbstractTestKernelPrecision {
 public:
  // Constructor : initialize data and kernel used in the tests
  TestKernelPrecisionP3M() {

    exactKernel2D.SetSources(Y2D);
    exactKernel2D.SetKernelWidth(kernel_width);
    exactKernel3D.SetSources(Y3D);
    exactKernel3D.SetKernelWidth(kernel_width);

    // Determine the bounding box and other quantities for P3M kernels
    MatrixType DataDomain(3, 2, 0.0f);
    for (int k = 0; k < X3D.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (X3D(k, d) < DataDomain(d, 0) ? X3D(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (X3D(k, d) > DataDomain(d, 1) ? X3D(k, d) : DataDomain(d, 1));
      }
    }

    for (int k = 0; k < Y3D.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (Y3D(k, d) < DataDomain(d, 0) ? Y3D(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (Y3D(k, d) > DataDomain(d, 1) ? Y3D(k, d) : DataDomain(d, 1));
      }
    }

    p3mKernel3D.SetDataDomain(DataDomain);
    p3mKernel3D.SetWorkingSpacingRatio(0.10);
    p3mKernel3D.SetPaddingFactor(3.0 * kernel_width);
    p3mKernel3D.SetSources(Y3D);
    p3mKernel3D.SetKernelWidth(kernel_width);
  }

 protected:

  ExactKernel<ScalarType, 2> exactKernel2D;
  ExactKernel<ScalarType, 3> exactKernel3D;

  P3MKernel<ScalarType, 3> p3mKernel3D;
  P3MKernel<ScalarType, 3> p3mKernel6D;

};

}
}
