/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestKernelPrecisionP3M.h"
#include <time.h>

namespace def {
namespace test {



// Convolve3D
TEST_F(TestKernelPrecisionP3M, p3m_vs_exact_Convolve_3) {

  p3mKernel3D.SetWeights(W3D);
  exactKernel3D.SetWeights(W3D);

  MatrixType result_made_by_p3m_kernel = p3mKernel3D.Convolve(X3D);
  MatrixType result_made_by_exact_kernel = exactKernel3D.Convolve(X3D);

  CompareAndDisp(result_made_by_exact_kernel, result_made_by_p3m_kernel, fuzzy_tol, "p3m", "exact");

}

// Convolve6D
TEST_F(TestKernelPrecisionP3M, p3m_vs_exact_Convolve_6) {

  p3mKernel3D.SetWeights(W6D);
  exactKernel3D.SetWeights(W6D);

  MatrixType result_made_by_p3m_kernel = p3mKernel3D.Convolve(X3D);
  MatrixType result_made_by_exact_kernel = exactKernel3D.Convolve(X3D);

  CompareAndDisp(result_made_by_exact_kernel, result_made_by_p3m_kernel, fuzzy_tol, "p3m", "exact");

}

// ConvolveGradient3D
TEST_F(TestKernelPrecisionP3M, p3m_vs_exact_ConvolveGradient_3) {

  p3mKernel3D.SetWeights(W3D);
  exactKernel3D.SetWeights(W3D);

  MatrixType result_made_by_p3m_kernel = p3mKernel3D.ConvolveGradient(X3D, Z3D);
  MatrixType result_made_by_exact_kernel = exactKernel3D.ConvolveGradient(X3D, Z3D);

  CompareAndDisp(result_made_by_exact_kernel, result_made_by_p3m_kernel, fuzzy_tol, "p3m", "exact");

}

// ConvolveGradient3D
TEST_F(TestKernelPrecisionP3M, p3m_vs_exact_ConvolveGradient1_3) {

  p3mKernel3D.SetWeights(W3D);
  exactKernel3D.SetWeights(W3D);

  MatrixType result_made_by_p3m_kernel0 = p3mKernel3D.ConvolveGradient(Y3D, int(0));
  MatrixType result_made_by_exact_kernel0 = exactKernel3D.ConvolveGradient(Y3D, int(0));

  CompareAndDisp(result_made_by_exact_kernel0, result_made_by_p3m_kernel0, fuzzy_tol, "p3m", "exact");

  MatrixType result_made_by_p3m_kernel1 = p3mKernel3D.ConvolveGradient(Y3D, int(1));
  MatrixType result_made_by_exact_kernel1 = exactKernel3D.ConvolveGradient(Y3D, int(1));

  CompareAndDisp(result_made_by_exact_kernel1, result_made_by_p3m_kernel1, fuzzy_tol, "p3m", "exact");

  MatrixType result_made_by_p3m_kernel2 = p3mKernel3D.ConvolveGradient(Y3D, int(2));
  MatrixType result_made_by_exact_kernel2 = exactKernel3D.ConvolveGradient(Y3D, int(2));

  CompareAndDisp(result_made_by_exact_kernel2, result_made_by_p3m_kernel2, fuzzy_tol, "p3m", "exact");
}

// ConvolveHessian3D
TEST_F(TestKernelPrecisionP3M, p3m_vs_exact_ConvolveSpecialHessian_3) {

  p3mKernel3D.SetWeights(W3D);
  exactKernel3D.SetWeights(W3D);

  MatrixType result_made_by_p3m_kernel = p3mKernel3D.ConvolveSpecialHessian(W3D);
  MatrixType result_made_by_exact_kernel = exactKernel3D.ConvolveSpecialHessian(W3D);

  CompareAndDisp(result_made_by_exact_kernel, result_made_by_p3m_kernel, fuzzy_tol, "p3m", "exact");

}

}
}