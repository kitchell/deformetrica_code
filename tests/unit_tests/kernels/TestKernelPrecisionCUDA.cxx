/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestKernelPrecisionCUDA.h"
#include <time.h>

namespace def {
namespace test {

// Convolve2D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_Convolve_2) {

exactKernel2D.SetWeights(W2D);
cudaExactKernel2D.SetWeights(W2D);

MatrixType result_made_by_exact_kernel = exactKernel2D.Convolve(X2D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel2D.Convolve(X2D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// Convolve4D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_Convolve_4) {

exactKernel2D.SetWeights(W4D);
cudaExactKernel2D.SetWeights(W4D);

MatrixType result_made_by_exact_kernel = exactKernel2D.Convolve(X2D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel2D.Convolve(X2D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// ConvolveGradient2D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradient_2) {

exactKernel2D.SetWeights(W2D);
cudaExactKernel2D.SetWeights(W2D);

MatrixType result_made_by_exact_kernel = exactKernel2D.ConvolveGradient(X2D, Z2D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel2D.ConvolveGradient(X2D, Z2D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// ConvolveGradient2D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradientVarlin_4) {

exactKernel2D.SetWeights(W4D);
cudaExactKernel2D.SetWeights(W4D);

MatrixListType result_made_by_exact_kernel = exactKernel2D.ConvolveGradient(X2D);
MatrixListType result_made_by_cuda_kernel = cudaExactKernel2D.ConvolveGradient(X2D);

for (int j = 0; j < result_made_by_exact_kernel.size(); j++)
CompareAndDisp(result_made_by_exact_kernel[j], result_made_by_cuda_kernel[j], eps_tol, "exact", "cuda");
}

// ConvolveGradient2D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradientVarlin_2) {

exactKernel2D.SetWeights(W2D);
cudaExactKernel2D.SetWeights(W2D);

MatrixListType result_made_by_exact_kernel = exactKernel2D.ConvolveGradient(X2D);
MatrixListType result_made_by_cuda_kernel = cudaExactKernel2D.ConvolveGradient(X2D);

for (int j = 0; j < result_made_by_exact_kernel.size(); j++)
CompareAndDisp(result_made_by_exact_kernel[j], result_made_by_cuda_kernel[j], eps_tol, "exact", "cuda");
}

// ConvolveGradient2D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradient1_2) {

exactKernel2D.SetWeights(W2D);
cudaExactKernel2D.SetWeights(W2D);

MatrixType result_made_by_exact_kernel0 = exactKernel2D.ConvolveGradient(Y2D, int(0));
MatrixType result_made_by_cuda_kernel0 = cudaExactKernel2D.ConvolveGradient(Y2D, int(0));

CompareAndDisp(result_made_by_exact_kernel0, result_made_by_cuda_kernel0, eps_tol, "exact", "cuda");

MatrixType result_made_by_exact_kernel1 = exactKernel2D.ConvolveGradient(Y2D, int(1));
MatrixType result_made_by_cuda_kernel1 = cudaExactKernel2D.ConvolveGradient(Y2D, int(1));

CompareAndDisp(result_made_by_exact_kernel1, result_made_by_cuda_kernel1, eps_tol, "exact", "cuda");
}

// ConvolveHessian2D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveSpecialHessian_2) {

exactKernel2D.SetWeights(W2D);
cudaExactKernel2D.SetWeights(W2D);

MatrixType result_made_by_exact_kernel = exactKernel2D.ConvolveSpecialHessian(W2D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel2D.ConvolveSpecialHessian(W2D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// Convolve3D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_Convolve_3) {

exactKernel3D.SetWeights(W3D);
cudaExactKernel3D.SetWeights(W3D);

MatrixType result_made_by_exact_kernel = exactKernel3D.Convolve(X3D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel3D.Convolve(X3D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// Convolve6D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_Convolve_6) {

exactKernel3D.SetWeights(W6D);
cudaExactKernel3D.SetWeights(W6D);

MatrixType result_made_by_exact_kernel = exactKernel3D.Convolve(X3D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel3D.Convolve(X3D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// ConvolveGradient3D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradientVarlin_6) {

exactKernel3D.SetWeights(W6D);
cudaExactKernel3D.SetWeights(W6D);

MatrixListType result_made_by_exact_kernel = exactKernel3D.ConvolveGradient(X3D);
MatrixListType result_made_by_cuda_kernel = cudaExactKernel3D.ConvolveGradient(X3D);

for (int j = 0; j < result_made_by_exact_kernel.size(); j++)
CompareAndDisp(result_made_by_exact_kernel[j], result_made_by_cuda_kernel[j], eps_tol, "exact", "cuda");
}

// ConvolveGradient3D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradientVarlin_3) {

exactKernel3D.SetWeights(W3D);
cudaExactKernel3D.SetWeights(W3D);

MatrixListType result_made_by_exact_kernel = exactKernel3D.ConvolveGradient(X3D);
MatrixListType result_made_by_cuda_kernel = cudaExactKernel3D.ConvolveGradient(X3D);

for (int j = 0; j < result_made_by_exact_kernel.size(); j++)
CompareAndDisp(result_made_by_exact_kernel[j], result_made_by_cuda_kernel[j], eps_tol, "exact", "cuda");
}

// ConvolveGradient3D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradient_3) {

exactKernel3D.SetWeights(W3D);
cudaExactKernel3D.SetWeights(W3D);

MatrixType result_made_by_exact_kernel = exactKernel3D.ConvolveGradient(X3D, Z3D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel3D.ConvolveGradient(X3D, Z3D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");
}

// ConvolveGradient3D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveGradient1_3) {

exactKernel3D.SetWeights(W3D);
cudaExactKernel3D.SetWeights(W3D);

MatrixType result_made_by_exact_kernel0 = exactKernel3D.ConvolveGradient(Y3D, int(0));
MatrixType result_made_by_cuda_kernel0 = cudaExactKernel3D.ConvolveGradient(Y3D, int(0));

CompareAndDisp(result_made_by_exact_kernel0, result_made_by_cuda_kernel0, eps_tol, "exact", "cuda");

MatrixType result_made_by_exact_kernel1 = exactKernel3D.ConvolveGradient(Y3D, int(1));
MatrixType result_made_by_cuda_kernel1 = cudaExactKernel3D.ConvolveGradient(Y3D, int(1));

CompareAndDisp(result_made_by_exact_kernel1, result_made_by_cuda_kernel1, eps_tol, "exact", "cuda");

MatrixType result_made_by_exact_kernel2 = exactKernel3D.ConvolveGradient(Y3D, int(2));
MatrixType result_made_by_cuda_kernel2 = cudaExactKernel3D.ConvolveGradient(Y3D, int(2));

CompareAndDisp(result_made_by_exact_kernel2, result_made_by_cuda_kernel2, eps_tol, "exact", "cuda");
}

// ConvolveHessian3D
TEST_F(TestKernelPrecisionCUDA, exact_vs_cuda_ConvolveSpecialHessian_3) {

exactKernel3D.SetWeights(W3D);
cudaExactKernel3D.SetWeights(W3D);

MatrixType result_made_by_exact_kernel = exactKernel3D.ConvolveSpecialHessian(W3D);
MatrixType result_made_by_cuda_kernel = cudaExactKernel3D.ConvolveSpecialHessian(W3D);

CompareAndDisp(result_made_by_exact_kernel, result_made_by_cuda_kernel, eps_tol, "exact", "cuda");

}

}
}

