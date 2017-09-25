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
#include <random>
#include "LinearAlgebra.h"

using namespace def::algebra;

namespace def {
namespace test {

class AbstractTestKernelPrecision : public ::testing::Test {
 public:
  // Constructor : initialize data and kernel used in the tests
  AbstractTestKernelPrecision() {
    kernel_width = 1.88;
    X_ROW = 501;
    Y_ROW = 201;
    W_ROW = Y_ROW;
    Z_ROW = X_ROW;

    X2D = MatrixType(X_ROW, 2);
    generate_random_matrix(X2D);
    Y2D = MatrixType(Y_ROW, 2);
    generate_random_matrix(Y2D);
    W2D = MatrixType(W_ROW, 2);
    generate_random_matrix(W2D);
    Z2D = MatrixType(Z_ROW, 2);
    generate_random_matrix(Z2D);

    X3D = MatrixType(X_ROW, 3);
    generate_random_matrix(X3D);
    Y3D = MatrixType(Y_ROW, 3);
    generate_random_matrix(Y3D);
    W3D = MatrixType(W_ROW, 3);
    generate_random_matrix(W3D);
    Z3D = MatrixType(Z_ROW, 3);
    generate_random_matrix(Z3D);

    W4D = MatrixType(W_ROW, 4);
    generate_random_matrix(W4D);
    W6D = MatrixType(W_ROW, 6);
    generate_random_matrix(W6D);

  }

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif

  void ASSERT_MATRIX_EQ(const MatrixType &X, const MatrixType &Y, const std::string &msg, const ScalarType error = 0.0);
  void EXPECT_MATRIX_EQ(const MatrixType &X, const MatrixType &Y, const std::string &msg, const ScalarType error = 0.0);
  void generate_random_matrix(MatrixType &X);
  void CompareAndDisp(MatrixType &X, MatrixType &Y, ScalarType precision, std::string msg, std::string msg2);

 protected:
  virtual void SetUp();

  // shared data initialized in constructor
  size_t X_ROW, Y_ROW, W_ROW, Z_ROW;
  ScalarType kernel_width;

  MatrixType X2D, Y2D, W2D, Z2D, X3D, Y3D, W3D, Z3D, W4D, W6D;

#ifdef USE_DOUBLE_PRECISION
  ScalarType eps_tol = 1e-10;
  ScalarType fuzzy_tol = 8e-2;
#else
  ScalarType eps_tol = 1e-6;
  ScalarType fuzzy_tol = 15e-2;
#endif
};

}
}
