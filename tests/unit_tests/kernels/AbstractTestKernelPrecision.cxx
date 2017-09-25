/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "AbstractTestKernelPrecision.h"
#include <time.h>

namespace def {
namespace test {

void AbstractTestKernelPrecision::SetUp() {
  Test::SetUp();
}

void AbstractTestKernelPrecision::generate_random_matrix(MatrixType &X) {
  std::mt19937 gen(1);

  for (int i = 0; i < X.rows(); ++i) {
    for (int j = 0; j < X.cols(); ++j) {
      X(i, j) = std::exp(-(ScalarType) gen() / (ScalarType) (gen.max()));
    }
  }
}

void AbstractTestKernelPrecision::ASSERT_MATRIX_EQ(const MatrixType &X,
                                           const MatrixType &Y,
                                           const std::string &msg,
                                           const ScalarType error) {
  ASSERT_EQ(X.rows(), Y.rows()) << msg;
  ASSERT_EQ(X.cols(), Y.cols()) << msg;

  for (int i = 0; i < X.rows(); ++i) {
    for (int j = 0; j < X.cols(); ++j) {
      ScalarType x = X(i, j);
      ScalarType y = Y(i, j);
      ASSERT_LE(std::min(std::abs(x - y) / std::abs(x), std::abs(x - y)), error)
                    << msg << " on position (" << i << "," << j << ")";
    }
  }
}

void AbstractTestKernelPrecision::EXPECT_MATRIX_EQ(const MatrixType &X,
                                           const MatrixType &Y,
                                           const std::string &msg,
                                           const ScalarType error) {
  ASSERT_EQ(X.rows(), Y.rows()) << msg;
  ASSERT_EQ(X.cols(), Y.cols()) << msg;

  for (int i = 0; i < X.rows(); ++i) {
    for (int j = 0; j < X.cols(); ++j) {
      ScalarType x = X(i, j);
      ScalarType y = Y(i, j);
      EXPECT_LE(std::min(std::abs(x - y) / std::abs(x), std::abs(x - y)), error)
                << msg << " on position (" << i << "," << j << ")  " << x << "  " << y;
    }
  }
}

void AbstractTestKernelPrecision::CompareAndDisp(MatrixType &results1,
                                         MatrixType &results2,
                                         ScalarType precision,
                                         std::string msg,
                                         std::string msg2) {
  EXPECT_MATRIX_EQ(results1, results2, msg + " VS " + msg2, precision);
}

}
}
