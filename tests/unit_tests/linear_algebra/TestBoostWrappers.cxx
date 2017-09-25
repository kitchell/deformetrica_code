/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestBoostWrappers.h"
#include "LinearAlgebra.h"

using namespace def::algebra;

namespace def {
namespace test {

void TestBoostWrappers::SetUp() {
  Test::SetUp();
}

TEST_F(TestBoostWrappers, TestLinearVariableAlgebra) {

  ScalarType sca = 1.0;
  VectorType vec(10, 2.0);
  MatrixType mat(5, 10, 3.0);
  MatrixListType ml(4);
  for (unsigned int k = 0; k < 4; ++k) { ml[k] = (k + 1) * mat; }

  LinearVariableType lv_sca_1 = sca;
  LinearVariableType lv_sca_2 = sca * 5.0;
  LinearVariableType lv_sca_plus = lv_sca_1 + lv_sca_2;
  LinearVariableType lv_sca_div = lv_sca_2 / 2.0;
  ASSERT_EQ(lv_sca_plus.sum_of_squares(), 36);
  ASSERT_EQ(lv_sca_2.sum_of_squares(), 25);
  ASSERT_EQ(lv_sca_div.sum_of_squares(), 6.25);

  LinearVariableType lv_vec_1 = vec;
  LinearVariableType lv_vec_2 = vec * 2.0;
  LinearVariableType lv_vec_plus = lv_vec_1 + lv_vec_2;
  LinearVariableType lv_vec_div = lv_vec_2 / 2.0;
  ASSERT_EQ(lv_vec_plus.sum_of_squares(), 36 * 10);
  ASSERT_EQ(lv_vec_2.sum_of_squares(), 16 * 10);
  ASSERT_EQ(lv_vec_div.sum_of_squares(), 4 * 10);

  LinearVariableType lv_mat_1 = mat;
  LinearVariableType lv_mat_2 = mat * 2.0;
  LinearVariableType lv_mat_plus = lv_mat_1 + lv_mat_2;
  LinearVariableType lv_mat_div = lv_mat_2 / 2.0;
  ASSERT_EQ(lv_mat_plus.sum_of_squares(), 81 * 50);
  ASSERT_EQ(lv_mat_2.sum_of_squares(), 36 * 50);
  ASSERT_EQ(lv_mat_div.sum_of_squares(), 9 * 50);

  LinearVariableType lv_ml_1 = ml;
  LinearVariableType lv_ml_2 = ml * 2.0;
  LinearVariableType lv_ml_plus = lv_ml_1 + lv_ml_2;
  LinearVariableType lv_ml_div = lv_ml_2 / 2.0;
  ASSERT_EQ(lv_ml_plus.sum_of_squares(), (9 * 9 + 18 * 18 + 27 * 27 + 36 * 36) * 50);
  ASSERT_EQ(lv_ml_2.sum_of_squares(), (36 + 144 + 324 + 576) * 50);
  ASSERT_EQ(lv_ml_div.sum_of_squares(), (9 + 36 + 81 + 144) * 50);

}

TEST_F(TestBoostWrappers, TestLinearVariablesAlgebra) {

  ScalarType sca = 1.0;
  VectorType vec(10, 2.0);
  MatrixType mat(5, 10, 3.0);
  MatrixListType ml(4);
  for (unsigned int k = 0; k < 4; ++k) { ml[k] = (k + 1) * mat; }

  std::vector<ScalarType> scav(2);
  scav[0] = sca;
  scav[1] = sca * 2;
  LinearVariablesType lv_sca = scav;
  LinearVariablesType lv_sca_plus = lv_sca + lv_sca;
  LinearVariablesType lv_sca_mult = lv_sca * 2.0;
  LinearVariablesType lv_sca_div = lv_sca_mult / 2.0;
  ASSERT_EQ(lv_sca.sum_of_squares(), 5);
  ASSERT_EQ(lv_sca_plus.sum_of_squares(), lv_sca_mult.sum_of_squares());
  ASSERT_EQ(lv_sca_div.sum_of_squares(), lv_sca.sum_of_squares());

  std::vector<VectorType> vecv(2);
  vecv[0] = vec;
  vecv[1] = vec * 2;
  LinearVariablesType lv_vec = vecv;
  LinearVariablesType lv_vec_plus = lv_vec + lv_vec;
  LinearVariablesType lv_vec_mult = lv_vec * 2.0;
  LinearVariablesType lv_vec_div = lv_vec_mult / 2.0;
  ASSERT_EQ(lv_vec.sum_of_squares(), 4 * 10 * 5);
  ASSERT_EQ(lv_vec_plus.sum_of_squares(), lv_vec_mult.sum_of_squares());
  ASSERT_EQ(lv_vec_div.sum_of_squares(), lv_vec.sum_of_squares());

  std::vector<MatrixType> matv(2);
  matv[0] = mat;
  matv[1] = mat * 2;
  LinearVariablesType lv_mat = matv;
  LinearVariablesType lv_mat_plus = lv_mat + lv_mat;
  LinearVariablesType lv_mat_mult = lv_mat * 2.0;
  LinearVariablesType lv_mat_div = lv_mat_mult / 2.0;
  ASSERT_EQ(lv_mat.sum_of_squares(), 9 * 50 * 5);
  ASSERT_EQ(lv_mat_plus.sum_of_squares(), lv_mat_mult.sum_of_squares());
  ASSERT_EQ(lv_mat_div.sum_of_squares(), lv_mat.sum_of_squares());

  std::vector<MatrixListType> mlv(2);
  mlv[0] = ml;
  mlv[1] = ml * 2;
  LinearVariablesType lv_ml = mlv;
  LinearVariablesType lv_ml_plus = lv_ml + lv_ml;
  LinearVariablesType lv_ml_mult = lv_ml * 2.0;
  LinearVariablesType lv_ml_div = lv_ml_mult / 2.0;
  ASSERT_EQ(lv_ml.sum_of_squares(), (9 + 36 + 81 + 144) * 50 * 5);
  ASSERT_EQ(lv_ml_plus.sum_of_squares(), lv_ml_mult.sum_of_squares());
  ASSERT_EQ(lv_ml_div.sum_of_squares(), lv_ml.sum_of_squares());

}

TEST_F(TestBoostWrappers, TestLinearVariableMapAlgebra) {

  ScalarType sca = 1.0;
  VectorType vec(10, 2.0);
  MatrixType mat(5, 10, 3.0);
  MatrixListType ml(4);
  for (unsigned int k = 0; k < 4; ++k) { ml[k] = (k + 1) * mat; }
  LinearVariableMapType lvm;

  lvm["Scalar"] = sca;
  ASSERT_EQ(lvm.sum_of_squares(), 1);

  lvm["Vector"] = vec;
  ASSERT_EQ(lvm.sum_of_squares(), 1 + 4 * 10);

  lvm["Matrix"] = mat;
  ASSERT_EQ(lvm.sum_of_squares(), 1 + 4 * 10 + 9 * 50);

  lvm["MatrixList"] = ml;
  ASSERT_EQ(lvm.sum_of_squares(), 1 + 4 * 10 + 9 * 50 + (9 + 36 + 81 + 144) * 50);

  LinearVariableMapType lvm_plus = lvm + lvm;
  LinearVariableMapType lvm_mult = lvm * 2.0;
  LinearVariableMapType lvm_div = lvm_mult / 2.0;
  ASSERT_EQ(lvm_plus.sum_of_squares(), lvm_mult.sum_of_squares());
  ASSERT_EQ(lvm_div.sum_of_squares(), lvm.sum_of_squares());

  lvm += lvm;
  ASSERT_EQ(lvm.sum_of_squares(), lvm_plus.sum_of_squares());
  lvm -= lvm_div;
  ASSERT_EQ(lvm.sum_of_squares(), lvm_div.sum_of_squares());
  lvm *= 2;
  ASSERT_EQ(lvm.sum_of_squares(), lvm_mult.sum_of_squares());

}

TEST_F(TestBoostWrappers, TestLinearVariablesMapAlgebra) {

  ScalarType sca = 1.0;
  VectorType vec(10, 2.0);
  MatrixType mat(5, 10, 3.0);
  MatrixListType ml(4);
  for (unsigned int k = 0; k < 4; ++k) { ml[k] = (k + 1) * mat; }
  LinearVariablesMapType lvm;

  std::vector<ScalarType> scav(2);
  scav[0] = sca;
  scav[1] = sca * 2;
  LinearVariablesType lv_sca = scav;

  std::vector<VectorType> vecv(2);
  vecv[0] = vec;
  vecv[1] = vec * 2;
  LinearVariablesType lv_vec = vecv;

  std::vector<MatrixType> matv(2);
  matv[0] = mat;
  matv[1] = mat * 2;
  LinearVariablesType lv_mat = matv;

  std::vector<MatrixListType> mlv(2);
  mlv[0] = ml;
  mlv[1] = ml * 2;
  LinearVariablesType lv_ml = mlv;

  lvm["Scalar"] = lv_sca;
  ASSERT_EQ(lvm["Scalar"].sum_of_squares(), 1 * 5);
  ASSERT_EQ(lvm.sum_of_squares(), 1 * 5);

  lvm["Vector"] = lv_vec;
  ASSERT_EQ(lvm["Vector"].sum_of_squares(), 4 * 10 * 5);
  ASSERT_EQ(lvm.sum_of_squares(), (1 + 4 * 10) * 5);

  lvm["Matrix"] = lv_mat;
  ASSERT_EQ(lvm["Matrix"].sum_of_squares(), 9 * 50 * 5);
  ASSERT_EQ(lvm.sum_of_squares(), (1 + 4 * 10 + 9 * 50) * 5);

  lvm["MatrixList"] = lv_ml;
  ASSERT_EQ(lvm["MatrixList"].sum_of_squares(), (9 + 36 + 81 + 144) * 50 * 5);
  ASSERT_EQ(lvm.sum_of_squares(), (1 + 4 * 10 + 9 * 50 + (9 + 36 + 81 + 144) * 50) * 5);

  LinearVariablesMapType lvm_plus = lvm + lvm;
  LinearVariablesMapType lvm_mult = lvm * 2.0;
  LinearVariablesMapType lvm_div = lvm_mult / 2.0;
  ASSERT_EQ(lvm_plus.sum_of_squares(), lvm_mult.sum_of_squares());
  ASSERT_EQ(lvm_div.sum_of_squares(), lvm.sum_of_squares());

  lvm += lvm;
  ASSERT_EQ(lvm.sum_of_squares(), lvm_plus.sum_of_squares());
  lvm -= lvm_div;
  ASSERT_EQ(lvm.sum_of_squares(), lvm_div.sum_of_squares());
  lvm *= 2;
  ASSERT_EQ(lvm.sum_of_squares(), lvm_mult.sum_of_squares());
}

TEST_F(TestBoostWrappers, TestLinearVariableVectorization) {

  ScalarType sca = 1.0;
  VectorType vec(10, 2.0);
  MatrixType mat(5, 10, 3.0);
  MatrixListType ml(4);
  for (unsigned int k = 0; k < 4; ++k) { ml[k] = (k + 1) * mat; }

  LinearVariableType lv_sca = sca;
  LinearVariableType lv_vec = vec;
  LinearVariableType lv_mat = mat;
  LinearVariableType lv_ml = ml;

  std::vector<unsigned int> struct_sca;
  VectorType lv_sca_vec = lv_sca.vectorize(struct_sca);
  ASSERT_EQ(lv_sca_vec, VectorType(1, sca));
  ASSERT_EQ(unvectorize(lv_sca_vec, struct_sca), lv_sca);
  ASSERT_EQ(unvectorize(lv_sca_vec, struct_sca).vectorize(), lv_sca_vec);

  std::vector<unsigned int> struct_vec;
  VectorType lv_vec_vec = lv_vec.vectorize(struct_vec);
  ASSERT_EQ(lv_vec_vec, vec);
  ASSERT_EQ(unvectorize(lv_vec_vec, struct_vec), lv_vec);
  ASSERT_EQ(unvectorize(lv_vec_vec, struct_vec).vectorize(), lv_vec_vec);

  std::vector<unsigned int> struct_mat;
  VectorType lv_mat_vec = lv_mat.vectorize(struct_mat);
  ASSERT_EQ(lv_mat_vec, mat.vectorize());
  ASSERT_EQ(unvectorize(lv_mat_vec, struct_mat), lv_mat);
  ASSERT_EQ(unvectorize(lv_mat_vec, struct_mat).vectorize(), lv_mat_vec);

  std::vector<unsigned int> struct_ml;
  VectorType lv_ml_vec = lv_ml.vectorize(struct_ml);
  ASSERT_EQ(lv_ml_vec, ml.vectorize());
  ASSERT_EQ(unvectorize(lv_ml_vec, struct_ml), lv_ml);
  ASSERT_EQ(unvectorize(lv_ml_vec, struct_ml).vectorize(), lv_ml_vec);
}

}
}

