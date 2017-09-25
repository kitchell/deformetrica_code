/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestSerialization.h"
#include <memory>
#include <cstdio>
#include "src/support/utilities/SerializeDeformationState.h"
#include <stdlib.h>

using namespace def::utils;
using namespace def::algebra::utils;

namespace def {
namespace test {

struct stmpfile {

  std::string file;

  stmpfile(std::string s=UNIT_TESTS_DIR"/serialize/data_empty_dir/ser.bin"):
  file(s)
  {}

  ~stmpfile() {
    std::remove(file.c_str());
  }
};

void TestSerialization::SetUp() {
  Test::SetUp();
}

TEST_F(TestSerialization, StoreAndLoad) {

  bool testBool1 = false;
  bool testBool2 = true;

  /* Define and Init Common use Variables */
  ScalarType testSca = 3.0015;
  VectorType testVec(10, 2.0364);
  MatrixType testMat(10, 15, 1.0014);
  MatrixListType testMatList(5);
  for (ScalarType k = 0 ; k < 5 ; ++k) {
    testMatList[k] = testMat*k;
  }

  LinearVariableType testLinVar = testVec;
  LinearVariableType testLinVar2 = testMat;
  LinearVariableMapType testLinVarMap;
  testLinVarMap["TestScalar"] = testSca;
  testLinVarMap["TestVector"] = testVec;
  testLinVarMap["TestMatrix"] = testMat;
  testLinVarMap["TestMatrixList"] = testMatList;

  LinearVariablesType testLinVars1({testVec, testMat});
  LinearVariablesType testLinVars2({testSca, testVec, testMat});
  LinearVariablesMapType testLinVarsMap;
  testLinVarsMap["var1"] = testLinVars1;
  testLinVarsMap["var2"] = testLinVars2;

  std::vector<ScalarType> testScas(5);
  std::vector<VectorType> testVecs(5);
  std::vector<MatrixType> testMats(5);
  std::vector<MatrixListType> testMatLiss(5);
  for (ScalarType k = 0 ; k < 5 ; ++k)
  {
    testScas[k] = testSca*k;
    testVecs[k] = testVec*k;
    testMats[k] = testMat*k;
    testMatLiss[k] = testMatList*k;
  }

  stmpfile tmpnam = stmpfile();

  /* Store all variables */
  DeformationState state_store;
  {
    state_store << testBool1 << testBool2;
    state_store << testSca << testVec << testMat << testMatList << testLinVar;
    state_store << testLinVarMap << testLinVars1 << testLinVarsMap;
    state_store << testScas << testVecs << testMats << testMatLiss;
    state_store.save(tmpnam.file);
  }

  /* Load from file all variables and check it*/
  DeformationState state_load;
  {
    bool b1,b2;
    ScalarType t;
    VectorType v;
    MatrixType m;


    MatrixListType ml;
    LinearVariableType lv;
    LinearVariableMapType lm;
    LinearVariablesType lvs;
    LinearVariablesMapType lvsm;
    std::vector<ScalarType> v_s;
    std::vector<VectorType> v_v;
    std::vector<MatrixType> v_m;
    std::vector<MatrixListType> v_ml;

    state_load.load(tmpnam.file);
    ASSERT_EQ(state_store, state_load);

    state_load >> b1 >> b2;
    state_load >> t >> v >> m >> ml >> lv >> lm >> lvs >> lvsm;
    state_load >> v_s >> v_v >> v_m >> v_ml;

    ASSERT_EQ(b1, testBool1);
    ASSERT_EQ(b2, testBool2);
    ASSERT_EQ(t, testSca);
    ASSERT_EQ(v, testVec);
    ASSERT_EQ(m, testMat);
    ASSERT_EQ(ml, testMatList);
    /// In order to fix boost 1.58 issue with "visitor"
    ASSERT_TRUE(lv == testLinVar);
    ASSERT_TRUE(lm == testLinVarMap);
    ASSERT_TRUE(lvs == testLinVars1);
    ASSERT_TRUE(lvsm == testLinVarsMap);
    ASSERT_TRUE(v_s == testScas);
    ASSERT_TRUE(v_v == testVecs);
    ASSERT_TRUE(v_m == testMats);
    ASSERT_TRUE(v_ml == testMatLiss);
  }

}

TEST_F(TestSerialization, SingletonUse) {
  stmpfile tmpnam = stmpfile();

  int a = 1;
  float b = 3.14;
  {
    serialize << a << b;
    serialize.save(tmpnam.file);
  }
  {
    int t_a;
    float t_b;
    serialize.load(tmpnam.file);
    serialize >> t_a >> t_b;
    ASSERT_EQ(t_a, a);
    ASSERT_EQ(t_b, b);
  }

  std::string c = "string";

  DeformationState *df = new DeformationState;
  *df << a << b << c;
  {
    serialize.deformation(std::move(*df));
    serialize.save(tmpnam.file);
  }
  {
    int t_a;
    float t_b;
    std::string t_c;
    serialize.load(tmpnam.file);
    serialize >> t_a >> t_b >> t_c;
    ASSERT_EQ(t_a, a);
    ASSERT_EQ(t_b, b);
    ASSERT_EQ(t_c, c);
  }

}

TEST_F(TestSerialization, Tolerance) {

  const float good_tolerance = 1e-5;
  const float bad_tolerance = 1e-9;
  const float add_tolerance = 1e-7;

  /* Define and Init Common use Variables */
  ScalarType testSca = 3.0015;
  ScalarType err_testSca = testSca + add_tolerance;

  VectorType testVec(10, 2.0364);
  VectorType err_testVec(10, 2.0364 + add_tolerance);

  MatrixType testMat(10, 15, 1.0014);
  MatrixType err_testMat(10, 15, 1.0014 + add_tolerance);

  MatrixListType testMatList(5);
  MatrixListType err_testMatList(5);

  for (ScalarType k = 0 ; k < 5 ; ++k) {
    testMatList[k] = testMat*k;
    testMatList[k] = testMat*k;
    err_testMatList[k] = err_testMat*k;
    err_testMatList[k] = err_testMat*k;
  }

  LinearVariableType testLinVar = testVec;
  LinearVariableType err_testLinVar = err_testVec;
  
  LinearVariableType testLinVar2 = testMat;
  LinearVariableType err_testLinVar2 = err_testMat;
  
  LinearVariableMapType testLinVarMap;
  LinearVariableMapType err_testLinVarMap;
  testLinVarMap["TestScalar"] = testSca;
  testLinVarMap["TestVector"] = testVec;
  testLinVarMap["TestMatrix"] = testMat;
  testLinVarMap["TestMatrixList"] = testMatList;
  err_testLinVarMap["TestScalar"] = err_testSca;
  err_testLinVarMap["TestVector"] = err_testVec;
  err_testLinVarMap["TestMatrix"] = err_testMat;
  err_testLinVarMap["TestMatrixList"] = err_testMatList;

  LinearVariablesType testLinVars1({testVec, testMat});
  LinearVariablesType err_testLinVars1({err_testVec, err_testMat});

  LinearVariablesType testLinVars2({testSca, testVec, testMat});
  LinearVariablesType err_testLinVars2({err_testSca, err_testVec, err_testMat});

  LinearVariablesMapType testLinVarsMap;
  LinearVariablesMapType err_testLinVarsMap;
  testLinVarsMap["var1"] = testLinVars1;
  testLinVarsMap["var2"] = testLinVars2;
  err_testLinVarsMap["var1"] = err_testLinVars1;
  err_testLinVarsMap["var2"] = err_testLinVars2;

  std::vector<ScalarType> testScas(5);
  std::vector<ScalarType> err_testScas(5);
  
  std::vector<VectorType> testVecs(5);
  std::vector<VectorType> err_testVecs(5);
  
  std::vector<MatrixType> testMats(5);
  std::vector<MatrixType> err_testMats(5);
  
  std::vector<MatrixListType> testMatLiss(5);
  std::vector<MatrixListType> err_testMatLiss(5);
  for (ScalarType k = 0 ; k < 5 ; ++k)
  {
    testScas[k] = testSca*k;
    testVecs[k] = testVec*k;
    testMats[k] = testMat*k;
    testMatLiss[k] = testMatList*k;
    err_testScas[k] = err_testSca*k;
    err_testVecs[k] = err_testVec*k;
    err_testMats[k] = err_testMat*k;
    err_testMatLiss[k] = err_testMatList*k;
  }

  {
    DeformationState df1,df2;
    df1 << testSca << testVec << testMat << testMatList << testLinVar;
    df2 << err_testSca << err_testVec << err_testMat << err_testMatList << err_testLinVar;

    df1 << testLinVarMap << testLinVars1 << testLinVarsMap;
    df2 << testLinVarMap << testLinVars1 << testLinVarsMap;

    df1 << testScas << testVecs << testMats << testMatLiss;
    df2 << testScas << testVecs << testMats << testMatLiss;

    const bool match = df1 == df2;
    ASSERT_FALSE(match);
    ASSERT_TRUE(df1.compare(df2, good_tolerance));
    ASSERT_FALSE(df1.compare(df2, bad_tolerance));
  }

}



}
}
