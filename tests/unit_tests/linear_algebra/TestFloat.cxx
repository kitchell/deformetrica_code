/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestFloat.h"

namespace def {
    namespace test {

        void TestFloat::SetUp() {
          Test::SetUp();
        }

        TEST_F(TestFloat, TestNorm) {
        VectorType vec(10, 0.0);
        vec(0) = 1.0;
        vec(1) = 2.0;
        vec(3) = 3.0;

        auto var = vec.squared_magnitude();
        float f;
        double d;

        std::cout << "var = " << typeid(var).name() << std::endl;
        std::cout << "float = " << typeid(f).name() << std::endl;
        std::cout << "double = " << typeid(d).name() << std::endl;

        ASSERT_EQ(var, (ScalarType)14.0);
    }

        TEST_F(TestFloat, TestSumOfSquares) {
            VectorType vec(10, 0.0);
            vec(0) = 1.0;
            vec(1) = 2.0;
            vec(2) = 3.0;

            auto var = vec.sum_of_squares();
            float f;
            double d;

            std::cout << "var = " << typeid(var).name() << std::endl;
            std::cout << "float = " << typeid(f).name() << std::endl;
            std::cout << "double = " << typeid(d).name() << std::endl;

            ASSERT_EQ(var, (ScalarType)14.0);
        }

        TEST_F(TestFloat, TestMagnitude) {
            VectorType vec(10, 0.0);
            vec(0) = 3.0;
            vec(3) = 4.0;

            auto var = vec.magnitude();
            float f;
            double d;

            std::cout << "var = " << typeid(var).name() << std::endl;
            std::cout << "float = " << typeid(f).name() << std::endl;
            std::cout << "double = " << typeid(d).name() << std::endl;

            ASSERT_EQ(var, (ScalarType)5.0);
        }

        TEST_F(TestFloat, TestArmaNorm) {
            VectorType vec(10, 0.0);
            vec(0) = 0.0;
            vec(1) = 3.0;
            vec(2) = 4.0;

            auto var = arma::norm(vec.toArmadillo(),2);
            float f;
            double d;

            std::cout << "var = " << typeid(var).name() << std::endl;
            std::cout << "float = " << typeid(f).name() << std::endl;
            std::cout << "double = " << typeid(d).name() << std::endl;

            ASSERT_EQ(var, (ScalarType)5.0);
        }

//        TEST_F(TestFloat, TestNormWithInputVectorFloat) {
//            VectorType vec = readMatrixDLM<float>(
//                UNIT_TESTS_DIR"/linear_algebra/data/testVector.txt").get_column(0);
//
//            auto var = vec.squared_magnitude();
//            float f;
//            double d;
//
//            std::cout << "var = " << typeid(var).name() << std::endl;
//            std::cout << "float = " << typeid(f).name() << std::endl;
//            std::cout << "double = " << typeid(d).name() << std::endl;
//
//            ASSERT_LE(var, (float)90.0);
//            ASSERT_GE(var, (float)70.0);
//        }
//
//        TEST_F(TestFloat, TestNormWithInputVectorDouble) {
//            VectorType vec = readMatrixDLM<double>(
//                UNIT_TESTS_DIR"/linear_algebra/data/testVector.txt").get_column(0);
//
//            auto var = vec.squared_magnitude();
//            float f;
//            double d;
//
//            std::cout << "var = " << typeid(var).name() << std::endl;
//            std::cout << "float = " << typeid(f).name() << std::endl;
//            std::cout << "double = " << typeid(d).name() << std::endl;
//
//            ASSERT_LE(var, (double)90.0);
//            ASSERT_GE(var, (double)70.0);
//        }

        TEST_F(TestFloat, TestNormWithInputVector) {
            VectorType vec = readMatrixDLM<ScalarType>(
                UNIT_TESTS_DIR"/linear_algebra/data/testVector.txt").get_column(0);

            auto var = vec.squared_magnitude();
            float f;
            double d;

            std::cout << "var = " << typeid(var).name() << std::endl;
            std::cout << "float = " << typeid(f).name() << std::endl;
            std::cout << "double = " << typeid(d).name() << std::endl;

            ASSERT_LE(var, (ScalarType)90.0);
            ASSERT_GE(var, (ScalarType)70.0);
        }


}
}
