/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestReadParametersXML.h"
#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"

namespace def {
namespace test {

void TestReadParametersXML::SetUp() {
    Test::SetUp();
}

TEST_F(TestReadParametersXML, SparseDiffeoParameters) {
    SparseDiffeoParameters::Pointer paramDiffeos;

    paramDiffeos = readSparseDiffeoParametersXML("fake_file");
    ASSERT_TRUE(paramDiffeos.IsNull());

    paramDiffeos = readSparseDiffeoParametersXML(UNIT_TESTS_DIR"/io/data/paramDiffeos.xml");
    ASSERT_FALSE(paramDiffeos.IsNull());

    ASSERT_EQ(paramDiffeos->GetKernelWidth(), 1.5);
    ASSERT_EQ(paramDiffeos->GetKernelType(), "Exact");
    ASSERT_EQ(paramDiffeos->GetInitialCPSpacing(), 0.3);
    ASSERT_EQ(paramDiffeos->GetNumberOfTimePoints(), 20);
    ASSERT_EQ(paramDiffeos->GetMaxIterations(), 200);
    ASSERT_EQ(paramDiffeos->GetMaxLineSearchIterations(), 20);
    ASSERT_EQ(paramDiffeos->GetStepExpand(), 2.0);
    ASSERT_EQ(paramDiffeos->GetStepShrink(), 0.5);
    ASSERT_EQ(paramDiffeos->GetInitialStepSize(), 0.001);
    ASSERT_EQ(paramDiffeos->GetAdaptiveTolerance(), 1e-4);
    ASSERT_EQ(paramDiffeos->GetNumberOfThreads(), 1);
    ASSERT_EQ(paramDiffeos->GetModelType(), "bayesian");
}

TEST_F(TestReadParametersXML, DeformableObjectParameters) {
    DeformableObjectParameters::Pointer paramObjects;

    paramObjects = readDeformableObjectParametersXML(UNIT_TESTS_DIR"/io/data/paramSurface.xml");
    ASSERT_FALSE(paramObjects.IsNull());

    ASSERT_EQ(paramObjects->GetDeformableObjectType(), "OrientedSurfaceMesh");
    ASSERT_EQ(paramObjects->GetDataSigma(), 0.01);
    ASSERT_EQ(paramObjects->GetKernelWidth(), 2.5);
    ASSERT_EQ(paramObjects->GetKernelType(), "Exact");
}
}
}
