//
// Created by LOUIS Maxime on 08/11/2016.
//

///Unit test for the parallel transport method of the diffeos class.


#include "TestParallelTransport.h"
#include <src/core/model_tools/deformations/Diffeos.h>
#include <src/support/kernels/KernelType.h>
#include <src/io/MatrixDLM.h>
#include <src/io/DeformableObjectReader.h>
#include <src/core/observations/deformable_objects/geometries/images/SSDImage.h>

namespace def {
    namespace test {

        void TestParallelTransport::SetUp() {
            Test::SetUp();
        }

        TEST_F(TestParallelTransport, CHECK_PARALLEL_TRANSPORT) {
            ///We initialize a diffeo object, feed it with some input, and assert the result is correct

            ///Filenames with input and ground truth output
            std::string initialCPRegressionFile = UNIT_TESTS_DIR"/parallel-transport/data/CPRegression.txt";
            std::string initialMomentaFile = UNIT_TESTS_DIR"/parallel-transport/data/MomRegression.txt";
            std::string initialCPToTransportFile = UNIT_TESTS_DIR"/parallel-transport/data/CPToTransport.txt";
            std::string initialMomentaToTransportFile = UNIT_TESTS_DIR"/parallel-transport/data/MomToTransport.txt";
            std::string groundTruthLastMomFile = UNIT_TESTS_DIR"/parallel-transport/data/GroundTruthMom.txt";
            std::string imageFile = UNIT_TESTS_DIR"/parallel-transport/data/Image.png";

            ///Reading data from the files
            MatrixType CPRegression = readMatrixDLM<double>(initialCPRegressionFile.c_str());
            MatrixType MomRegression = readMatrixDLM<double>(initialMomentaFile.c_str());
            MatrixType CPMatching = readMatrixDLM<double>(initialCPToTransportFile.c_str());
            MatrixType MomMatching = readMultipleMatrixDLM<double>(initialMomentaToTransportFile.c_str())[0];
            MatrixType groundTruthMom = readMatrixDLM<double>(groundTruthLastMomFile.c_str());

            ///Reading a deformable object (used only so that the Diffeo does not complain)
            std::shared_ptr<SSDImage<double,2>> objectImg = std::make_shared<SSDImage<double,2>>();
            ImageTypePointer img;
            typename ImageReaderType ::Pointer imgReader = ReaderType::New();
            imgReader->SetFileName(imageFile.c_str());
            img = imgReader->GetOutput();
            objectImg->SetImage(img);
            objectImg->Update();
            AbstractGeometryList objectList(1.);
            objectList[0] = objectImg;
            std::shared_ptr<DeformableMultiObject<double,2>> multiObject = std::make_shared<DeformableMultiObject<double,2>>();
            multiObject->SetObjectList(objectList);

            ///Instatiating the diffeo object.
            std::shared_ptr<Diffeos<double,2>> def;
            def = std::make_shared<Diffeos<double,2>>();
            def->SetKernelWidth(15.);
            def->SetKernelType(Exact);
            def->SetNumberOfTimePoints(30);
            def->SetStartMomentas(MomRegression);
            def->SetStartPositions(CPRegression);
            def->SetDeformableMultiObject(multiObject);
            MatrixType boundingBox = objectImg->GetBoundingBox();
            def->SetDataDomain(boundingBox);
            def->UseImprovedEuler();
            def->Update();

            std::vector<double> targetTimes(10);
            for (unsigned int i = 0; i < 10; ++i)
              targetTimes[i] = i * 1. / 9.;

            ///Computing parallel transport.
            MatrixListType velocities;
            MatrixListType transportedMomentas = def->ParallelTransport(MomMatching, CPMatching, 0., targetTimes, velocities);
            MatrixType lastMomentas = transportedMomentas.at(transportedMomentas.size()-1);


            for (int j = 0;j<lastMomentas.rows();j++)
            {
                for (int k=0;k<lastMomentas.cols();k++) {
                  ASSERT_LE(std::abs(lastMomentas(j, k) - groundTruthMom(j, k)), 1e-4);
                }
            }
        }
    }
}



