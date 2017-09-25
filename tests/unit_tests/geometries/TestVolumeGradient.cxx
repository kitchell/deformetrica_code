
#include "TestVolumeGradient.h"

using namespace def::algebra;

namespace def {
    namespace test {

        void TestVolumeGradient::SetUp() {
            Test::SetUp();
        }

        template <unsigned int Dimension> std::shared_ptr<OrientedVolumeMesh<ScalarType, Dimension>> TestVolumeGradient::ReadOrientedVolumeMesh(const char *filePath) {
            std::shared_ptr<OrientedVolumeMesh<ScalarType, Dimension>> out = std::make_shared<OrientedVolumeMesh<ScalarType, Dimension>>();
            if (Dimension==3) {
              vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
              reader->SetFileName(filePath);
              reader->Update();
              auto * unstructuredData = reader->GetOutput();
              out->SetPolyData((vtkPointSet *) unstructuredData);
            }
            else {
              vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
              reader->SetFileName(filePath);
              reader->Update();
              auto *unstructuredData = reader->GetOutput();
              out->SetPolyData((vtkPointSet *) unstructuredData);
            }
            out->SetAnatomicalCoordinateSystem("LPS");
            out->SetKernelType(KernelEnumType::Exact);
            out->SetKernelWidth(1.0);
            return out;
        }

        MatrixType
        TestVolumeGradient::GetNumericalGradient3(MatrixType &sourcePoints, std::shared_ptr<OrientedVolumeMesh<ScalarType, 3>> source,
                                           const std::shared_ptr<OrientedVolumeMesh<ScalarType, 3>> target) {
            ScalarType h = 0.00001f;
            ScalarType match = source->ComputeMatch(target);
            MatrixType numericallyEstimatedMatchGradient(sourcePoints.rows(), sourcePoints.cols());
            //For all the points :
            for (int i = 0; i < sourcePoints.rows(); i++)
                //For all the coordinates :
                for (int d = 0; d < 3; d++) {
                    //Add h to the d-th coord of the i-th points
                    ScalarType aux = sourcePoints(i, d);
                    sourcePoints(i, d) = sourcePoints(i, d) + h;
                    //Create a deformed object.
                    OrientedVolumeMesh<ScalarType, 3> slightlyDeformedSource = OrientedVolumeMesh<ScalarType, 3>(
                            *source, sourcePoints);
                    slightlyDeformedSource.Update();
                    //f' = ((f(x+h) - f(x))/h
                    numericallyEstimatedMatchGradient(i, d) = (slightlyDeformedSource.ComputeMatch(target) - match) / h;
                    sourcePoints(i, d) = aux;
                }
            return numericallyEstimatedMatchGradient;
        }

        MatrixType
        TestVolumeGradient::GetNumericalGradient2(MatrixType &sourcePoints, std::shared_ptr<OrientedVolumeMesh<ScalarType, 2>> source,
                                           const std::shared_ptr<OrientedVolumeMesh<ScalarType, 2>> target) {
            ScalarType h = 0.00001f;
            ScalarType match = source->ComputeMatch(target);
            MatrixType numericallyEstimatedMatchGradient(sourcePoints.rows(), sourcePoints.cols());
            //For all the points :
            for (int i = 0; i < sourcePoints.rows(); i++)
                //For all the coordinates :
                for (int d = 0; d < 2; d++) {
                    //Add h to the d-th coord of the i-th points
                    ScalarType aux = sourcePoints(i, d);
                    sourcePoints(i, d) = sourcePoints(i, d) + h;
                    //Create a deformed object.
                    OrientedVolumeMesh<ScalarType, 2> slightlyDeformedSource = OrientedVolumeMesh<ScalarType, 2>(
                            *source, sourcePoints);
                    slightlyDeformedSource.Update();
                    //f' = ((f(x+h) - f(x))/h
                    numericallyEstimatedMatchGradient(i, d) = (slightlyDeformedSource.ComputeMatch(target) - match) / h;
                    sourcePoints(i, d) = aux;
                }
            return numericallyEstimatedMatchGradient;
        }

        TEST_F(TestVolumeGradient, CHECK_GRADIENT_VOLUME_MESH) {
          ///Dimension 3
            std::shared_ptr<OrientedVolumeMesh<ScalarType, 3>> source3;
            source3 = ReadOrientedVolumeMesh<3>(UNIT_TESTS_DIR"/geometries/data/VolumeCube.vtk");
            std::shared_ptr<OrientedVolumeMesh<ScalarType, 3>> target3 = ReadOrientedVolumeMesh<3>(
                    UNIT_TESTS_DIR"/geometries/data/VolumeBall.vtk");
            source3->Update();
            target3->Update();
            MatrixType sourcePoints = source3->GetPointCoordinates();
            MatrixType analyticalGradient = source3->ComputeMatchGradient(target3);
            MatrixType numericalGradient = GetNumericalGradient3(sourcePoints, source3, target3);
            for (int i = 0; i < analyticalGradient.rows(); i++)
                for (int j = 0; j < analyticalGradient.cols(); j++) {
                  ASSERT_LE(std::abs(analyticalGradient(i, j) - numericalGradient(i, j)),1e-5);
            }

            ///Dimension 2
            std::shared_ptr<OrientedVolumeMesh<ScalarType, 2>> source2;
            source2 = ReadOrientedVolumeMesh<2>(UNIT_TESTS_DIR"/geometries/data/VolumeSquare.vtk");
            std::shared_ptr<OrientedVolumeMesh<ScalarType, 2>> target2 = ReadOrientedVolumeMesh<2>(
                    UNIT_TESTS_DIR"/geometries/data/VolumeDisk.vtk");
            source2->Update();
            target2->Update();
            sourcePoints = source2->GetPointCoordinates();
            analyticalGradient = source2->ComputeMatchGradient(target2);
            numericalGradient = GetNumericalGradient2(sourcePoints, source2, target2);
            for (int i = 0; i < analyticalGradient.rows(); i++)
                for (int j = 0; j < analyticalGradient.cols(); j++) {
                    ASSERT_LE(std::abs(analyticalGradient(i,j)- numericalGradient(i,j)),1e-5);
                }
        }
}
}



