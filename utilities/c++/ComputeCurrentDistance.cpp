/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "DeformetricaConfig.h"

#include <src/core/observations/deformable_objects/geometries/landmarks/OrientedPolyLine.h>
#include <src/core/observations/deformable_objects/geometries/landmarks/OrientedSurfaceMesh.h>

#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <src/support/linear_algebra/LinearAlgebra.h>
#include <src/support/kernels/KernelFactory.h>
#include <src/support/kernels/KernelType.h>



int main(int argc, char **argv) {

  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " dimension kernelWidth pathToMesh1 pathToMesh2"
              << std::endl;
    return -1;
  }

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif

  unsigned int dimension = atoi(argv[1]);

  ScalarType kernelWidth = atof(argv[2]);
  assert(kernelWidth>0);

  char *pathToMesh1 = argv[3];
  char *pathToMesh2 = argv[4];

  /// Typedefs.
  typedef OrientedPolyLine<ScalarType, 2> OrientedPolyLineType;
  typedef OrientedSurfaceMesh<ScalarType, 3> OrientedSurfaceMeshType;

  /// Main.
  if (dimension == 2) {

    /// Reading the poly line data.
    std::shared_ptr<OrientedPolyLineType> pl1 = std::make_shared<OrientedPolyLineType>();
    vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
    reader1->SetFileName(pathToMesh1);
    reader1->Update();
    vtkSmartPointer<vtkPolyData> polyData1 = reader1->GetOutput();
    pl1->SetAnatomicalCoordinateSystem("LPS");
    pl1->SetPolyData(polyData1);
    pl1->SetKernelType(Exact);
    pl1->SetKernelWidth(kernelWidth);
    pl1->Update();

    std::shared_ptr<OrientedPolyLineType> pl2 = std::make_shared<OrientedPolyLineType>();
    vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(pathToMesh2);
    reader2->Update();
    vtkSmartPointer<vtkPolyData> polyData2 = reader2->GetOutput();
    pl2->SetAnatomicalCoordinateSystem("LPS");
    pl2->SetPolyData(polyData2);
    pl2->SetKernelType(Exact);
    pl2->SetKernelWidth(kernelWidth);
    pl2->Update();

    /// Compute and display match value.
    ScalarType match = pl1->ComputeMatch(pl2);

    cout.precision(10);
    std::cout << match << std::endl;

  } else if (dimension == 3) {

    /// Reading the surface meshes.
    std::shared_ptr<OrientedSurfaceMeshType> sm1 = std::make_shared<OrientedSurfaceMeshType>();
    vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
    reader1->SetFileName(pathToMesh1);
    reader1->Update();
    vtkSmartPointer<vtkPolyData> polyData1 = reader1->GetOutput();
    sm1->SetAnatomicalCoordinateSystem("LPS");
    sm1->SetPolyData(polyData1);
    sm1->SetReorient();
    sm1->SetKernelType(Exact);
    sm1->SetKernelWidth(kernelWidth);
    sm1->Update();

    std::shared_ptr<OrientedSurfaceMeshType> sm2 = std::make_shared<OrientedSurfaceMeshType>();
    vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(pathToMesh2);
    reader2->Update();
    vtkSmartPointer<vtkPolyData> polyData2 = reader2->GetOutput();
    sm2->SetAnatomicalCoordinateSystem("LPS");
    sm2->SetPolyData(polyData2);
    sm2->SetReorient();
    sm2->SetKernelType(Exact);
    sm2->SetKernelWidth(kernelWidth);
    sm2->Update();

    /// Compute and display match value.
    ScalarType match = sm1->ComputeMatch(sm2);

    cout.precision(10);
    std::cout << match << std::endl;

  } else { std::cerr << "Dimension must be either 2 or 3." << std::endl; }

  return EXIT_SUCCESS;
}


