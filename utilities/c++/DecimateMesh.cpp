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

#include <vtkVersion.h>
#include <vtkPolyData.h>
#include <vtkDecimatePro.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmoothPolyDataFilter.h>

using namespace std;

int main(int argc, char **argv) {

  /// Initialization.
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " sourceFilename outFilename reductionRate" << endl;
    return -1;
  }

  char *sourceFilename= argv[1];
  char *outFilename  = argv[2];
  double reductionRate = atof(argv[3]);
  double featureAngle = 0.0;

  /// Reading source file.
  vtkSmartPointer<vtkPolyDataReader> inputReader = vtkSmartPointer<vtkPolyDataReader>::New();
  inputReader->SetFileName(sourceFilename);
  inputReader->Update();
  vtkSmartPointer<vtkPolyData> input = inputReader->GetOutput();

  /// Decimation.
  vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
#if VTK_MAJOR_VERSION <= 5
  decimate->SetInputConnection(input->GetProducerPort());
#else
  decimate->SetInputData(input);
#endif
  decimate->SetTargetReduction(reductionRate);
  decimate->SetFeatureAngle(featureAngle);
  decimate->Update();
  vtkSmartPointer<vtkPolyData> decimated = vtkSmartPointer<vtkPolyData>::New();
  decimated->ShallowCopy(decimate->GetOutput());
  std::cout << input->GetNumberOfPoints() << " / " << input->GetNumberOfPoints() << " / "
            << decimated->GetNumberOfPoints() << " points (before / intermediate / after)." << std::endl;
  std::cout << input->GetNumberOfPolys() << " / " << input->GetNumberOfPolys() << " / "
            << decimated->GetNumberOfPolys() << " polygons (before / intermediate / after)." << std::endl;

  /// Saving result.
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(outFilename);
#if (VTK_MAJOR_VERSION <= 5)
  writer->SetInput(decimate->GetOutput());
#elseif (VTK_MAJOR_VERSION < 6)
  writer->SetInputConnection(input->GetProducerPort());
#else
  writer->SetInputData(decimate->GetOutput());
#endif
  writer->Update();

  return EXIT_SUCCESS;
}