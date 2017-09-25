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
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " sourceFilename outFilename" << endl;
    return -1;
  }

  char *sourceFilename = argv[1];
  char *outFilename = argv[2];

  unsigned int numberOfIterations = 0;

  /// Reading source file.
  vtkSmartPointer<vtkPolyDataReader> inputReader = vtkSmartPointer<vtkPolyDataReader>::New();
  inputReader->SetFileName(sourceFilename);
  inputReader->Update();
  vtkSmartPointer<vtkPolyData> input = inputReader->GetOutput();

  // Smoothing :
  vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  smoothFilter->SetInputData(input);
  smoothFilter->SetNumberOfIterations(numberOfIterations);
  smoothFilter->Update();

  /// Saving result.
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(outFilename);

#if (VTK_MAJOR_VERSION <= 5)
  writer->SetInput(smoothFilter->GetOutput());
#elseif (VTK_MAJOR_VERSION < 6)
  writer->SetInputConnection(input->GetProducerPort());
#else
  writer->SetInputData(smoothFilter->GetOutput());
#endif

  writer->Update();

  return EXIT_SUCCESS;
}