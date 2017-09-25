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

#include "myvtkPolyDataNormals.h"

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

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

  /// Reorienting.
  vtkSmartPointer<myvtkPolyDataNormals> normalf = vtkSmartPointer<myvtkPolyDataNormals>::New();
  normalf->SetInputData(input);
  normalf->ConsistencyOn();
  normalf->AutoOrientNormalsOn(); // Should have closed surface
  normalf->FlipNormalsOn();
  normalf->Update();

  vtkSmartPointer<vtkPolyData> output = normalf->GetOutput();
  output->BuildLinks();

  /// Saving result.
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(outFilename);

#if (VTK_MAJOR_VERSION <= 5)
  writer->SetInput(output);
#elseif (VTK_MAJOR_VERSION < 6)
  writer->SetInputConnection(output);
#else
  writer->SetInputData(output);
#endif

  writer->Update();

  return EXIT_SUCCESS;
}