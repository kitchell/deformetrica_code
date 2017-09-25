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

/// Input-output files.
#include "DeformableObjectReader.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkPolyDataNormals.h"

/// Support files.
#include <src/support/linear_algebra/LinearAlgebra.h>


int main(int argc, char **argv) {
  if (argc != 4) {
      std::cerr << "Usage: " << argv[0]
                << "distance (L2, MutualInformation, EQLA or LCC) pathToImage1 pathToImage2"
                << std::endl;
      return -1;
  }

#ifdef USE_DOUBLE_PRECISION
  typedef double ScalarType;
#else
  typedef float ScalarType;
#endif
  std::string distType = std::string(argv[1]);
  char *pathToImg1 = argv[2];
  char *pathToImg2 = argv[3];

//  std::cout << distType << pathToImg1 << pathToImg2 << std::endl;

  DeformableObjectReader<double, 2> objectReader = DeformableObjectReader<double, 2>();

  DeformableObjectParameters::Pointer parameters = DeformableObjectParameters::New();

  if (distType == "L2")
    parameters->SetDeformableObjectType("SSDImage");
  else if (distType == "LCC")
  {
    parameters->SetDeformableObjectType("LCCImage");
    parameters->SetKernelWidth(5.);
  }
  else if (distType == "EQLA")
  {
    parameters->SetDeformableObjectType("EQLAImage");
    parameters->SetKernelWidth(5.);
  }
  else
    parameters->SetDeformableObjectType("MutualInformationImage");

  objectReader.SetObjectParameters(parameters);

  parameters->SetFilename(pathToImg1);
  objectReader.SetFileName(pathToImg1);
  objectReader.Update();
  std::shared_ptr<AbstractGeometry<double,2>> img1 = objectReader.GetOutput();

  parameters->SetFilename(pathToImg2);
  objectReader.SetFileName(pathToImg2);
  objectReader.Update();
  std::shared_ptr<AbstractGeometry<double,2>> img2 = objectReader.GetOutput();

  ScalarType match = img1->ComputeMatch(img2);

  std::cout << match << std::endl;

  return EXIT_SUCCESS;
}


