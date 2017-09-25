/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformetricaConfig.h"

/// Input-output files.
#include <src/io/MatrixDLM.h>

/// Support files.
#include <src/support/linear_algebra/LinearAlgebra.h>
#include <src/support/kernels/KernelFactory.h>
#include <src/support/kernels/ExactKernel.h>

/// Libraries files.
#include <iostream>
#include <cstring>

using namespace def::algebra;


template<class ScalarType, unsigned int Dimension>
void main_aux(ScalarType const &kernelWidth,
              MatrixType const &controlPoints1,
              MatrixType const &momenta1,
              MatrixType const &controlPoints2,
              MatrixType const &momenta2) {

  /// Number of control points.
  const unsigned int nbControlPoints1 = controlPoints1.rows();
  const unsigned int nbControlPoints2 = controlPoints2.rows();

  /// Initialize the kernel.
  ExactKernel<ScalarType, Dimension> kernel;
  kernel.SetKernelWidth(kernelWidth);

  ScalarType norm1;
  kernel.SetSources(controlPoints1);
  kernel.SetWeights(momenta1);
  MatrixType Velocity = kernel.Convolve(controlPoints1);
  norm1 = std::sqrt(dot_product(Velocity, Velocity));



  MatrixType allControlPoints(nbControlPoints1+nbControlPoints2,Dimension,0.), allMomentas(nbControlPoints1+nbControlPoints2,Dimension,0.);
  for (unsigned int i=0;i<nbControlPoints1;++i)
  {
    allControlPoints.set_row(i, controlPoints1.get_row(i));
    allMomentas.set_row(i, momenta1.get_row(i));
  }
  for (unsigned int i=0;i<nbControlPoints2;++i)
  {
    allControlPoints.set_row(i+nbControlPoints1, controlPoints2.get_row(i));
    allMomentas.set_row(i+nbControlPoints1, -1.*momenta2.get_row(i));
  }

  kernel.SetSources(allControlPoints);
  kernel.SetWeights(allMomentas);
  Velocity = kernel.Convolve(allControlPoints);

  ScalarType relativeError = std::sqrt(dot_product(Velocity, Velocity))/norm1;

  std::cout << relativeError << std::endl;

}


int main(int argc, char **argv) {

  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " dimension kernelWidth pathToControlPoints1 pathToMomenta1 pathToControlPoints2 pathToMomenta2"
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
  char *pathToControlPoints1 = argv[3];
  char *pathToMomenta1 = argv[4];
  char *pathToControlPoints2 = argv[5];
  char *pathToMomenta2 = argv[6];

  /// Read input files.
  const MatrixType controlPoints1 = readMatrixDLM<ScalarType>(pathToControlPoints1);
  const MatrixType momenta1 = readMatrixDLM<ScalarType>(pathToMomenta1);
  const MatrixType controlPoints2 = readMatrixDLM<ScalarType>(pathToControlPoints2);
  const MatrixType momenta2 = readMatrixDLM<ScalarType>(pathToMomenta2);

  /// Dimension.
  assert(controlPoints1.cols() == dimension);
  assert(controlPoints2.cols() == dimension);

  /// Launch.
  if (dimension == 2) { main_aux<ScalarType, 2>(kernelWidth, controlPoints1, momenta1, controlPoints2, momenta2); }
  else if (dimension == 3) { main_aux<ScalarType, 3>(kernelWidth, controlPoints1, momenta1, controlPoints2, momenta2); }

}


