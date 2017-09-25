/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "SSDImage.h"

#include "KernelFactory.h"

#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "GridFunctions.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
SSDImage<ScalarType, Dimension>
::SSDImage() : Superclass()
{
  this->SetSSDImageType();
}



template <class ScalarType, unsigned int Dimension>
SSDImage<ScalarType, Dimension>
::SSDImage(const SSDImage& o) : Superclass(o)
{

}

template <class ScalarType, unsigned int Dimension>
SSDImage<ScalarType, Dimension>
::SSDImage(const SSDImage& ex, const MatrixType& IP) : Superclass(ex, IP)
{
  this->SetSSDImageType();
  this->Update();
}


template <class ScalarType, unsigned int Dimension>
SSDImage<ScalarType, Dimension>
::~SSDImage()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, unsigned int Dimension>
void
SSDImage<ScalarType, Dimension>
::Update()
{
  if (this->IsModified())
  {
    Superclass::Update();
  }

  this->UnSetModified();
}

template<class ScalarType, unsigned int Dimension>
ScalarType
SSDImage<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<AbstractGeometryType> target)
{
//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->GetType() != target->GetType())
    throw std::runtime_error("Abstract Geometries types mismatched");

  const std::shared_ptr<const SSDImage> targ = std::static_pointer_cast<const SSDImage>(target);

  this->Update();
  // target->Update();

  typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
  VectorType I0 = GridFunctionsType::VectorizeImage(this->GetImage());
  VectorType I1 = GridFunctionsType::VectorizeImage(targ->GetImage());

  // SSD norm between images
  if (I0.size() != I1.size())
    throw std::runtime_error("image sizes mismatch");

  VectorType D = I0 - I1;
  ScalarType match = D.squared_magnitude();

//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//	auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
//	std::cout << ">>>> SSDImage::ComputeMatch() method took " << d << " ms." << std::endl;

  return match;
}



// compute 2 (I_0(y(0)) - I_1) nabla_{y(0)}I_0
template<class ScalarType, unsigned int Dimension>
MatrixType
SSDImage<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target)
{
//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->GetType() != target->GetType())
    std::cerr << "Abstract Geometries types mismatched: " << this->GetType() << " and " << target->GetType() << "\n";

  const std::shared_ptr<const SSDImage> targ = std::static_pointer_cast<const SSDImage>(target);

  this->Update();
  // target->Update();

  typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;

//	MatrixType Yfinal = GridFunctionsType::UpsampleImagePoints(this->GetImage(), this->GetDownSampledWorkingImage(), Superclass::m_DownSampledY1);

  VectorType I0 = GridFunctionsType::VectorizeImage(this->GetImage());
  VectorType I1 = GridFunctionsType::VectorizeImage(targ->GetImage());
  VectorType D = I0 - I1;

  MatrixType gradMatch(Superclass::m_NumberOfVoxels, Dimension);
  for (unsigned int dim = 0; dim < Dimension; dim++)
  {
    VectorType gradI_d = GridFunctionsType::VectorizeImage(Superclass::m_GradientImages[dim]);
    for (unsigned int i = 0; i < Superclass::m_NumberOfVoxels; i++)
      gradI_d[i] *= D[i];

    gradI_d *= Superclass::m_FlipAxes[dim];
    gradMatch.set_column(Superclass::m_PermutationAxes[dim], gradI_d);
  }


  gradMatch *= 2.0f;


  return gradMatch;
}

template class SSDImage<double,2>;
template class SSDImage<double,3>;

