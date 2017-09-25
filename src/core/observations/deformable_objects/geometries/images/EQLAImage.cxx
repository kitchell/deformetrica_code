/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "EQLAImage.h"

#include "KernelFactory.h"

#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkSqrtImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkAddImageFilter.h"

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
EQLAImage<ScalarType, Dimension>
::EQLAImage() : Superclass()
{
	this->SetEQLAImageType();
	m_EQLAKernelWidth = 5.0;
	// m_SmallestImageDimension = 0;
}



template <class ScalarType, unsigned int Dimension>
EQLAImage<ScalarType, Dimension>
::EQLAImage(const EQLAImage& o) : Superclass(o)
{

	m_EQLAKernelWidth = o.m_EQLAKernelWidth;
	m_LocalMeanImage = o.m_LocalMeanImage;
	m_LocalVarianceImage = o.m_LocalVarianceImage;
	// m_SmallestImageDimension = o.m_SmallestImageDimension;

}


template <class ScalarType, unsigned int Dimension>
EQLAImage<ScalarType, Dimension>
::EQLAImage(const EQLAImage& ex, const MatrixType& IP) : Superclass(ex, IP)
{
	this->SetEQLAImageType();
	m_EQLAKernelWidth = ex.m_EQLAKernelWidth;
	this->Update();
}



template <class ScalarType, unsigned int Dimension>
EQLAImage<ScalarType, Dimension>
::~EQLAImage()
{
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation methods
////////////////////////////////////////////////////////////////////////////////////////////////////


template <class ScalarType, unsigned int Dimension>
void
EQLAImage<ScalarType, Dimension>
::Update()
{
	if (this->IsModified())
	{
		Superclass::Update();
		
		m_LocalMeanImage = this->ComputeLocalMeanImage(this->GetImage());
		ImageTypePointer LocalVarianceImage = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, this->GetImage(), m_LocalMeanImage);

		// small perturbation to avoid "division by zero" issue in image regions of constant intensity	
		typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddImageFilterType;
		typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
		addFilter->SetInput1(LocalVarianceImage);
		addFilter->SetConstant2(pow(10,-12));
		addFilter->Update();

		m_LocalVarianceImage = addFilter->GetOutput();
	}
	
	this->UnSetModified();
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
ScalarType EQLAImage<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<AbstractGeometryType> target)
{
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	const std::shared_ptr<const EQLAImage> targ = std::static_pointer_cast<const EQLAImage>(target);
	
	this->Update();

	ImageTypePointer LocalCovariance = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, targ->GetImage(), targ->GetLocalMeanImage());

	// Vectorize images
	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
	VectorType LC = GridFunctionsType::VectorizeImage(LocalCovariance);
	VectorType Var1 = GridFunctionsType::VectorizeImage(m_LocalVarianceImage);
	VectorType Var2 = GridFunctionsType::VectorizeImage(targ->GetLocalVarianceImage());
	
	ScalarType match = 0.0;
	for (int i = 0; i < LC.size(); i++)
	{
		match += Var1.get(i) - ( LC.get(i) * LC.get(i) / Var2.get(i) );
	}
	

	return match;
}



template<class ScalarType, unsigned int Dimension>
MatrixType EQLAImage<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target)
{
	if (this->GetType() != target->GetType())
		std::cerr << "Deformable objects types mismatched: " << this->GetType() << " and " << target->GetType() << "\n";

	const std::shared_ptr<const EQLAImage> targ = std::static_pointer_cast<const EQLAImage>(target);
	
	this->Update();
	// target->Update();
	
	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;

	ImageTypePointer LocalCovariance = this->ComputeLocalCovarianceImage(this->GetImage(), m_LocalMeanImage, targ->GetImage(), targ->GetLocalMeanImage());

	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyImageFilterType;
	typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractImageFilterType;
	typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DivideImageFilterType;
	
	// Compute Cov(Idef,targ)/Var(targ)
	typename DivideImageFilterType::Pointer divideFilter = DivideImageFilterType::New();
	divideFilter->SetInput1(LocalCovariance);
	divideFilter->SetInput2(targ->GetLocalVarianceImage());
	divideFilter->Update();
	
	// Compute W*( Cov(Idef,targ)/Var(targ) )
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetSigma( m_EQLAKernelWidth );
	gaussf->SetInput(divideFilter->GetOutput());
	gaussf->Update();
	
	ImageTypePointer Convolution1 = gaussf->GetOutput();
	
	// Compute W*( targMean * Cov(Idef,targ)/Var(targ) )
	typename MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
	multiplyFilter->SetInput1(targ->GetLocalMeanImage());
	multiplyFilter->SetInput2(divideFilter->GetOutput());
	multiplyFilter->Update();

	typename GaussianFilterType::Pointer gaussf2 = GaussianFilterType::New();
	gaussf2->SetSigma( m_EQLAKernelWidth );
	gaussf2->SetInput(multiplyFilter->GetOutput());
	gaussf2->Update();
	
	ImageTypePointer Convolution2 = gaussf2->GetOutput();
	
	// Compute W*( IdefMean)
	typename GaussianFilterType::Pointer gaussf3 = GaussianFilterType::New();
	gaussf3->SetSigma( m_EQLAKernelWidth );
	gaussf3->SetInput(m_LocalMeanImage);
	gaussf3->Update();

	ImageTypePointer IdefMeanMean = gaussf3->GetOutput();
	
	
	// Vectorize required images to compute the gradient
	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
	VectorType VtargImg = GridFunctionsType::VectorizeImage(targ->GetImage());
	VectorType VdefImg = GridFunctionsType::VectorizeImage(this->GetImage());
	VectorType VdefImgMeanMean = GridFunctionsType::VectorizeImage(IdefMeanMean);
	VectorType VConv1 = GridFunctionsType::VectorizeImage(Convolution1);
	VectorType VConv2 = GridFunctionsType::VectorizeImage(Convolution2);
	
	// Compute the gradient of the deformed source image
	MatrixType gradMatch(Superclass::m_NumberOfVoxels, Dimension, 0.0);
	MatrixType gradI(Superclass::m_NumberOfVoxels, Dimension);

	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		VectorType gradI_d = GridFunctionsType::VectorizeImage(Superclass::m_GradientImages[dim]);
		gradI_d *= Superclass::m_FlipAxes[dim];
		gradI.set_column(Superclass::m_PermutationAxes[dim], gradI_d);
	}

	
	for (int i = 0; i < Superclass::m_NumberOfVoxels; i++)
	{
		ScalarType val = VdefImg.get(i) - VdefImgMeanMean.get(i) - VtargImg.get(i) * VConv1.get(i) + VConv2.get(i);
		gradMatch.set_row(i, val * gradI.get_row(i) );
	}

	gradMatch *= 2.0f;
	return gradMatch;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class ScalarType, unsigned int Dimension>
typename EQLAImage<ScalarType, Dimension>::ImageTypePointer
EQLAImage<ScalarType, Dimension>
::ComputeLocalMeanImage(const ImageType* img) const
{
	// convolution between input image and a gaussian filter: W*I
	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetInput(img);
	gaussf->SetSigma( m_EQLAKernelWidth );
	gaussf->Update();

	return gaussf->GetOutput();
}



template<class ScalarType, unsigned int Dimension>
typename EQLAImage<ScalarType, Dimension>::ImageTypePointer
EQLAImage<ScalarType, Dimension>
::ComputeLocalCovarianceImage(const ImageType* img1, const ImageType* LocalMeanImg1,
		const ImageType* img2, const ImageType* LocalMeanImg2) const
{
	// multiply img1 and img2 voxel-wise in multiplyFilter->GetOutput()
	typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyImageFilterType;
	typename MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
	multiplyFilter->SetInput1(img1);
	multiplyFilter->SetInput2(img2);
	multiplyFilter->Update();

	// convolve the multiplied image: W*(I.J) in gaussf->GetOutput()
	typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
	typename GaussianFilterType::Pointer gaussf = GaussianFilterType::New();
	gaussf->SetSigma( m_EQLAKernelWidth );
	gaussf->SetInput(multiplyFilter->GetOutput());
	gaussf->Update();

	// multiply the local mean images (W*I).(W*J) in multiplyFilter->GetOutput()
	typename MultiplyImageFilterType::Pointer multiplyFilter2 = MultiplyImageFilterType::New();
	multiplyFilter2->SetInput1(LocalMeanImg1);
	multiplyFilter2->SetInput2(LocalMeanImg2);
	multiplyFilter2->Update();
	
	// subtract the last two images: W*(I.J) - (W*I).(W*J)
	typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractImageFilterType;
	typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
	subtractFilter->SetInput1(gaussf->GetOutput());
	subtractFilter->SetInput2(multiplyFilter2->GetOutput());
	subtractFilter->Update();
	
	return subtractFilter->GetOutput();

}

template class EQLAImage<double,2>;
template class EQLAImage<double,3>;

