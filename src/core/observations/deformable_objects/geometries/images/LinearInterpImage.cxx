/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/


#include "LinearInterpImage.h"

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
#include "itkImageDuplicator.h"

#include "itkOrientImageFilter.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
LinearInterpImage<ScalarType, Dimension>
::LinearInterpImage() : Superclass()
{
	m_MinIntensityOutput = 0;
	m_MaxIntensityOutput = (Dimension==2)?255:1.0;

	std::vector<unsigned int> p(Dimension);
	std::vector<int> f(Dimension);
	for (unsigned int d=0; d<Dimension; d++)
	{
		p[d] = d;
		f[d] = 1;
	}
	m_PermutationAxes = p;
	m_FlipAxes = f;
	
}



template <class ScalarType, unsigned int Dimension>
LinearInterpImage<ScalarType, Dimension>
::LinearInterpImage(const LinearInterpImage& o) : Superclass(o)
{
    typedef typename itk::ImageDuplicator< ImageType > DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

	if (!o.m_Image.IsNull()) {
	    duplicator->SetInputImage(o.m_Image);
	    duplicator->Update();
	    m_Image = duplicator->GetOutput(); }

	if (!o.m_DownSampledImage.IsNull()) {
		duplicator->SetInputImage(o.m_DownSampledImage);
		duplicator->Update();
		m_DownSampledImage = duplicator->GetOutput(); }
	
	m_GradientImages.resize(Dimension);
	for (int dim = 0; dim < Dimension; dim++)
	{
		if (!o.m_GradientImages[dim].IsNull())
		{
			duplicator->SetInputImage(o.m_GradientImages[dim]);
			duplicator->Update();
			m_GradientImages[dim] = duplicator->GetOutput();
		}
	}
	
	m_MinIntensityOutput = o.m_MinIntensityOutput;
	m_MaxIntensityOutput = o.m_MaxIntensityOutput;

	m_PermutationAxes = o.m_PermutationAxes;
	m_FlipAxes = o.m_FlipAxes;
	
	m_NumberOfVoxels = o.m_NumberOfVoxels;
	
}


template <class ScalarType, unsigned int Dimension>
LinearInterpImage<ScalarType, Dimension>
::
LinearInterpImage(const LinearInterpImage& example, const MatrixType& DownSampledImageMap) : Superclass(example, DownSampledImageMap)
{
	if (example.m_Image.IsNull() || example.m_DownSampledImage.IsNull())
		throw std::runtime_error("LinearInterpImage : Cannot re-sample image if no image has been set.");
	
    typedef typename itk::ImageDuplicator< ImageType > DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
	
    duplicator->SetInputImage(example.m_Image);
    duplicator->Update();
    m_Image = duplicator->GetOutput();		
	
	duplicator->SetInputImage(example.m_DownSampledImage);
	duplicator->Update();
	m_DownSampledImage = duplicator->GetOutput();
	
	m_GradientImages.resize(Dimension);
	for (int dim = 0; dim < Dimension; dim++)
	{
		duplicator->SetInputImage(example.m_GradientImages[dim]);
		duplicator->Update();
		m_GradientImages[dim] = duplicator->GetOutput();
	}
		
	m_MinIntensityOutput = example.m_MinIntensityOutput;
	m_MaxIntensityOutput = example.m_MaxIntensityOutput;

	m_PermutationAxes = example.m_PermutationAxes;
	m_FlipAxes = example.m_FlipAxes;
	
	m_NumberOfVoxels = example.m_NumberOfVoxels;
	
	// Upsample ImagePoints to be at the same resolution as the reference image
	if (example.m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels() != DownSampledImageMap.rows())
		throw std::runtime_error("Number of voxels in downsampled image maps mismatch in LinearInterpImage constructor");

	/// For bug-tracking.
//	writeMatrixDLM<ScalarType>("0_Y.txt", DownSampledImageMap);

	MatrixType Yup = this->UpSampleImageMap(DownSampledImageMap);
	VectorType I1 = GridFunctionsType::Interpolate(Yup, m_Image);
	
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType it(m_Image, m_Image->GetLargestPossibleRegion());

    unsigned int r = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      it.Set(I1[r++]);
    }

	for (int dim = 0; dim < Dimension; dim++)
	{
		VectorType I1g = GridFunctionsType::Interpolate(Yup, m_GradientImages[dim]);
	
	    IteratorType itg(m_GradientImages[dim], m_GradientImages[dim]->GetLargestPossibleRegion());

	    unsigned int r = 0;
	    for (itg.GoToBegin(); !itg.IsAtEnd(); ++itg)
	    {
	      itg.Set(I1g[r++]);
	    }
		
	}

	// DownsampledReferenceImage is not re-sampled since only the positions of its voxels is used.
	// Call to Update() is done in children classes
}




template <class ScalarType, unsigned int Dimension>
LinearInterpImage<ScalarType, Dimension>
::~LinearInterpImage()
{}


template <class ScalarType, unsigned int Dimension>
void
LinearInterpImage<ScalarType, Dimension>
::SetImageAndDownSamplingFactor(ImageType* img, int ds)
{
	if (ds<0)
		throw std::runtime_error("Negative downsampling factor not allowed");
	
	m_Image = img;
	m_NumberOfVoxels = img->GetLargestPossibleRegion().GetNumberOfPixels();
	
	if (ds==0 || ds == 1)
	{
		m_DownSampledImage = m_Image; // I don't think we need to duplicate
	}
	else
	{
		m_DownSampledImage = GridFunctionsType::DownsampleImage(m_Image, ds);
	}
	
	
	m_GradientImages.resize(Dimension);
	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
		typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
		derivf->SetInput(m_Image);
		derivf->SetDirection(dim);
		derivf->SetOrder(1);
		derivf->SetUseImageSpacingOn();
		derivf->Update();
	
		m_GradientImages[dim] = derivf->GetOutput();
	}
	
	this->SetModified();
}


template <class ScalarType, unsigned int Dimension>
void
LinearInterpImage<ScalarType, Dimension>
::UpdateImageIntensity(const MatrixType& I)
{
	if (m_Image.IsNull())
		throw std::runtime_error("ITK image should have been set before setting new intensities");

	if (I.rows() != m_NumberOfVoxels)
		throw std::runtime_error("image size and number of pixels mismatched");
	
	if (I.cols() != 1)
		throw std::runtime_error("only scalar images allowed");
		
	int i=0;
	
	typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
	ImageIterator it(m_Image, m_Image->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		it.Set( I(i++,0) );

	this->SetModified(); 
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void 
LinearInterpImage<ScalarType, Dimension>
::Update()
{
	if (m_Image.IsNull() || m_DownSampledImage.IsNull())
		throw std::runtime_error("Reference image and downsampled version of it should be set in LinearInterpImage");

	if (this->IsModified())
	{
			
		this->UpdateBoundingBox();

		// // If m_DownSampledY1 has been set, it is used to re-sample the reference image, otherwise one assumes m_DownSampledY1 to be a regular sampling of m_DownSampledReferenceImage (i.e. no deformation)	
		// if (m_DownSampledY1.rows() == 0)
		// {
		// 	m_DownSampledY1.set_size(m_DownSampledReferenceImage->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension);
		// 	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
		// 	m_DownSampledY1 = GridFunctionsType::ImageToPoints(m_DownSampledReferenceImage);
		// }
		// 
		// this->UpdateWorkingImageFromDownSampledImagePoints();
	}
	// this->UnSetModified() not necessary, since this method could only be called by children classes, which will take care of this.
}



template<class ScalarType, unsigned int Dimension>
MatrixType
LinearInterpImage<ScalarType, Dimension>
::SplatDifferenceImage(const std::shared_ptr<LinearInterpImage<ScalarType, Dimension>> target,
					   const MatrixType& DownSampledImageMap)
{
	MatrixType Y0 = UpSampleImageMap(DownSampledImageMap);
	VectorType I0 = GridFunctionsType::Interpolate(Y0, m_Image);
	VectorType I1 = GridFunctionsType::VectorizeImage(target->GetImage());
	VectorType Residual = I0 - I1;

	ImageTypePointer splat = GridFunctionsType::SplatToImage(m_Image, Y0, Residual, 0);
	
	MatrixType out(m_NumberOfVoxels,1);
	
	out.set_column(0 ,GridFunctionsType::VectorizeImage(splat));

	return out;
}




template <class ScalarType, unsigned int Dimension>
void LinearInterpImage<ScalarType, Dimension>
::WriteObject(std::string str) const
{
	// Rescale and convert intensities of the working image, then save it.
	// If initial intensity range is not included in [0,1], intensity values are converted to unsigned short.
	// Otherwise, Tscalar (float or double) is used
	if ( (m_MinIntensityOutput >= 0.0) && (m_MaxIntensityOutput > 1.0) )
	{
		typedef itk::Image<unsigned char, Dimension> OutImageType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
		typename CasterType::Pointer castf = CasterType::New();
		castf->SetInput(m_Image);
		castf->SetOutputMinimum(m_MinIntensityOutput);
		castf->SetOutputMaximum(m_MaxIntensityOutput);
		castf->Update();

		// std::cout << "Defined cast filter" << std::endl;
		// If original image (> 2D) was not oriented according to LPS reference, switch back to original coordinate system
		std::string lps_label("LPS");

		if ((Dimension ==3) && (lps_label.compare(Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel()) != 0))
		{
			typedef itk::Image<unsigned char, 3> OutImageType_dim_3;
			typename itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::Pointer orienter = itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::New();
			orienter->UseImageDirectionOn();
			orienter->SetDesiredCoordinateOrientation(Superclass::m_AnatomicalOrientation.GetITKAnatomicalCoordinateSystemLabel());
			orienter->SetInput(dynamic_cast<OutImageType_dim_3*>( castf->GetOutput() ));
			orienter->Update();

			// std::cout << "Orientation filter" << std::endl;
			typedef itk::ImageFileWriter<OutImageType_dim_3> WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput(orienter->GetOutput()); //
			writer->SetFileName(str.c_str());
			writer->Update();
			// std::cout << "Writing filter for oriented image" << std::endl;
		}
		else
		{
			// write output image
			typedef itk::ImageFileWriter<OutImageType> WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput(castf->GetOutput());
			writer->SetFileName(str.c_str());
			writer->Update();
			// std::cout << "Writing filter" << std::endl;
		}
	}
	else
	{
		typedef itk::Image<ScalarType, Dimension> OutImageType;
		typedef itk::RescaleIntensityImageFilter<ImageType, OutImageType> CasterType;
		typename CasterType::Pointer castf = CasterType::New();
		castf->SetInput(m_Image);
		castf->SetOutputMinimum(m_MinIntensityOutput);
		castf->SetOutputMaximum(m_MaxIntensityOutput);
		castf->Update();

		// If original image was not oriented according to LPS reference, switch back to original coordinate system
		std::string lps_label("LPS");

		if ((Dimension ==3) && (lps_label.compare(Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel()) != 0))
		{
			std::cout << "Reorienting output image from LPS orientation to the original " << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel() << " orientation." << std::endl;

			typedef itk::Image<ScalarType, 3> OutImageType_dim_3;
			typename itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::Pointer orienter = itk::OrientImageFilter<OutImageType_dim_3,OutImageType_dim_3>::New();
			orienter->UseImageDirectionOn();
			orienter->SetDesiredCoordinateOrientation(Superclass::m_AnatomicalOrientation.GetITKAnatomicalCoordinateSystemLabel());
			orienter->SetInput(dynamic_cast<OutImageType_dim_3*>( castf->GetOutput() ));
			orienter->Update();

			typedef itk::ImageFileWriter<OutImageType_dim_3> WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput(orienter->GetOutput()); //
			writer->SetFileName(str.c_str());
			writer->Update();
		}
		else
		{
	    // write output image
		typedef itk::ImageFileWriter<OutImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput(castf->GetOutput()); //
		writer->SetFileName(str.c_str());
		writer->Update();
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void 
LinearInterpImage<ScalarType, Dimension>
::UpdateBoundingBox()
{
	// Not really efficient... but the safiest given possible weird coordinate systems...
	MatrixType Yex = GridFunctionsType::ImageToPoints(m_Image);
	
	for (unsigned int d=0; d<Dimension; d++)
	{
		Superclass::m_BoundingBox(d,0) = Yex.get_column(d).min_value(); // - 0.5;
		Superclass::m_BoundingBox(d,1) = Yex.get_column(d).max_value(); // + 0.5;
	}
}



template <class ScalarType, unsigned int Dimension>
MatrixType
LinearInterpImage<ScalarType, Dimension>
::UpSampleImageMap(const MatrixType& Y)
{
	if (Y.rows() != m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels())
		throw std::runtime_error("Number of points in image map and size of downsampled image mismatch");
	
	return GridFunctionsType::UpsampleImagePoints(m_Image, m_DownSampledImage, Y);
}





template class LinearInterpImage<double,2>;
template class LinearInterpImage<double,3>;
