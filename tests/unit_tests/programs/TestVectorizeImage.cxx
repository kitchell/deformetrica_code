#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "DeformableObject.h"
#include "Landmark.h"
#include "LinearInterpImage.h"
#include "SSDImage.h"

#include "GridFunctions.h"

#include <cmath>
#include <string>



const unsigned int Dimension = 3;

typedef float ScalarType;

typedef vnl_vector<ScalarType> VNLVectorType;
typedef vnl_matrix<ScalarType> VNLMatrixType;



VNLVectorType vectorizeImage(const char* sourceFileName);



using namespace std;



int main(int argc, char** argv)
{
	char* sourceFileName1 = NULL;
	char* sourceFileName2 = NULL;
	VNLVectorType I1, I2;

	if (argc != 3) {
		cerr << "Usage :" << endl;
		cerr << argv[0] << " sourceFileName1 sourceFileName2" << endl;
		return -1;
	}

	sourceFileName1 = argv[1];
	sourceFileName2 = argv[2];

	I1 = vectorizeImage(sourceFileName1);
	I2 = vectorizeImage(sourceFileName2);

	VNLVectorType delta = I1 - I2;
	cout << "||I1 - I2||^2 = " << delta.squared_magnitude() << endl;


	return 0;
}



VNLVectorType vectorizeImage(const char* sourceFileName) {

	typedef itk::Image<ScalarType, Dimension> ImageType;
	typedef typename ImageType::Pointer ImageTypePointer;

	// Reading source file :
	typedef LinearInterpImage<ScalarType, Dimension> LinearInterpImageType;
	typedef SSDImage<ScalarType, Dimension> SSDImageType;
	LinearInterpImageType* sourceImage = new SSDImageType();
	{
		typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
		typedef itk::ImageFileReader<ImageType> ReaderType;

		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(sourceFileName);
		reader->Update();
		ImageTypePointer img = reader->GetOutput();

		typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
		typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
		minMaxFilter->SetInput(img);
		minMaxFilter->Update();
		double minI = minMaxFilter->GetMinimum();
		double maxI = minMaxFilter->GetMaximum();

		sourceImage->SetImage(img);
		sourceImage->RescaleImage(0.0, 1.0);
		sourceImage->SetMinIntensityOutput(minI);
		sourceImage->SetMaxIntensityOutput(maxI);

		typename ImageType::PointType minAux = img->GetOrigin();
		typename ImageType::PointType maxAux = minAux;
		typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
		typename ImageType::SpacingType spacing = img->GetSpacing();

		// cout << "Image read: origin = " << img->GetOrigin() << " size = " << size << " spacing = " << spacing << " min value = " << minI << " max value = " << maxI << ", rescaled between 0 and 1" << endl;
		// cout << "orientation =\n" << img->GetDirection() << endl;

		itk::Matrix<double,Dimension,Dimension> dir = img->GetDirection();
		vector<int> permutation(Dimension);
		vector<int> flip_axes(Dimension);
		for (unsigned int d=0; d < Dimension; d++)
		{
			int ind = 0;
			double val = fabs(dir(0,d));
			for (unsigned int k=1; k < Dimension; k++)
			{
				if (fabs(dir(k,d))>val)
				{
					ind = k;
					val = fabs(dir(k,d));
				}
			}
			permutation[d] = ind;
			flip_axes[d] = (dir(ind,d)>0.0)?1:-1;
		}
		sourceImage->SetPermutationAxes(permutation);
		sourceImage->SetFlipAxes(flip_axes);

		// cout << "permutation\n" << permutation[0] << permutation[1] << permutation[2] << endl;
		// cout << "flip axes\n" << flip_axes[0] << flip_axes[1] << flip_axes[2] << endl;

		ScalarType ImageGridDownsampling = 1.0; //m_ParamObject->GetImageGridDownsampling();
		ImageTypePointer ExampleDownsampledImage;
		if (ImageGridDownsampling <= 1.0)
			ExampleDownsampledImage = img;
		else
			ExampleDownsampledImage = GridFunctionsType::DownsampleImage(img, ImageGridDownsampling);

		// vnl_matrix<ScalarType> DownsampledY1 = GridFunctionsType::ImageToPoints(ExampleDownsampledImage);

		// std::cout << "DownsampledImage = origin " <<  ExampleDownsampledImage->GetOrigin() << " size = " << ExampleDownsampledImage->GetLargestPossibleRegion().GetSize() << " spacing = " << ExampleDownsampledImage->GetSpacing() << std::endl;
		sourceImage->SetDownSampledWorkingImage(ExampleDownsampledImage);

	}

	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
	return GridFunctionsType::VectorizeImage(sourceImage->GetImage());

}
