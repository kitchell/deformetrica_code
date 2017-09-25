/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformableObjectReader.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
DeformableObjectReader<ScalarType, Dimension>
::DeformableObjectReader() : m_IsTemplate(false) {}

template<class ScalarType, unsigned int Dimension>
DeformableObjectReader<ScalarType, Dimension>
::~DeformableObjectReader() {}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
DeformableObjectReader<ScalarType, Dimension>
::Update() {
  typedef Landmark<ScalarType, Dimension> LandmarkType;
  typedef PointCloud<ScalarType, Dimension> PointCloudType;
  typedef OrientedPolyLine<ScalarType, Dimension> OrientedPolyLineType;
  typedef NonOrientedPolyLine<ScalarType, Dimension> NonOrientedPolyLineType;
  typedef OrientedSurfaceMesh<ScalarType, Dimension> OrientedSurfaceMeshType;
  typedef NonOrientedSurfaceMesh<ScalarType, Dimension> NonOrientedSurfaceMeshType;
  typedef LinearInterpImage<ScalarType, Dimension> LinearInterpImageType;
  typedef ParametricImage<ScalarType, Dimension> ParametricImageType;
  typedef SSDImage<ScalarType, Dimension> SSDImageType;
  typedef LCCImage<ScalarType, Dimension> LCCImageType;
  typedef EQLAImage<ScalarType, Dimension> EQLAImageType;
  typedef MutualInformationImage<ScalarType, Dimension> MutualInformationImageType;
  typedef OrientedVolumeMesh<ScalarType, Dimension> OrientedVolumeMeshType;

  typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;

  std::string ObjectTypeStr = m_ParamObject->GetDeformableObjectType();
  const char *ObjectType = ObjectTypeStr.c_str();

  if (itksys::SystemTools::Strucmp(ObjectType, "Landmark") == 0) {
    std::shared_ptr<LandmarkType> objectLandmark = std::make_shared<LandmarkType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

      objectLandmark->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectLandmark->SetPolyData(PolyData);
    }
    m_Object = objectLandmark;
  } else if (itksys::SystemTools::Strucmp(ObjectType, "PointCloud") == 0) {
    std::shared_ptr<PointCloudType> objectPointCloud = std::make_shared<PointCloudType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

      objectPointCloud->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectPointCloud->SetPolyData(PolyData);

      objectPointCloud->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectPointCloud->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectPointCloud;
  } else if (itksys::SystemTools::Strucmp(ObjectType, "OrientedPolyLine") == 0) {
    std::shared_ptr<OrientedPolyLineType> objectOrientedPolyLine = std::make_shared<OrientedPolyLineType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

      objectOrientedPolyLine->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectOrientedPolyLine->SetPolyData(PolyData);

      objectOrientedPolyLine->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectOrientedPolyLine->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectOrientedPolyLine;
  } else if (itksys::SystemTools::Strucmp(ObjectType, "NonOrientedPolyLine") == 0) {
    std::shared_ptr<NonOrientedPolyLineType> objectNonOrientedPolyLine = std::make_shared<NonOrientedPolyLineType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

      objectNonOrientedPolyLine->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectNonOrientedPolyLine->SetPolyData(PolyData);

      objectNonOrientedPolyLine->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectNonOrientedPolyLine->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectNonOrientedPolyLine;
  } else if (itksys::SystemTools::Strucmp(ObjectType, "OrientedSurfaceMesh") == 0) {
    bool ReOrientSurface = m_ParamObject->ReOrient();

    std::shared_ptr<OrientedSurfaceMeshType> objectOrientedSurfaceMesh = std::make_shared<OrientedSurfaceMeshType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

      objectOrientedSurfaceMesh->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectOrientedSurfaceMesh->SetPolyData(PolyData);
      (ReOrientSurface) ? objectOrientedSurfaceMesh->SetReorient() : objectOrientedSurfaceMesh->UnSetReorient();

      objectOrientedSurfaceMesh->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectOrientedSurfaceMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectOrientedSurfaceMesh;
  } else if (itksys::SystemTools::Strucmp(ObjectType, "NonOrientedSurfaceMesh") == 0) {
    std::shared_ptr<NonOrientedSurfaceMeshType>
        objectNonOrientedSurfaceMesh = std::make_shared<NonOrientedSurfaceMeshType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

      objectNonOrientedSurfaceMesh->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectNonOrientedSurfaceMesh->SetPolyData(PolyData);

      objectNonOrientedSurfaceMesh->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectNonOrientedSurfaceMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectNonOrientedSurfaceMesh;
  }

    //Eventhough we see 2-current in R^2 and 3-current in R^ as the same object, we need to adapt the reader to vtkPolyData or vtkUnstructuredGrid
    //hence the condition on the dimension.
  else if (itksys::SystemTools::Strucmp(ObjectType, "OrientedVolumeMesh") == 0 && Dimension == 3) {
    std::shared_ptr<OrientedVolumeMeshType> objectOrientedVolumeMesh = std::make_shared<OrientedVolumeMeshType>();
    {
      vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
//      vtkSmartPointer<vtkUnstructuredGrid> unstructuredData = reader->GetOutput();
      auto * unstructuredData = reader->GetOutput();

      objectOrientedVolumeMesh->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectOrientedVolumeMesh->SetPolyData((vtkPointSet *) unstructuredData);

      objectOrientedVolumeMesh->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectOrientedVolumeMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectOrientedVolumeMesh;
  } else if (itksys::SystemTools::Strucmp(ObjectType, "OrientedVolumeMesh") == 0 && Dimension == 2) {
    std::shared_ptr<OrientedVolumeMeshType> objectOrientedVolumeMesh = std::make_shared<OrientedVolumeMeshType>();
    {
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(m_FileName);
      reader->Update();
      vtkSmartPointer<vtkPolyData> unstructuredData = reader->GetOutput();
      objectOrientedVolumeMesh->SetAnatomicalCoordinateSystem(m_ParamObject->GetAnatomicalCoordinateSystem());
      objectOrientedVolumeMesh->SetPolyData(unstructuredData);
      objectOrientedVolumeMesh->SetKernelType(this->StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
      objectOrientedVolumeMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
    }
    m_Object = objectOrientedVolumeMesh;
  } else if ((itksys::SystemTools::Strucmp(ObjectType, "SSDImage") == 0) ||
      (itksys::SystemTools::Strucmp(ObjectType, "LCCImage") == 0) ||
      (itksys::SystemTools::Strucmp(ObjectType, "EQLAImage") == 0) ||
      (itksys::SystemTools::Strucmp(ObjectType, "ParametricImage") == 0) ||
      (itksys::SystemTools::Strucmp(ObjectType, "MutualInformationImage") == 0) ) {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ImageTypePointer ExampleImg;

    std::shared_ptr<LinearInterpImageType> objectImg;
    if (itksys::SystemTools::Strucmp(ObjectType, "SSDImage") == 0 ||
        itksys::SystemTools::Strucmp(ObjectType, "ParametricImage") == 0) {
      std::shared_ptr<SSDImageType> objectImgAux = std::make_shared<SSDImageType>();
      objectImg = objectImgAux;
    } else if (itksys::SystemTools::Strucmp(ObjectType, "EQLAImage") == 0) {
      std::shared_ptr<EQLAImageType> objectImgAux = std::make_shared<EQLAImageType>();
      objectImgAux->SetEQLAKernelWidth(m_ParamObject->GetKernelWidth());
      objectImg = objectImgAux;
    } else if (itksys::SystemTools::Strucmp(ObjectType, "MutualInformationImage") == 0) {
      std::shared_ptr<MutualInformationImageType> objectImgAux = std::make_shared<MutualInformationImageType>();
      objectImg = objectImgAux;
    } else {
      std::shared_ptr<LCCImageType> objectImgAux = std::make_shared<LCCImageType>();
      objectImgAux->SetLCCKernelWidth(m_ParamObject->GetKernelWidth());
      objectImg = objectImgAux;
    }
/*
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(m_FileName);

		typedef itk::NiftiImageIO NiftiIOType;
		typename NiftiIOType::Pointer NiftiIOPointer = NiftiIOType::New();

		ImageTypePointer img;
		if ( NiftiIOPointer->CanReadFile(m_FileName) )
		{
			NiftiIOPointer->SetFileName(m_FileName);
			reader->SetImageIO(NiftiIOPointer);
			reader->Update();
			img = reader->GetOutput();
			// std::cout << "Read " << m_FileName << " as Nifti image" << std::endl;
		}
		else
		{
			reader->Update();
			img = reader->GetOutput();
		}

		typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
		typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
		minMaxFilter->SetInput(img);
		minMaxFilter->Update();
		double minI = minMaxFilter->GetMinimum();
		double maxI = minMaxFilter->GetMaximum();

		typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
		typename RescalerType::Pointer resf = RescalerType::New();
		resf->SetInput(img);
		resf->SetOutputMinimum(0.0);
		resf->SetOutputMaximum(1.0);
		resf->Update();
		ImageTypePointer rescaledImage = resf->GetOutput();

		ScalarType ImageGridDownsampling = m_ParamObject->GetImageGridDownsampling();

		objectImg->SetImageAndDownSamplingFactor(rescaledImage, ImageGridDownsampling);
		objectImg->SetMinIntensityOutput(minI);
		objectImg->SetMaxIntensityOutput(maxI);

		typename ImageType::PointType minAux = img->GetOrigin();
		typename ImageType::PointType maxAux = minAux;
		typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
		typename ImageType::SpacingType spacing = img->GetSpacing();

		std::cout << "Image read: origin = " << img->GetOrigin() << " size = " << size << " spacing = " << spacing << " min value = " << minI << " max value = " << maxI << ", rescaled between 0 and 1" << std::endl;
		std::cout << "orientation =\n" << img->GetDirection() << std::endl;
		// minMaxFilter->SetInput(objectImg->GetImage());
		// minI = minMaxFilter->GetMinimum();
		// maxI = minMaxFilter->GetMaximum();
		// std::cout << "after renormalization " << minI << " " << maxI << std::endl;

		itk::Matrix<double,Dimension,Dimension> dir = img->GetDirection();
		std::vector<int> permutation(Dimension);
		std::vector<int> flip_axes(Dimension);
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
		objectImg->SetPermutationAxes(permutation);
		objectImg->SetFlipAxes(flip_axes);

		// std::cout << "permutation\n" << permutation[0] << permutation[1] << permutation[2] << std::endl;
		// std::cout << "flip axes\n" << flip_axes[0] << flip_axes[1] << flip_axes[2] << std::endl;

		// MatrixType Yex = GridFunctionsType::ImageToPoints(img);
		// for (unsigned int d=0; d<Dimension; d++)
		// {
		// 	m_Max[d] = Yex.get_column(d).max_value();
		// 	m_Min[d] = Yex.get_column(d).min_value();
		// 	// wrong if axes are flipped
		// 	// m_Max[d] = maxAux[d] + size[d]*spacing[d];
		// 	// m_Min[d] = minAux[d];
		// }

		// if (m_IsTemplate)
		// {
		// ScalarType ImageGridDownsampling = m_ParamObject->GetImageGridDownsampling();
		// ImageTypePointer ExampleDownsampledImage;
		// if (ImageGridDownsampling <= 1.0)
		// 	ExampleDownsampledImage = img;
		// else
		// 	ExampleDownsampledImage = GridFunctionsType::DownsampleImage(img, ImageGridDownsampling);
		//
		// 	// MatrixType DownsampledY1 = GridFunctionsType::ImageToPoints(ExampleDownsampledImage);
		//
		// 	// std::cout << "DownsampledImage = origin " <<  ExampleDownsampledImage->GetOrigin() << " size = " << ExampleDownsampledImage->GetLargestPossibleRegion().GetSize() << " spacing = " << ExampleDownsampledImage->GetSpacing() << std::endl;
		// objectImg->SetDownSampledReferenceImage(ExampleDownsampledImage);
		// // }


		m_Object = objectImg;*/

    {

      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(m_FileName);

      typedef itk::NiftiImageIO NiftiIOType;
      typename NiftiIOType::Pointer NiftiIOPointer = NiftiIOType::New();

      ImageTypePointer img, inputimg;
      // img = working image, which can be the same as the original image, if no orientation issues.
      if (NiftiIOPointer->CanReadFile(m_FileName)) {
        //ImageTypePointer inputimg;
        NiftiIOPointer->SetFileName(m_FileName);
        reader->SetImageIO(NiftiIOPointer);
        reader->Update();
        inputimg = reader->GetOutput();

        // Reorient image if necessary, saving original coordinate system info in the final corresponding deformable object
        itk::Matrix<double, Dimension, Dimension> dir = inputimg->GetDirection();
        itk::Matrix<double, Dimension, Dimension> id;
        id.SetIdentity();

        // If not LPS orientation
        if (dir != id) {
          if (Dimension != 3) {
            std::cout << "Warning: cannot reorient images that are not 3D" << std::endl;
          } else {
            MatrixType orient_mtx(Dimension, Dimension, 0.0);
            for (unsigned int i = 0; i < Dimension; i++)
              for (unsigned int j = 0; j < Dimension; j++)
                orient_mtx(i, j) = dir(i, j);

            //std::cout << "Orientation matrix " << orient_mtx << std::endl;
            objectImg->SetAnatomicalCoordinateSystem(orient_mtx);
            std::cout << "Reorienting input image from " << objectImg->GetAnatomicalCoordinateSystemLabel()
                      << " orientation to the internal LPS standard." << std::endl;

            //std::cout << "Deformable Object is an image originally oriented in the " << objectImg->GetAnatomicalCoordinateSystemLabel() << " coordinate system." << std::endl;
            //std::cout << "Reorienting image to match the LPS standard." << std::endl;

            // reorient image (LPS orientation / RAI code for ITK filter)
            typedef itk::Image<ScalarType, 3> ImageType_dim_3;
            typename itk::OrientImageFilter<ImageType_dim_3, ImageType_dim_3>::Pointer
                orienter = itk::OrientImageFilter<ImageType_dim_3, ImageType_dim_3>::New();
            orienter->UseImageDirectionOn();
            orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
            orienter->SetInput(dynamic_cast<ImageType_dim_3 *>( inputimg.GetPointer()));
            orienter->Update();
            img = dynamic_cast<ImageType *>( orienter->GetOutput());

            /* typedef itk::ImageFileWriter<ImageType> WriterType;
                        typename WriterType::Pointer writer = WriterType::New();
                        writer->SetInput(img); //
                        writer->SetFileName(string(m_FileName).append("LPS.nii"));
                        writer->Update(); */

            //std::cout << "Reoriented image matrix (LPS): " << img << std::endl; // << " " << (*outimg)->GetDirection() << std::endl;
          }
        }
          //std::cout << "Read " << m_FileName << " as Nifti image" << std::endl;
        else
          img = inputimg;
      } else {
        reader->Update();
        img = reader->GetOutput();
      }

      //std::cout << "Working image = " << img << std::endl;

      typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
      typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
      minMaxFilter->SetInput(img);
      minMaxFilter->Update();
      double minI = minMaxFilter->GetMinimum();
      double maxI = minMaxFilter->GetMaximum();

      typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
      typename RescalerType::Pointer resf = RescalerType::New();
      resf->SetInput(img);
      resf->SetOutputMinimum(0.0);
      resf->SetOutputMaximum(1.0);
      resf->Update();
      ImageTypePointer rescaledImage = resf->GetOutput();

      ScalarType ImageGridDownsampling = m_ParamObject->GetImageGridDownsampling();

      objectImg->SetImageAndDownSamplingFactor(rescaledImage, ImageGridDownsampling);
      objectImg->SetMinIntensityOutput(minI);
      objectImg->SetMaxIntensityOutput(maxI);

//			typename ImageType::PointType minAux = img->GetOrigin();
//			typename ImageType::PointType maxAux = minAux;
      typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
      typename ImageType::SpacingType spacing = img->GetSpacing();

//      std::cout << "Image read: origin = " << img->GetOrigin() << " size = " << size << " spacing = " << spacing
//                << " min value = " << minI << " max value = " << maxI << ", rescaled between 0 and 1" << std::endl;

      // At this stage, img->GetDirection should be the identity matrix for a 3D Nifti image. The following could be useful for other types of images.
      itk::Matrix<double, Dimension, Dimension> dir = img->GetDirection();
      std::vector<unsigned int> permutation(Dimension);
      std::vector<int> flip_axes(Dimension);
      for (unsigned int d = 0; d < Dimension; d++) {
        unsigned int ind = 0;
        double val = fabs(dir(0, d));
        for (unsigned int k = 1; k < Dimension; k++) {
          if (fabs(dir(k, d)) > val) {
            ind = k;
            val = fabs(dir(k, d));
          }
        }
        permutation[d] = ind;
        flip_axes[d] = (dir(ind, d) > 0.0) ? 1 : -1;
      }
      objectImg->SetPermutationAxes(permutation);
      objectImg->SetFlipAxes(flip_axes);

      // std::cout << "permutation\n" << permutation[0] << permutation[1] << permutation[2] << std::endl;
      // std::cout << "flip axes\n" << flip_axes[0] << flip_axes[1] << flip_axes[2] << std::endl;

/*			MatrixType Yex = GridFunctionsType::ImageToPoints(img);
			for (unsigned int d=0; d<Dimension; d++)
			{
				m_Max[d] = Yex.get_column(d).max_value();
				m_Min[d] = Yex.get_column(d).min_value();
				// wrong if axes are flipped
				// m_Max[d] = maxAux[d] + size[d]*spacing[d];
				// m_Min[d] = minAux[d];
			}
			*/

    }
    m_Object = objectImg;
  } else {
    std::cerr << "Unknown object type: " << m_ParamObject->GetDeformableObjectType() <<
              "\nCurrently available types are: SSDImage, LCCImage, EQLAImage, ParametricImage, MutualInformationImage, Landmark, PointCloud, OrientedPolyLine, "
                  "NonOrientedPolyLine, OrientedSurfaceMesh, NonOrientedSurfaceMesh, OrientedVolumeMesh." << std::endl;
    return;
  }

  if (m_IsTemplate && itksys::SystemTools::Strucmp(ObjectType, "ParametricImage") == 0) {
    std::shared_ptr<LinearInterpImageType> rawTemplateImage = std::static_pointer_cast<LinearInterpImageType>(m_Object);
    std::shared_ptr<ParametricImageType> parametricTemplateImage = std::make_shared<ParametricImageType>();

    parametricTemplateImage
        ->SetPhotometricKernelType(StringToKernelEnumType(m_ParamObject->GetKernelType().c_str()));
    parametricTemplateImage->SetPhotometricKernelWidth(m_ParamObject->GetKernelWidth());
    parametricTemplateImage->SetPhotometricCPSpacing(m_ParamObject->GetPhotometricCPSpacing());
    parametricTemplateImage->Initialize(rawTemplateImage);

    m_Object = parametricTemplateImage;
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
KernelEnumType
DeformableObjectReader<ScalarType, Dimension>
::StringToKernelEnumType(const char *kernelType) {
  KernelEnumType result = null;
  if (itksys::SystemTools::Strucmp(kernelType, "p3m") == 0) { result = P3M; }
#ifdef USE_CUDA
  else if (itksys::SystemTools::Strucmp(kernelType, "cudaexact") == 0) { result = CUDAExact; }
#endif
  else {
    if (itksys::SystemTools::Strucmp(kernelType, "exact") != 0)
      std::cerr << "Unknown kernel type for the deformable object : defaulting to exact" << std::endl;
    result = Exact;
  }
  return result;
}

template class DeformableObjectReader<double,2>;
template class DeformableObjectReader<double,3>;
//#endif
