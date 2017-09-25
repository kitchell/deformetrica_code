/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformableMultiObject.h"

#include "AbstractGeometry.h"
#include "Landmark.h"
#include "LinearInterpImage.h"

#include "GridFunctions.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#include "KernelFactory.h"

#include "LinearAlgebra.h"

using namespace def::algebra;


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////



template <class ScalarType, unsigned int Dimension>
DeformableMultiObject<ScalarType, Dimension>
::DeformableMultiObject()
 {
	m_ImageIndex = -1;

	m_NumberOfImageKindObjects = 0;
	m_NumberOfLandmarkKindObjects = 0;
	m_NumberOfObjects = 0;
	
	m_NumberOfLandmarkPoints = 0;
	m_NumberOfImagePoints = 0;
	
	m_ObjectList.resize(0);
	m_NumberOfPoints.resize(0);
	
 }



template <class ScalarType, unsigned int Dimension>
DeformableMultiObject<ScalarType, Dimension>
::~DeformableMultiObject()
{}



template <class ScalarType, unsigned int Dimension>
DeformableMultiObject<ScalarType, Dimension>
::DeformableMultiObject(const DeformableMultiObject& other)
{
	// Deep copy so we can safely modify the objects
	m_ObjectList.resize(other.m_ObjectList.size());
	for (unsigned int i = 0; i < other.m_ObjectList.size(); i++)
		m_ObjectList[i] = other.m_ObjectList[i]->Clone();

	m_NumberOfObjects = other.m_NumberOfObjects;
	m_NumberOfImageKindObjects = other.m_NumberOfImageKindObjects;
	m_NumberOfLandmarkKindObjects = other.m_NumberOfLandmarkKindObjects;

	m_NumberOfLandmarkPoints = other.m_NumberOfLandmarkPoints;
	m_NumberOfImagePoints = other.m_NumberOfImagePoints;

	m_NumberOfPoints.resize(m_ObjectList.size());
	for (unsigned int i = 0; i < m_ObjectList.size(); i++)
		m_NumberOfPoints[i] = other.m_NumberOfPoints[i];

	m_ImageIndex = other.m_ImageIndex;

	m_BoundingBox = other.m_BoundingBox;
}


template <class ScalarType, unsigned int Dimension>
DeformableMultiObject<ScalarType, Dimension>
::DeformableMultiObject(const DeformableMultiObject& Example, const MatrixType& LandmarkPoints, const MatrixType& DownSampledImageMap)
{
	if (Example.m_NumberOfLandmarkPoints != LandmarkPoints.rows() || Example.m_NumberOfImagePoints != DownSampledImageMap.rows() )
		throw std::runtime_error("Number of Landmark and/or Image points mismatch in copy/update DeformableMultiObject constructor ");

	MatrixListType List(Example.m_NumberOfObjects);
	int counter = 0;
	for (int i = 0; i < Example.m_NumberOfObjects; i++)
	{
		if (Example.m_ObjectList[i]->IsOfLandmarkKind())
		{
			List[i] = LandmarkPoints.get_n_rows(counter, Example.m_NumberOfPoints[i]);
			counter += Example.m_NumberOfPoints[i];
		}
		else if (Example.m_ObjectList[i]->IsOfImageKind())
		{
			List[i] = DownSampledImageMap;
		}
		else
			throw std::runtime_error("unknown object type");
	}
			
	m_ObjectList.resize(Example.m_ObjectList.size());
	for (unsigned int i = 0; i < Example.m_ObjectList.size(); i++)
		m_ObjectList[i] = Example.m_ObjectList[i]->DeformedObject(List[i]);
		
	this->Update();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::SetObjectList(AbstractGeometryList& objList)
{
	// Deep copy so we can safely modify the objects
	m_ObjectList.resize(objList.size());
	for (unsigned int i = 0; i < objList.size(); i++)
		m_ObjectList[i] = objList[i]->Clone();
}


template <class ScalarType, unsigned int Dimension>
std::vector<bool>
DeformableMultiObject<ScalarType, Dimension>
::IsOfImageKind() const
{
	std::vector<bool> result(m_NumberOfObjects);
	for (unsigned int i = 0; i < m_NumberOfObjects; i++)
		result[i] = m_ObjectList[i]->IsOfImageKind();
	return result;
}


template <class ScalarType, unsigned int Dimension>
std::vector<bool>
DeformableMultiObject<ScalarType, Dimension>
::IsOfLandmarkKind() const
{
	std::vector<bool> result(m_NumberOfObjects);
	for (unsigned int i = 0; i < m_NumberOfObjects; i++)
		result[i] = m_ObjectList[i]->IsOfLandmarkKind();
	return result;
}


template <class ScalarType, unsigned int Dimension>
MatrixType
DeformableMultiObject<ScalarType, Dimension>
::GetLandmarkPoints() const
{
	if(m_NumberOfLandmarkKindObjects == 0)
		return MatrixType(0,0);


	MatrixType X0(m_NumberOfLandmarkPoints,Dimension);
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			MatrixType Xaux = std::static_pointer_cast<LandmarkType>(m_ObjectList[i])->GetPointCoordinates();
			for (int r = 0; r < m_NumberOfPoints[i] ; r++)
				X0.set_row(counter + r, Xaux.get_row(r));
			counter += m_NumberOfPoints[i];
		}
	return X0;
}

template <class ScalarType, unsigned int Dimension>
MatrixType
DeformableMultiObject<ScalarType, Dimension>
::GetDownSampledImageMap() const
{
	if (m_NumberOfImageKindObjects == 0)
		return MatrixType(0,0);

	if (m_NumberOfImageKindObjects > 1)
		throw std::runtime_error("In DeformableMultiObject::GetDownSampledImageMap() - There are more than one image (it should not be).");

	if (m_ObjectList[m_ImageIndex]->IsLinearInterpImage())
	{
		std::shared_ptr<LIImageType> tempLII = std::static_pointer_cast<LIImageType>(m_ObjectList[m_ImageIndex]);
		return tempLII->GetDownSampledImageMap();
	}
	else // Parametric image.
	{
		std::shared_ptr<ParametricImageType> tempPI = std::static_pointer_cast<ParametricImageType>(m_ObjectList[m_ImageIndex]);
		return tempPI->GetDownSampledImageMap();
	}
}

template <class ScalarType, unsigned int Dimension>
typename DeformableMultiObject<ScalarType, Dimension>::ImageTypePointer
DeformableMultiObject<ScalarType, Dimension>
::GetDownSampledImage() const
{
	if(m_NumberOfImageKindObjects == 0)
		return NULL;

	if(m_NumberOfImageKindObjects > 1)
		throw std::runtime_error("In DeformableMultiObject::GetDownSampledImage() - There are more than one image (it should not be).");

	if (m_ObjectList[m_ImageIndex]->IsLinearInterpImage())
	{
		std::shared_ptr<LIImageType> tempLII = std::static_pointer_cast<LIImageType>(m_ObjectList[m_ImageIndex]);
		return tempLII->GetDownSampledImage();
	}
	else // Parametric image.
	{
		std::shared_ptr<ParametricImageType> tempPI = std::static_pointer_cast<ParametricImageType>(m_ObjectList[m_ImageIndex]);
		return tempPI->GetDownSampledImage();
	}
}

template<class ScalarType, unsigned int Dimension>
typename DeformableMultiObject<ScalarType, Dimension>::ImageTypePointer
DeformableMultiObject<ScalarType, Dimension>
::GetImage() const
{
	if (m_NumberOfImageKindObjects == 0)
		return NULL;

	if(m_NumberOfImageKindObjects > 1)
		throw std::runtime_error("In DeformableMultiObject::GetImage() - There are more than one image (it should not be).");

	if (m_ObjectList[m_ImageIndex]->IsLinearInterpImage())
	{
		std::shared_ptr<LIImageType> tempLII = std::static_pointer_cast<LIImageType>(m_ObjectList[m_ImageIndex]);
		return tempLII->GetImage();
	}
	else // Parametric image.
	{
		std::shared_ptr<ParametricImageType> tempPI = std::static_pointer_cast<ParametricImageType>(m_ObjectList[m_ImageIndex]);
		return tempPI->GetImage();
	}
}


template <class ScalarType, unsigned int Dimension>
typename std::vector<unsigned long>
DeformableMultiObject<ScalarType, Dimension>
::GetDimensionOfDiscretizedObjects() const
 {
	 std::vector<unsigned long> out(m_NumberOfObjects);
	 for (unsigned int i = 0; i < m_NumberOfObjects; i++)
		 out[i] = m_ObjectList[i]->GetDimensionOfDiscretizedObject();
	 
	 return out;
 }

template <class ScalarType, unsigned int Dimension>
MatrixListType
DeformableMultiObject<ScalarType, Dimension>
::GetImageIntensityAndLandmarkPointCoordinates() const
{
	MatrixListType Y(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			std::shared_ptr<LandmarkType> obj = std::static_pointer_cast<LandmarkType>(m_ObjectList[i]);
			Y[i] = obj->GetPointCoordinates();
		}
		else if ( m_ObjectList[i]->IsLinearInterpImage() )
		{
			std::shared_ptr<LIImageType> obj = std::static_pointer_cast<LIImageType>(m_ObjectList[i]);
			MatrixType Yaux(obj->GetImage()->GetLargestPossibleRegion().GetNumberOfPixels(), 1);
			Yaux.set_column(0, GridFunctionsType::VectorizeImage(obj->GetImage()) );
			Y[i] = Yaux;
		}
		else if ( m_ObjectList[i]->IsParametricImage() )
		{
			std::shared_ptr<ParametricImageType> obj = std::static_pointer_cast<ParametricImageType>(m_ObjectList[i]);
			Y[i] = obj->GetPhotometricWeights();
		}
		else
			throw std::runtime_error("Unknown object type.");
	}

	return Y;
}

template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::UpdateImageIntensityAndLandmarkPointCoordinates(const MatrixListType& Y)
 {
	if (Y.size() != m_NumberOfObjects)
		throw std::runtime_error("size of list and number of objects mismatch");

	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			std::shared_ptr<LandmarkType> obj = std::static_pointer_cast<LandmarkType>(m_ObjectList[i]);
			obj->UpdatePolyDataWithPointCoordinates(Y[i]);
		}
		else if ( m_ObjectList[i]->IsLinearInterpImage() )
		{
			std::shared_ptr<LIImageType> obj = std::static_pointer_cast<LIImageType>(m_ObjectList[i]);
			obj->UpdateImageIntensity(Y[i]);
		}
		else if ( m_ObjectList[i]->IsParametricImage() )
		{
			std::shared_ptr<ParametricImageType> obj = std::static_pointer_cast<ParametricImageType>(m_ObjectList[i]);
			obj->UpdateImageIntensity(Y[i].get_column(0));
		}
		else
			throw std::runtime_error("Unknown object type");
	}
 }

template <class ScalarType, unsigned int Dimension>
MatrixType
DeformableMultiObject<ScalarType, Dimension>
::SplatDifferenceImage(const std::shared_ptr<DeformableMultiObject> targ, const MatrixType& downSampledImageMap) const
{

	if (m_NumberOfImageKindObjects)
	{
        std::shared_ptr<LIImageType> targImage
            = std::static_pointer_cast<LIImageType>(targ->GetObjectList()[targ->GetImageIndex()]);
        MatrixType splat;

        if (m_ObjectList[m_ImageIndex]->IsParametricImage())
        {
            std::shared_ptr<ParametricImageType> tempImage
                = std::static_pointer_cast<ParametricImageType>(m_ObjectList[m_ImageIndex]);
            splat = tempImage->SplatDifferenceImage(targImage, downSampledImageMap);
        }
        else // Linear interpolation template.
        {
            std::shared_ptr<LIImageType> tempImage = std::static_pointer_cast<LIImageType>(m_ObjectList[m_ImageIndex]);
            splat = tempImage->SplatDifferenceImage(targImage, downSampledImageMap);
        }
        return splat;
	}
	return MatrixType(0,0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::Update()
{
	m_NumberOfObjects = m_ObjectList.size();

	m_NumberOfImageKindObjects = 0;
	m_NumberOfLandmarkKindObjects = 0;
	m_NumberOfLandmarkPoints = 0;
	m_NumberOfImagePoints = 0;

	m_NumberOfPoints.resize(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		// Update each object in the list
		if (m_ObjectList[i] == NULL)
			throw std::runtime_error("All deformable objects have not been set");
		
		m_ObjectList[i]->Update();
		
		int nb = m_ObjectList[i]->GetNumberOfPoints();
		m_NumberOfPoints[i] = nb;

		if ( m_ObjectList[i]->IsOfImageKind() )
		{
			m_NumberOfImageKindObjects += 1;
			m_NumberOfImagePoints = nb;
			m_ImageIndex = i;
		}

		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{	
			m_NumberOfLandmarkPoints += nb;
			m_NumberOfLandmarkKindObjects += 1;
		}
	}

	if ( m_NumberOfImageKindObjects > 1 )
		throw std::runtime_error("In DeformableMultiObject::SetObjectList() - There are more than one image (it should not be).");
	// in case of multi-modalities, it would be better to add a new class multi-image, since the same image domain should be shared by all images
	// and the eta maps in TransportAlongGeodesic could be added, as well as the initial eta0 maps (provided they are weighted by their respective dataSigma)
	
	this->UpdateBoundingBox();
	
}

template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::UpdateBoundingBox()
{
	MatrixType BB(Dimension,2);
	
	BB = m_ObjectList[0]->GetBoundingBox();
	
	for (int i = 1; i < m_NumberOfObjects; i++)
	{
		for (int d = 0; d < Dimension; d++)
		{
			MatrixType BBaux = m_ObjectList[i]->GetBoundingBox();
			BB(d,0) = (BB(d,0)<BBaux(d,0)?BB(d,0):BBaux(d,0));
			BB(d,1) = (BB(d,1)>BBaux(d,1)?BB(d,1):BBaux(d,1));
		}
	}

	m_BoundingBox = BB;
}



template <class ScalarType, unsigned int Dimension>
typename std::vector<ScalarType>
DeformableMultiObject<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<DeformableMultiObject> target)
{
	if (m_NumberOfObjects != target->GetNumberOfObjects())
		throw std::runtime_error("number of objects mismatched");

	 AbstractGeometryList targetList = target->GetObjectList();

	std::vector<ScalarType> match(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		match[i] = m_ObjectList[i]->ComputeMatch(targetList[i]);
		
	return match;
}



template <class ScalarType, unsigned int Dimension>
MatrixListType
DeformableMultiObject<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<DeformableMultiObject> target)
{
	if (m_NumberOfObjects != target->GetNumberOfObjects())
		throw std::runtime_error("number of objects mismatched");

 	AbstractGeometryList targetList = target->GetObjectList();
	MatrixListType gradmatch(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		gradmatch[i] = m_ObjectList[i]->ComputeMatchGradient(targetList[i]);

	return gradmatch;
}



template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::ListToMatrices(const MatrixListType& L, MatrixType& MLandmark, MatrixType& MImage) const
{
	if (L.size() != m_NumberOfObjects)
		throw std::runtime_error("Number of objects mismatch in DeformableMultiObject::ListToMatrices.");
	
	int dim = L[0].columns();
	
	MLandmark.set_size(m_NumberOfLandmarkPoints, dim);
	MImage.set_size(m_NumberOfImagePoints, dim);
	
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if (m_ObjectList[i]->IsOfLandmarkKind())
		{
			for (int r = 0; r < m_NumberOfPoints[i]; r++)
			{
				MLandmark.set_row(counter + r, L[i].get_row(r));
			}
			counter += m_NumberOfPoints[i];
		}
		else if (m_ObjectList[i]->IsOfImageKind())
		{
			MImage = L[i];
		}
		else
			throw std::runtime_error("Unknown object type.");
	}
}


template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::MatricesToList(MatrixListType& L, const MatrixType& MLandmark, const MatrixType& MImage) const
{
	L.resize(m_NumberOfObjects);
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if (m_ObjectList[i]->IsOfLandmarkKind())
		{
			L[i] = MLandmark.get_n_rows(counter, m_NumberOfPoints[i]);
			counter += m_NumberOfPoints[i];
		}
		else if (m_ObjectList[i]->IsOfImageKind())
		{
			L[i] = MImage;
		}
		else
			throw std::runtime_error("Unknown object type.");
	}
}



template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::WriteMultiObject(std::vector<std::string>& outfn) const
{
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->WriteObject(outfn[i]);
}


template <class ScalarType, unsigned int Dimension>
void
DeformableMultiObject<ScalarType, Dimension>
::WriteMultiObject(std::vector<std::string>& outfn, const MatrixListType& velocity) const
{
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->WriteObject(outfn[i], velocity);
}

template class DeformableMultiObject<double,2>;
template class DeformableMultiObject<double,3>;
