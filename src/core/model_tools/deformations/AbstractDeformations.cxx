/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/


#include "AbstractDeformations.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
AbstractDeformations<ScalarType, Dimension>
::AbstractDeformations() :
    m_Type(null), m_DeformableMultiObject(NULL), m_Modified(true), m_DeformableObjectModified(true) {}

template<class ScalarType, unsigned int Dimension>
AbstractDeformations<ScalarType, Dimension>
::~AbstractDeformations() {}

template<class ScalarType, unsigned int Dimension>
AbstractDeformations<ScalarType, Dimension>
::AbstractDeformations(const AbstractDeformations &other) {
  m_Type = other.m_Type;

  if (other.m_DeformableMultiObject != NULL) {
    m_DeformableMultiObject = other.m_DeformableMultiObject->Clone();
  } else {
    m_DeformableMultiObject = NULL;
  }

  m_LandmarkPoints = other.m_LandmarkPoints;
  m_ImagePoints = other.m_ImagePoints;


  m_IsLandmarkPoints = other.m_IsLandmarkPoints;
  m_IsImagePoints = other.m_IsImagePoints;

  m_Modified = other.m_Modified;
  m_DeformableObjectModified = other.m_DeformableObjectModified;
}



template<class ScalarType, unsigned int Dimension>
void
AbstractDeformations<ScalarType, Dimension>
::SetDeformableMultiObject(std::shared_ptr<DeformableMultiObjectType> DMO) {
  m_DeformableMultiObject = DMO;
  m_LandmarkPoints = DMO->GetLandmarkPoints();
  m_ImagePoints = DMO->GetDownSampledImageMap();
  m_DownSampledImage = DMO->GetDownSampledImage();
  m_Image = DMO->GetImage();

  m_IsLandmarkPoints = (DMO->GetNumberOfLandmarkKindObjects() != 0);
  m_IsImagePoints = (DMO->GetNumberOfImageKindObjects() == 1);

  m_DeformableObjectModified = true;
}

template class AbstractDeformations<double, 2>;
template class AbstractDeformations<double, 3>;
