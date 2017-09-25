/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformableObjectParameters.h"

#include "itksys/SystemTools.hxx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

DeformableObjectParameters
::DeformableObjectParameters() {
  m_DeformableObjectType = "unknown";

  m_DataSigma = 0.0;
  m_KernelType = "null";
  // m_P3MWorkingSpacingRatio = 0.2; // 1/5 gives a relative approximation error of about 5%, use 0.3 for increased speed
  // m_P3MPaddingFactor = 3.0f; // enlarge grids by 3 x kernelwidth to avoid side effects (FFTs have circular boundary conditions). It is also used to define a bounding box

  m_KernelWidth = 0.0;

  m_ImageGridDownsampling = 1.0;

  m_reOrient = true;

  m_DataSigma_Normalized_Hyperparameter = 0.2;
  m_DataSigma_Prior = 1.0;

  m_AnatomicalCoordinateSystem = "LPS";


//	m_UseParametricTemplateImage = false;
//	m_PhotometricKernelWidth = 0.0;
//	m_PhotometricKernelType = "exact";
  m_PhotometricCPSpacing = 0.0;
}

DeformableObjectParameters
::~DeformableObjectParameters() {

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
DeformableObjectParameters
::Update() {
  if (m_KernelWidth < 1e-20) {
    std::cout << "Warning : no value set for the deformable object kernel width. Defaulted to 1 (arbitrary). "
        "Ignore this warning if no kernel width is needed (e.g. SSDImage, Landmark objects)."
              << std::endl;
    m_KernelWidth = 1;
  }
  if (m_PhotometricCPSpacing < 1e-20) { m_PhotometricCPSpacing = m_KernelWidth; }
}

bool
DeformableObjectParameters
::CheckValues() {
  if (m_DataSigma < 1e-20)
    return false;

  if ((itksys::SystemTools::Strucmp(m_DeformableObjectType.c_str(), "Landmarks") != 0)
      && (itksys::SystemTools::Strucmp(m_DeformableObjectType.c_str(), "SSDImage") != 0)) {
    if (m_KernelWidth < 1e-20)
      return false;
  }

  if (m_ImageGridDownsampling < 1.0)
    return false;

  if (!(m_AnatomicalCoordinateSystem.size() == 0 || m_AnatomicalCoordinateSystem.size() == 3))
    return false;

  if (m_DataSigma_Normalized_Hyperparameter < 0.0)
    return false;

  if (m_DataSigma_Prior < 0.0)
    return false;

  if (!(m_AnatomicalCoordinateSystem.size() == 0 || m_AnatomicalCoordinateSystem.size() == 3))
    return false;

  return true;
}

void
DeformableObjectParameters
::PrintSelf(std::ostream &os) {
  os << "Deformable object type: " << m_DeformableObjectType << std::endl;
  os << "Data sigma = " << m_DataSigma << std::endl;
  os << "Data sigma prior = " << m_DataSigma_Prior << std::endl;
  os << "Data sigma normalized hyperparameter = " << m_DataSigma_Normalized_Hyperparameter << std::endl;
//	os << "Template smoothness prior = " << m_SmoothnessPrior << std::endl;
  os << std::endl;
  os << "Kernel type = " << m_KernelType << std::endl;
  os << "Kernel width = " << m_KernelWidth << std::endl;
  os << std::endl;
//	os << "P3M working spacing ratio = " << m_P3MWorkingSpacingRatio << std::endl;
//	os << "P3M padding factor = " << m_P3MPaddingFactor << std::endl;
  os << "Image grid downsampling = " << m_ImageGridDownsampling << std::endl;
  os << "Reorient normals: " << (m_reOrient ? "On" : "Off") << std::endl;
  os << "Anatomical Coordinate System: " << m_AnatomicalCoordinateSystem << std::endl;

//	os << "Use a parametric template (for images) = " << (m_UseParametricTemplateImage?"On":"Off") << std::endl;
//	os << "Photometric kernel width (for parametric template image) = " << m_PhotometricKernelWidth << std::endl;
//	os << "Photometric kernel type (for parametric template image) = " << m_PhotometricKernelType << std::endl;
  os << "Photometric CP spacing (for parametric template image) = " << m_PhotometricCPSpacing << std::endl;
  os << std::endl;
}
