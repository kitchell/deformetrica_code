/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

/// Core files.
#include "AbstractGeometry.h"

/// Input-output files.
#include "DeformableObjectParametersXMLFile.h"

/// Support files.
#include "KernelType.h"
#include "LinearAlgebra.h"

/// Non-core files.
#include "itkImage.h"
#include <cstring>
#include <iostream>
#include <sstream>


#include "LinearInterpImage.h"
#include "ParametricImage.h"
#include "SSDImage.h"
#include "LCCImage.h"
#include "EQLAImage.h"
#include "MutualInformationImage.h"
#include "Landmark.h"
#include "PointCloud.h"
#include "OrientedPolyLine.h"
#include "NonOrientedPolyLine.h"
#include "OrientedSurfaceMesh.h"
#include "NonOrientedSurfaceMesh.h"
#include "OrientedVolumeMesh.h"

#include "KernelFactory.h"

#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNiftiImageIO.h"
#include "itksys/SystemTools.hxx"
#include "itkOrientImageFilter.h"

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGridReader.h"

using namespace def::algebra;


template<class ScalarType, unsigned int Dimension>
class DeformableObjectReader {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// ITK Image type
  typedef itk::Image<ScalarType, Dimension> ImageType;
  /// ITK Image Pointer type
  typedef typename ImageType::Pointer ImageTypePointer;

  /// Deformable object type.
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  DeformableObjectReader();

  ~DeformableObjectReader();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the file name to \e fn.
  void SetFileName(char *fn) { m_FileName = fn; }
  void SetFileName(const std::string &fn) { m_FileName = (char *) fn.c_str(); }

  void SetObjectParameters(DeformableObjectParameters::Pointer param) { m_ParamObject = param; }

  void SetTemplateType() { m_IsTemplate = true; }

  /// Returns the ouput i.e. the deformable object.
  std::shared_ptr<AbstractGeometryType> GetOutput() { return m_Object; }

  // /// Returns the bounding box.
  // MatrixType GetBoundingBox() { return m_Bbox; }

  // /// Returns data-sigma squared.
  // double GetDataSigmaSquared() const { return (m_ParamObject->GetDataSigma() * m_ParamObject->GetDataSigma()); }
  // /************************************** CHANGE BEGIN ******************************************/
  // double GetWeightData() const { return m_ParamObject->GetWeightData(); }
  // double GetPriorData() const { return m_ParamObject->GetPriorData(); }
  // double GetKernelWidth() const { return m_ParamObject->GetKernelWidth(); }
  // /************************************** CHANGE END ******************************************/


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  void Update();

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  KernelEnumType StringToKernelEnumType(const char *kernelType);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Name of the file where the deformable object will be read.
  char *m_FileName;

  DeformableObjectParameters::Pointer m_ParamObject;

  /// Our deformable object which will be extracted from the file.
  std::shared_ptr<AbstractGeometryType> m_Object;


  // MatrixType m_Bbox;

  bool m_IsTemplate;

  // VectorType m_Min;
  // VectorType m_Max;


};

