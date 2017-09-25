/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "AbstractGeometry.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
AbstractGeometry<ScalarType, Dimension>
::AbstractGeometry() {
  m_Type = null;
  m_Modified = true;
  m_BoundingBox.set_size(Dimension, 2);
  m_AnatomicalOrientation.SetAnatomicalCoordinateSystemLabel("LPS");

}

template<class ScalarType, unsigned int Dimension>
AbstractGeometry<ScalarType, Dimension>
::AbstractGeometry(const AbstractGeometry &other) {
  m_Type = other.m_Type;
  m_Modified = other.m_Modified;
  m_BoundingBox = other.m_BoundingBox;
  m_AnatomicalOrientation.SetAnatomicalCoordinateSystemLabel(other.m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel());

// m_Timepoint = other.m_Timepoint;
// m_TimeIndex = other.m_TimeIndex;
}

template<class ScalarType, unsigned int Dimension>
AbstractGeometry<ScalarType, Dimension>
::AbstractGeometry(const AbstractGeometry &example, const MatrixType &LandmarkOrImagePoints) {
  m_Type = example.m_Type;
  m_Modified = true;
  m_BoundingBox = example.m_BoundingBox;
  m_AnatomicalOrientation = example.m_AnatomicalOrientation;
}

template<class ScalarType, unsigned int Dimension>
AbstractGeometry<ScalarType, Dimension>
::~AbstractGeometry() {

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
MatrixType
AbstractGeometry<ScalarType, Dimension>
::UnionBoundingBox(MatrixType BB1, MatrixType BB2) {
  MatrixType Union(Dimension, 2);
  for (int d = 0; d < Dimension; d++) {
    Union(d, 0) = (BB1(d, 0) < BB2(d, 0) ? BB1(d, 0) : BB2(d, 0));
    Union(d, 1) = (BB1(d, 1) > BB2(d, 1) ? BB1(d, 1) : BB2(d, 1));
  }
  return Union;
}

template class AbstractGeometry<double,2>;
template class AbstractGeometry<double,3>;

