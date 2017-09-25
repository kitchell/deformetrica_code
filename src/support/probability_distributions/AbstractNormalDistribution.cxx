/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "AbstractNormalDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
AbstractNormalDistribution<ScalarType>
::AbstractNormalDistribution() {
  m_Mean.set_size(1);
  m_Mean(0) = sqrt(-1); // Nan.
}

template<class ScalarType>
AbstractNormalDistribution<ScalarType>
::~AbstractNormalDistribution() {}

template<class ScalarType>
AbstractNormalDistribution<ScalarType>
::AbstractNormalDistribution(const AbstractNormalDistribution &other) {
  m_Mean = other.m_Mean;
}

template class AbstractNormalDistribution<double>;
