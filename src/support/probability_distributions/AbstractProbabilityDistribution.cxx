/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "AbstractProbabilityDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType>
AbstractProbabilityDistribution<ScalarType>
::AbstractProbabilityDistribution()
{}


template <class ScalarType>
AbstractProbabilityDistribution<ScalarType>
::~AbstractProbabilityDistribution()
{}


template <class ScalarType>
AbstractProbabilityDistribution<ScalarType>
::AbstractProbabilityDistribution(const AbstractProbabilityDistribution& other)
{}

template class AbstractProbabilityDistribution<double>;
