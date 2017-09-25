/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "AbstractSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
AbstractSampler<ScalarType, Dimension>
::AbstractSampler() {}

template<class ScalarType, unsigned int Dimension>
AbstractSampler<ScalarType, Dimension>
::~AbstractSampler() {}

template<class ScalarType, unsigned int Dimension>
AbstractSampler<ScalarType, Dimension>
::AbstractSampler(const AbstractSampler &other) {}

template
class AbstractSampler<double, 2>;
template
class AbstractSampler<double, 3>;

