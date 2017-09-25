/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "DirichletDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType>
DirichletDistribution<ScalarType>
::DirichletDistribution()
{
    m_ConcentrationParameters.set_size(1);
    m_ConcentrationParameters(0) = 1.0;
}


template <class ScalarType>
DirichletDistribution<ScalarType>
::~DirichletDistribution()
{}


template <class ScalarType>
DirichletDistribution<ScalarType>
::DirichletDistribution(const DirichletDistribution& other)
{
    m_ConcentrationParameters = other.m_ConcentrationParameters;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


template <class ScalarType>
ScalarType
DirichletDistribution<ScalarType>
::ComputeLogLikelihood(LinearVariableType const& obs) const
{
    return dot_product(m_ConcentrationParameters - 1, obs.vectorize());
}

template <class ScalarType>
ScalarType
DirichletDistribution<ScalarType>
::ComputeLogLikelihood(VectorType const& obs) const
{
    return dot_product(m_ConcentrationParameters - 1, obs);
}

template class DirichletDistribution<double>;
