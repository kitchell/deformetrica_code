/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "InverseWishartDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType>
InverseWishartDistribution<ScalarType>
::InverseWishartDistribution()
{
    m_DegreesOfFreedom = 1;

    m_ScaleMatrix.set_size(1, 1);
    m_ScaleMatrix(0,0) = 1.0;

    m_ScaleMatrixInverse.set_size(1, 1);
    m_ScaleMatrixInverse(0,0) = 1.0;
}


template <class ScalarType>
InverseWishartDistribution<ScalarType>
::~InverseWishartDistribution()
{}


template <class ScalarType>
InverseWishartDistribution<ScalarType>
::InverseWishartDistribution(const InverseWishartDistribution& other)
{
    m_DegreesOfFreedom = other.m_DegreesOfFreedom;
    m_ScaleMatrix = other.m_ScaleMatrix;
    m_ScaleMatrixInverse = other.m_ScaleMatrixInverse;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType>
ScalarType
InverseWishartDistribution<ScalarType>
::ComputeLogLikelihood(LinearVariableType const& obs) const // Input : inverse matrix.
{
    const MatrixType aux = recast<MatrixType>(obs);

    return 0.5 * m_DegreesOfFreedom * (log_det(aux) - trace(aux.transpose() * m_ScaleMatrix));
}

template class InverseWishartDistribution<double>;
