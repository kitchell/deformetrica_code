/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _InverseWishartDistribution_h
#define _InverseWishartDistribution_h

/// Class file.
#include "AbstractProbabilityDistribution.h"


/**
 *  \brief      InverseWishartDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a normal distribution, determined by its mean and variance.
 */
template <class ScalarType>
class InverseWishartDistribution : public AbstractProbabilityDistribution<ScalarType>
{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Abstract probability distribution type.
    typedef AbstractProbabilityDistribution<ScalarType> Superclass;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Constructors.
    InverseWishartDistribution();

    /// Copy constructor.
    InverseWishartDistribution(const InverseWishartDistribution& other);

    /// Makes a copy of the object.
    virtual InverseWishartDistribution* Clone() { return new InverseWishartDistribution(*this); }

    /// Destructor.
    virtual ~InverseWishartDistribution();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Returns the degrees of freedom of the inverse wishart distribution.
    ScalarType GetDegreesOfFreedom() const { return m_DegreesOfFreedom; }
    /// Sets the degrees of freedom of the inverse wishart distribution.
    void SetDegreesOfFreedom(ScalarType const& a) { m_DegreesOfFreedom = a; }


    /// Returns the scale matrix of the inverse wishart distribution.
    MatrixType GetScaleMatrix() const { return m_ScaleMatrix; }
    /// Sets the scale matrix of the inverse wishart distribution.
    void SetScaleMatrix(MatrixType const& psi) {
        m_ScaleMatrix = psi;
        m_ScaleMatrixInverse = inverse_sympd(psi); }

    /// Returns the inverse of the scale matrix.
    MatrixType GetScaleMatrixInverse() const { return m_ScaleMatrixInverse; }
    /// Sets the inverse scale matrix.
    void SetScaleMatrixInverse(MatrixType const& psi) {
        m_ScaleMatrixInverse = psi;
        m_ScaleMatrix = inverse_sympd(psi); }



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other public method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Samples from the inverse wishart distribution to generate a realization. TODO.
    virtual VectorType Sample() const {
        std::cout << "Warning : in InverseWishartDistribution, the Sample() method is not implemented yet." << std::endl;
        return VectorType(1, 0.0); }

    /// Computes the log-likelihood of an observation. Input : inverse matrix.
    virtual ScalarType ComputeLogLikelihood(LinearVariableType const& obs) const;


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Protected attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Degrees of freedom of the inverse wishart distribution.
    ScalarType m_DegreesOfFreedom;

    /// Scale matrix of the inverse wishart distribution.
    MatrixType m_ScaleMatrix;

    /// Inverse of the scale matrix.
    MatrixType m_ScaleMatrixInverse;



}; /* class InverseWishartDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "InverseWishartDistribution.txx"
//#endif


#endif /* _InverseWishartDistribution_h */
