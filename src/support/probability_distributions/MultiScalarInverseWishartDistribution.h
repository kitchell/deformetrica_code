/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _MultiScalarInverseWishartDistribution_h
#define _MultiScalarInverseWishartDistribution_h

/// Class file.
#include "AbstractProbabilityDistribution.h"


/**
 *  \brief      MultiScalarInverseWishartDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a normal distribution, determined by its mean and variance.
 */
template <class ScalarType>
class MultiScalarInverseWishartDistribution : public AbstractProbabilityDistribution<ScalarType>
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
    MultiScalarInverseWishartDistribution();

    /// Copy constructor.
    MultiScalarInverseWishartDistribution(const MultiScalarInverseWishartDistribution& other);

    /// Makes a copy of the object.
    virtual MultiScalarInverseWishartDistribution* Clone() { return new MultiScalarInverseWishartDistribution(*this); }

    /// Destructor.
    virtual ~MultiScalarInverseWishartDistribution();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Returns the degrees of freedom of the inverse wishart distribution.
    VectorType GetDegreesOfFreedom() const { return m_DegreesOfFreedom; }
    /// Sets the degrees of freedom of the inverse wishart distribution.
    void SetDegreesOfFreedom(VectorType const& a) { m_DegreesOfFreedom = a; }

    /// Returns the scale vector of the multi scalar inverse wishart distribution.
    VectorType GetScaleVector() const { return m_ScaleVector; }
    /// Sets the scale vector of the multi scalar inverse wishart distribution.
    void SetScaleVector(VectorType const& psi) { m_ScaleVector = psi; }



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other public method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Samples from the multi scalar inverse wishart distribution to generate a realization. TODO.
    virtual VectorType Sample() const {
        std::cout << "Warning : in MultiScalarInverseWishartDistribution, the Sample() method is not "
            "implemented yet." << std::endl;
        return VectorType(1, 0.0); }

    /// Computes the log-likelihood of an observation. Input : variances.
    virtual ScalarType ComputeLogLikelihood(LinearVariableType const& obs) const;


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Protected attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Degrees of freedom of the multi scalar inverse wishart distribution.
    VectorType m_DegreesOfFreedom;

    /// Scale vector of the multi scalar inverse wishart distribution.
    VectorType m_ScaleVector;



}; /* class MultiScalarInverseWishartDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "MultiScalarInverseWishartDistribution.txx"
//#endif


#endif /* _MultiScalarInverseWishartDistribution_h */
