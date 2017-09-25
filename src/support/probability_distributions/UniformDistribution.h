/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _UniformDistribution_h
#define _UniformDistribution_h

/// Class file.
#include "AbstractProbabilityDistribution.h"


/**
 *  \brief      UniformDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a normal distribution, determined by its mean and variance.
 */
template <class ScalarType>
class UniformDistribution : public AbstractProbabilityDistribution<ScalarType>
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
    UniformDistribution(); // Standard uniform distribution.
    UniformDistribution(VectorType const& lowers, VectorType const& uppers);

    /// Copy constructor.
    UniformDistribution(const UniformDistribution& other);

    /// Makes a copy of the object.
    virtual UniformDistribution* Clone() { return new UniformDistribution(*this); }

    /// Destructor.
    virtual ~UniformDistribution();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other public method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Samples from the normal distribution to generate a realization.
    virtual VectorType Sample() const;

    /// Computes the log-likelihood of an observation.
    virtual ScalarType ComputeLogLikelihood(LinearVariableType const& obs) const;


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Protected attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Bounds of the uniform distribution.
    VectorType m_LowerBounds;
    VectorType m_UpperBounds;

}; /* class UniformDistribution */

#endif /* _UniformDistribution_h */
