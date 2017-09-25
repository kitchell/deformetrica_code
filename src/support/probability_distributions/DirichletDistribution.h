/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DirichletDistribution_h
#define _DirichletDistribution_h

/// Class file.
#include "AbstractProbabilityDistribution.h"


/**
 *  \brief      DirichletDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a Dirichlet distribution.
 */
template <class ScalarType>
class DirichletDistribution : public AbstractProbabilityDistribution<ScalarType>
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
    DirichletDistribution();

    /// Copy constructor.
    DirichletDistribution(const DirichletDistribution& other);

    /// Makes a copy of the object.
    virtual DirichletDistribution* Clone() { return new DirichletDistribution(*this); }

    /// Destructor.
    virtual ~DirichletDistribution();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Returns the concentration parameters of the Dirichlet distribution.
    VectorType GetConcentrationParameters() const { return m_ConcentrationParameters; }
    /// Sets the concentration parameters of the Dirichlet distribution, for different input types.
    void SetConcentrationParameters(VectorType const& cp) { m_ConcentrationParameters = cp; }
    void SetConcentrationParameters(LinearVariableType const& cp) { m_ConcentrationParameters = vectorize(cp); }


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other public method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Samples from the Dirichlet distribution to generate a realization. TODO.
    virtual VectorType Sample() const {
        std::cout << "Exception. DirichletDistribution::Sample() method is not available yet." << std::endl;
        return VectorType(0, 0.0); }

    /// Computes the log-likelihood of an observation. Vectorizes the linear variable in a first step.
    virtual ScalarType ComputeLogLikelihood(LinearVariableType const& obs) const;
    /// Computes the log-likelihood of an observation.
    ScalarType ComputeLogLikelihood(VectorType const& obs) const;


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Protected attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Concentration parameters of the Dirichlet distribution.
    VectorType m_ConcentrationParameters;


}; /* class DirichletDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "DirichletDistribution.txx"
//#endif


#endif /* _DirichletDistribution_h */
