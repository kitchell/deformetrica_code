/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

/// Class file.
#include "LongitudinalDataSet.h"

/**
 *	\brief      CrossSectionalDataSet object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A CrossSectional data set is a collection of MultiDeformableObjects for a series of subjects.
 */
template <class ScalarType, unsigned int Dimension>
class CrossSectionalDataSet : public LongitudinalDataSet<ScalarType, Dimension>
{
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// LongitudinalDataSet as Superclass
    typedef LongitudinalDataSet<ScalarType, Dimension> Superclass;

    /// DeformableMultiObjectType
    typedef typename Superclass::DeformableMultiObjectType DeformableMultiObjectType;

    /// Dataset indices type (See AbstractStatisticalModel).
    typedef unsigned int IndicesType;
    /// Dataset residuals type (See AbstractStatisticalModel).
    typedef std::vector<std::vector<ScalarType>> ResidualsType;



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    CrossSectionalDataSet();
    /// Copy constructor.
    CrossSectionalDataSet(const CrossSectionalDataSet& other);

    /// Makes a copy of the object.
    CrossSectionalDataSet* Clone() { return new CrossSectionalDataSet(*this); }

    virtual ~CrossSectionalDataSet();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Gets the deformable multi objects.
    std::vector<std::shared_ptr<DeformableMultiObjectType>> GetDeformableMultiObjects() const;
    /// Sets the deformable multi objects to \e DeformableMultiObjects.
    void SetDeformableMultiObjects(const std::vector<std::shared_ptr<DeformableMultiObjectType>> obj);

    /// Get DeformableMultiObject for a given subject \e s.
    std::shared_ptr<DeformableMultiObjectType> GetDataForSubject(unsigned int s) const
    {
        if (Superclass::m_NumberOfSubjects <= s)
            throw std::runtime_error("Subject's index in CrossSectionalDataSet::GetDataForSubject is out of bounds.");

        return Superclass::m_DeformableMultiObjects[s][0];
    }

    virtual void GetDataForSubject(unsigned int s, std::vector<std::shared_ptr<DeformableMultiObjectType>> objects, std::vector<unsigned int> timeIndices)
    {
        throw std::runtime_error("The method GetDataForSubject(s, objects, time-points) does not make sense for a cross-sectional data set.");
    }



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other public method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Check consistency of vector dimensions, and other sanity checks.
    virtual void Update();

}; /* class CrossSectionalDataSet */

