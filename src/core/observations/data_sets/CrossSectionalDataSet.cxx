/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "CrossSectionalDataSet.h"


template <class ScalarType, unsigned int Dimension>
CrossSectionalDataSet<ScalarType, Dimension>
::CrossSectionalDataSet() : Superclass()
{
    this->SetCrossSectionalDataSetType();
}

template <class ScalarType, unsigned int Dimension>
CrossSectionalDataSet<ScalarType, Dimension>
::CrossSectionalDataSet(const CrossSectionalDataSet& other) : Superclass(other)
{}


template <class ScalarType, unsigned int Dimension>
CrossSectionalDataSet<ScalarType, Dimension>
::~CrossSectionalDataSet()
{}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
std::vector<std::shared_ptr<typename CrossSectionalDataSet<ScalarType, Dimension>::DeformableMultiObjectType>>
CrossSectionalDataSet<ScalarType, Dimension>
::GetDeformableMultiObjects() const
{
    const unsigned int nbSubjects = Superclass::m_NumberOfSubjects;

    std::vector<std::shared_ptr<DeformableMultiObjectType>> obj;
    obj.resize(nbSubjects);

    for (unsigned int s = 0; s < nbSubjects; s++)
    {
        obj[s] = Superclass::m_DeformableMultiObjects[s][0];
    }

    return obj;
}


template <class ScalarType, unsigned int Dimension>
void
CrossSectionalDataSet<ScalarType, Dimension>
::SetDeformableMultiObjects(const std::vector<std::shared_ptr<DeformableMultiObjectType>> obj)
{
    unsigned int nbSubjects = obj.size();

    Superclass::m_DeformableMultiObjects.resize(nbSubjects);
    Superclass::m_NumberOfSubjects = nbSubjects;

    for (int i = 0; i < nbSubjects; ++i)
    {
        Superclass::m_DeformableMultiObjects[i].resize(1);
        Superclass::m_DeformableMultiObjects[i][0] = obj[i];
    }
}


template<class ScalarType, unsigned int Dimension>
void
CrossSectionalDataSet<ScalarType, Dimension>
::Update()
{
    Superclass::Update();
}

template class CrossSectionalDataSet<double,2>;
template class CrossSectionalDataSet<double,3>;

