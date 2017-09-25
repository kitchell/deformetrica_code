/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TimeSeriesDataSet.h"


template <class ScalarType, unsigned int Dimension>
TimeSeriesDataSet<ScalarType, Dimension>
::TimeSeriesDataSet() : Superclass()
{
    Superclass::SetTimeSeriesDataSetType();
}


template <class ScalarType, unsigned int Dimension>
TimeSeriesDataSet<ScalarType, Dimension>
::TimeSeriesDataSet(const TimeSeriesDataSet& other) : Superclass(other)
{}


template <class ScalarType, unsigned int Dimension>
TimeSeriesDataSet<ScalarType, Dimension>
::~TimeSeriesDataSet()
{}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
std::vector<std::shared_ptr<typename TimeSeriesDataSet<ScalarType, Dimension>::DeformableMultiObjectType>>
TimeSeriesDataSet<ScalarType, Dimension>
::GetDeformableMultiObjects() const
{
    const unsigned int nbTimePoints = Superclass::m_NumberOfTimePointsPerSubjects[0];

    std::vector<std::shared_ptr<DeformableMultiObjectType>> obj;
    obj.resize(nbTimePoints);

    for (unsigned int t = 0; t < nbTimePoints; t++)
    {
        obj[t] = Superclass::m_DeformableMultiObjects[0][t];
    }

    return obj;
}


template <class ScalarType, unsigned int Dimension>
void
TimeSeriesDataSet<ScalarType, Dimension>
::SetDeformableMultiObjects(std::vector<std::shared_ptr<DeformableMultiObjectType>> obj)
{
    const unsigned int nbTimePoints = obj.size();

    Superclass::m_DeformableMultiObjects.resize(1);
    Superclass::m_DeformableMultiObjects[0].resize(nbTimePoints);
    for (int i = 0; i < nbTimePoints; ++i)
    {
        Superclass::m_DeformableMultiObjects[0][i] = obj[i];
    }
}


template <class ScalarType, unsigned int Dimension>
void
TimeSeriesDataSet<ScalarType, Dimension>
::SetTimeIndices(std::vector< unsigned int > ti)
{
    Superclass::m_TimeIndices.resize(1);
    Superclass::m_TimeIndices[0] = ti;
}


template <class ScalarType, unsigned int Dimension>
void
TimeSeriesDataSet<ScalarType, Dimension>
::GetDataTimeSeries(std::vector<std::shared_ptr<DeformableMultiObjectType>>& objects, IndicesType& timeIndices) const
{
    objects = Superclass::m_DeformableMultiObjects[0];
    timeIndices = Superclass::m_TimeIndices[0];
}


template <class ScalarType, unsigned int Dimension>
std::shared_ptr<typename TimeSeriesDataSet<ScalarType, Dimension>::DeformableMultiObjectType>
TimeSeriesDataSet<ScalarType, Dimension>
::GetDataForTimeIndex(unsigned int t) const
{
    if (Superclass::m_NumberOfTimePointsPerSubjects[0] <= t)
        throw std::runtime_error("Time index in TimeSeriesDataSet::GetDataForTimeIndex is out of bounds.");

    return Superclass::m_DeformableMultiObjects[0][t];
}


template <class ScalarType, unsigned int Dimension>
void
TimeSeriesDataSet<ScalarType, Dimension>
::SetDataForTimeIndex(std::shared_ptr<DeformableMultiObjectType> const obj, const unsigned int t)
{
    if (Superclass::m_NumberOfTimePointsPerSubjects[0] <= t)
        Superclass::m_DeformableMultiObjects[0].resize(t);

    Superclass::m_DeformableMultiObjects[0][t] = obj;
}


template<class ScalarType, unsigned int Dimension>
void
TimeSeriesDataSet<ScalarType, Dimension>
::Update()
{
    Superclass::Update();

    if (Superclass::m_NumberOfSubjects != 1)
        throw std::runtime_error("Time series data set should have only one subjects");

}

template class TimeSeriesDataSet<double,2>;
template class TimeSeriesDataSet<double,3>;

