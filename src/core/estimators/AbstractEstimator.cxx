/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "AbstractEstimator.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
AbstractEstimator<ScalarType, Dimension>
::AbstractEstimator() : m_StatisticalModel(0), m_DataSet(0), m_MaxIterations(100), m_CurrentIteration(0),
                             m_PrintEveryNIters(1), m_SaveEveryNIters(100) {}

template<class ScalarType, unsigned int Dimension>
AbstractEstimator<ScalarType, Dimension>
::~AbstractEstimator() {}

template<class ScalarType, unsigned int Dimension>
AbstractEstimator<ScalarType, Dimension>
::AbstractEstimator(const AbstractEstimator &other) {
  m_StatisticalModel = other.m_StatisticalModel->Clone();
  m_DataSet = other.m_DataSet->Clone();
  m_MaxIterations = other.m_MaxIterations;
  m_PrintEveryNIters = other.m_PrintEveryNIters;
  m_SaveEveryNIters = other.m_SaveEveryNIters;
  m_CurrentIteration = other.m_CurrentIteration;
  m_PopulationRER = other.m_PopulationRER;
  m_IndividualRER = other.m_IndividualRER;
}

template
class AbstractEstimator<double, 2>;
template
class AbstractEstimator<double, 3>;
