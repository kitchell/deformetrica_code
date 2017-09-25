/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "AbstractStatisticalModel.h"

/// Initialization of the static count of un-named instances
template<class ScalarType, unsigned int Dimension>
int AbstractStatisticalModel<ScalarType, Dimension>
    ::count = 0;


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
AbstractStatisticalModel<ScalarType, Dimension>
::AbstractStatisticalModel() : m_Type(null) {
  count++;
  std::ostringstream name;
  name << "Model_";
  name << count;
  m_Name = name.str();

}

template<class ScalarType, unsigned int Dimension>
AbstractStatisticalModel<ScalarType, Dimension>
::~AbstractStatisticalModel() {}

template<class ScalarType, unsigned int Dimension>
AbstractStatisticalModel<ScalarType, Dimension>
::AbstractStatisticalModel(const AbstractStatisticalModel &other) {
  count++;
  std::ostringstream name;
  name << "Model ";
  name << count;
  m_Name = name.str();

  m_Type = other.m_Type;
}

template
class AbstractStatisticalModel<double, 2>;
template
class AbstractStatisticalModel<double, 3>;