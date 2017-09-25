/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "LongitudinalDataSet.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
LongitudinalDataSet<ScalarType, Dimension>
::LongitudinalDataSet() {
  m_Type = Longitudinal;
}

template<class ScalarType, unsigned int Dimension>
LongitudinalDataSet<ScalarType, Dimension>
::LongitudinalDataSet(const LongitudinalDataSet &other) {
  m_Type = other.m_Type;
}

template<class ScalarType, unsigned int Dimension>
LongitudinalDataSet<ScalarType, Dimension>
::~LongitudinalDataSet() {}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalDataSet<ScalarType, Dimension>
::GetDataForSubject(unsigned int s,
                    std::vector<std::shared_ptr<DeformableMultiObjectType>> &objects,
                    std::vector<unsigned int> &timeIndices) const {
  if (m_NumberOfSubjects <= s)
    throw std::runtime_error("subject's index in LongitudinalDataSet::GetDataForSubject out of bound");

  objects = m_DeformableMultiObjects[s];
  timeIndices = m_TimeIndices[s];
}

template<class ScalarType, unsigned int Dimension>
void
LongitudinalDataSet<ScalarType, Dimension>
::GetDataForSubject(unsigned int s,
                    std::vector<std::shared_ptr<DeformableMultiObjectType>> &objects,
                    std::vector<ScalarType> &times) const {
  if (m_NumberOfSubjects <= s)
    throw std::runtime_error("subject's index in LongitudinalDataSet::GetDataForSubject out of bound");

  objects = m_DeformableMultiObjects[s];
  times = m_Times[s];
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
LongitudinalDataSet<ScalarType, Dimension>
::Update() {
  m_NumberOfSubjects = m_DeformableMultiObjects.size();
  if (m_NumberOfSubjects == 0)
    throw std::runtime_error("No DeformationMultiObjects given in LongitudinalDataSet");

  m_NumberOfObjects = m_DeformableMultiObjects[0][0]->GetNumberOfObjects();

  m_NumberOfTimePointsPerSubjects.resize(m_NumberOfSubjects);
  unsigned int totalNumberOfObservations = 0;
  for (unsigned int i = 0; i < m_NumberOfSubjects; ++i) {
    unsigned int nbTimePoints = m_DeformableMultiObjects[i].size();
    totalNumberOfObservations += nbTimePoints;
    m_NumberOfTimePointsPerSubjects[i] = nbTimePoints;

    /// Checks if the number of data in the subject time series matches the number of time indices given.
    if ((m_Type == TimeSeries) && nbTimePoints != m_TimeIndices[i].size())
      throw std::runtime_error("Number of time points and number of multiobjects mismatch in TimeSeries dataset.");

    /// Similar check in the longitudinal case.
    if ((m_Type == Longitudinal) && nbTimePoints != m_Times[i].size())
      throw std::runtime_error("Number of time points and number of multiobjects mismatch in TimeSeries dataset.");

    for (unsigned int j = 0; j < nbTimePoints; j++) {
      /// Check if the number of objects is identical in every DeformableMultiObject in the data set
      if (m_DeformableMultiObjects[i][j]->GetNumberOfObjects() != m_NumberOfObjects)
        throw std::runtime_error("Number of objects mismatch in Longitudinal Data Set.");
    }
  }
  m_TotalNumberOfObservations = totalNumberOfObservations;
}

template class LongitudinalDataSet<double,2>;
template class LongitudinalDataSet<double,3>;
