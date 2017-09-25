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

/// Core files.
#include "DeformableMultiObject.h"

/// Non-core files.
#include "LinearAlgebra.h"

/// Standard files.
#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>

/**
 *	\brief      LongitudinalDataSet object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A longitudinal data set is a collection of MultiDeformableObjects for a series of subjects at multiple time-points.
 *              If there is only one subject and several time-points, the data set is called a time-series data set (see child class TimeSeriesDataSet.h).
 *              If there is several subjects with only one time-point, the data set is called a cross-sectional data set (see child class CrossSectionalDataSet.h).
 */
template<class ScalarType, unsigned int Dimension>
class LongitudinalDataSet {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Possible type of Atlas.
  typedef enum {
    null,            /*!< Null value. */
    Longitudinal,    /*!< Longitudinal Data Set (the most general case) */
    TimeSeries,      /*!< Time series data set (only one subject, several time-points */
    CrossSectional,  /*!< Cross-sectional data set (several subjects, only one time-point */
  } DataSetType;

  /// Multi-Deformable object type.
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  LongitudinalDataSet();

  /// Copy constructor.
  LongitudinalDataSet(const LongitudinalDataSet &other);

  /// Makes a copy of the object.
  LongitudinalDataSet *Clone() { return new LongitudinalDataSet(*this); }

  /// Destructor.
  virtual ~LongitudinalDataSet();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns true if the data set is longitudinal.
  bool IsLongitudinal() { return (m_Type == Longitudinal); }
  /// Returns true if data set is of TimeSeries Data Set Type.
  bool IsTimeSeriesDataSetType() { return (m_Type == TimeSeries); }
  /// Sets TimeSeries Data Set Type.
  void SetTimeSeriesDataSetType() { m_Type = TimeSeries; }
  /// Returns true if data set is of CrossSectional Data Set Type.
  bool IsCrossSectionalDataSetType() { return (m_Type == CrossSectional); }
  /// Sets CrossSectional Data Set Type.
  void SetCrossSectionalDataSetType() { m_Type = CrossSectional; }

  /// Gets the number of subjects.
  unsigned int GetNumberOfSubjects() const { return m_NumberOfSubjects; }

  /// Gets the total number of observations.
  unsigned int GetTotalNumberOfObservations() const { return m_TotalNumberOfObservations; }

  /// Gets the deformable multi objects.
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> GetDeformableMultiObjects() const {
    return m_DeformableMultiObjects;
  }
  /// Sets the deformable multi objects to \e DeformableMultiObjects.
  void SetDeformableMultiObjects(const std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> obj) {
    m_DeformableMultiObjects = obj;
  }

  /// Sets the time indices corresponding to each observation of each subject.
  void SetTimeIndices(const std::vector<std::vector<unsigned int>> ti) { m_TimeIndices = ti; }

  /// Gets the times corresponding to each observation of each subject.
  std::vector<std::vector<ScalarType>> GetTimes() const { return m_Times; }
  /// Sets the times corresponding to each observation of each subject.
  void SetTimes(const std::vector<std::vector<ScalarType>> &times) { m_Times = times; }

  /// Gets the subject ids list.
  std::vector<std::string> GetSubjectIds() const { return m_SubjectIds; }
  /// Sets the subject ids list.
  void SetSubjectIds(const std::vector<std::string> &ids) { m_SubjectIds = ids; }

  /// Get time series DeformableMultiObject and associated time indices for a given subject \e s.
  virtual void GetDataForSubject(unsigned int s,
                                 std::vector<std::shared_ptr<DeformableMultiObjectType>> &objects,
                                 std::vector<unsigned int> &timeIndices) const;
  /// Get time series DeformableMultiObject and associated times for a given subject \e s.
  void GetDataForSubject(unsigned int s,
                         std::vector<std::shared_ptr<DeformableMultiObjectType>> &objects,
                         std::vector<ScalarType> &times) const;

  /// Getter for the classes.
  std::vector<unsigned int> GetClasses() const { return m_Classes; }
  /// Setter for the classes.
  void SetClasses(std::vector<unsigned int> classes) { m_Classes = classes; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Check consistency of vector dimensions, and other sanity checks
  virtual void Update();

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected Method(s):
  ////////////////////////////////////////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s):
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Type of data set (Longitudinal, CrossSectional, TimeSeries).
  DataSetType m_Type;

  /// m_DeformableMultiObjects[i][j]: the jth data of the ith subject.
  std::vector<std::vector<std::shared_ptr<DeformableMultiObjectType>>> m_DeformableMultiObjects;
  /// m_TimeIndices[i][j] the time index corresponding of the jth data of the ith subject.
  std::vector<std::vector<unsigned int>> m_TimeIndices;
  /// m_Times[i][j] the time corresponding of the jth data of the ith subject.
  std::vector<std::vector<ScalarType>> m_Times;

  /// An ordered vector containing the id (string) of each subject.
  std::vector<std::string> m_SubjectIds;

  /// Number of subjects.
  unsigned int m_NumberOfSubjects;
  /// Total number of observation.
  unsigned int m_TotalNumberOfObservations;
  /// Number of objects: should be identical for each data regardless of the subject and time-point.
  unsigned int m_NumberOfObjects;

  /// Number of class for each subject, for supervised algorithms.
  std::vector<unsigned int> m_Classes;

  /// m_NumberOfTimePointsPerSubjects[i]: number of time points for the ith subject
  std::vector<unsigned int> m_NumberOfTimePointsPerSubjects;

}; /* class LongitudinalDataSet */

