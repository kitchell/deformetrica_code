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
 *	\brief      TimeSeriesDataSet object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A TimeSeries data set is a collection of MultiDeformableObjects for a series of subjects.
 */
template<class ScalarType, unsigned int Dimension>
class TimeSeriesDataSet : public LongitudinalDataSet<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// LongitudinalDataSet as Superclass
  typedef LongitudinalDataSet<ScalarType, Dimension> Superclass;

  /// DeformableMultiObjectType
  typedef typename Superclass::DeformableMultiObjectType DeformableMultiObjectType;

  /// Dataset indices type (See AbstractStatisticalModel).
  typedef std::vector<unsigned int> IndicesType;
  /// Dataset residuals type (See AbstractStatisticalModel).
  typedef std::vector<std::vector<ScalarType>> ResidualsType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TimeSeriesDataSet();
  /// Copy constructor.
  TimeSeriesDataSet(const TimeSeriesDataSet &other);

  /// Makes a copy of the object.
  TimeSeriesDataSet *Clone() { return new TimeSeriesDataSet(*this); }

  virtual ~TimeSeriesDataSet();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Gets the deformable multi objects.
  std::vector<std::shared_ptr<DeformableMultiObjectType>> GetDeformableMultiObjects() const;
  /// Sets the deformable multi objects to \e DeformableMultiObjects.
  void SetDeformableMultiObjects(const std::vector<std::shared_ptr<DeformableMultiObjectType>> obj);

  /// Gets timeIndices corresponding to each observation in the time series
  std::vector<unsigned int> GetTimeIndices() const { return Superclass::m_TimeIndices[0]; }
  /// Sets timeIndices corresponding to each observation in the time series
  void SetTimeIndices(const std::vector<unsigned int> ti);

  /// Gets the subject times.
  std::vector<ScalarType> GetTimes() const { return Superclass::m_Times[0]; }

  /// Gets the subject id.
  std::string GetSubjectId() const { return Superclass::m_SubjectIds[0]; }

  /// Gets the time series data in the appropriate format.
  void GetDataTimeSeries(std::vector<std::shared_ptr<DeformableMultiObjectType>> &objects,
                         IndicesType &timeIndices) const;

  /// Gets DeformableMultiObject for a given time point \e t.
  std::shared_ptr<DeformableMultiObjectType> GetDataForTimeIndex(unsigned int t) const;
  /// Sets DeformableMultiObject for a given time point /e t/
  void SetDataForTimeIndex(std::shared_ptr<DeformableMultiObjectType> const obj, const unsigned int t);

  /// Overwrites the parent method, unfit in the particular case of a time series dataset.
  virtual void GetDataForSubject(unsigned int s,
                                 std::vector<std::shared_ptr<DeformableMultiObjectType>> objects,
                                 std::vector<unsigned int> timeIndices) {
    throw std::runtime_error("the method GetDataForSubject(s, objects, time-points does not make sense "
                                 "for a cross-sectional data set");
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Check consistency of vector dimensions, and other sanity checks
  virtual void Update();

}; /* class TimeSeriesDataSet */
