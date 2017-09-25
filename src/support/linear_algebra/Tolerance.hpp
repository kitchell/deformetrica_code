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

namespace def {
namespace algebra {
namespace utils {

#include <cmath>

/**
 * Define a common interface to compare matrix,vector, list of vectors...
 */
class ScalarCompare {

  float tolerance;
 public:

  ScalarCompare(const float t)
      : tolerance(t) {}

  const bool operator()(const float &t1, const float &t2) const {
    return std::fabs((t1 - t2)) <= tolerance;
  }

  const bool operator()(const double &t1, const double &t2) const {
    return std::fabs((float)(t1 - t2)) <= tolerance;
  }

  const bool operator()(const std::vector<float> &t1, const std::vector<float> &t2) const {
    return std::equal(t1.begin(), t1.end(), t2.begin(), ScalarCompare(tolerance));
  }

  const bool operator()(const std::vector<double> &t1, const std::vector<double> &t2) const {
    return std::equal(t1.begin(), t1.end(), t2.begin(), ScalarCompare(tolerance));
  }

};

template <typename T>
class Tolerance {
 public:

  virtual const bool compare(const T&, const float) = 0;

  template<class InputIt1, class InputIt2>
  bool equal(InputIt1 first1, InputIt1 last1,
             InputIt2 first2, const float tolerance) const
  {
    for (; first1 != last1; ++first1, ++first2) {
      if (!(first1->compare(*first2, tolerance))) {
        return false;
      }
    }
    return true;
  }

};

}
}
}