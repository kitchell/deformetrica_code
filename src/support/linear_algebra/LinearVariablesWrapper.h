/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _LinearVariablesWrapper_h
#define _LinearVariablesWrapper_h


/// Support files.
#include "MatrixListWrapper.h"
#include "LinearVariableWrapper.h"

/// Librairies files.
#include "boost/variant.hpp"
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

template<class ScalarType>
class LinearVariablesWrapper;

template<class ScalarType, typename VariantType>
VariantType recast(std::vector<LinearVariablesWrapper<ScalarType>> const &vec);

/**
 *  \brief      LinearVariablesWrapper object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    // TODO
 *
 */
template<class ScalarType>
class LinearVariablesWrapper : public Tolerance<LinearVariablesWrapper<ScalarType>>  {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Friend methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Converts the vector of linear variable back to a vector of either scalars, vectors, matrices or matrix lists.
  template<typename VariantType>
  friend std::vector<VariantType> recast(LinearVariablesWrapper<ScalarType> const &vec) {
    const unsigned long length = vec.size();
    std::vector<VariantType> result(length);
    for (unsigned int k = 0; k < length; ++k)
      result[k] = recast<VariantType>(vec[k]);
    return result;
  }

 public :

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Linear variable type.
  typedef LinearVariableWrapper<ScalarType> LinearVariableType;

  /// Linear variables type.
  typedef std::vector<LinearVariableType> RawLinearVariablesType;

  /// Iterators.
  typedef typename RawLinearVariablesType::iterator iterator;
  typedef typename RawLinearVariablesType::const_iterator const_iterator;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  inline LinearVariablesWrapper() {}

  /// Copy contructor.
  inline LinearVariablesWrapper(const LinearVariablesWrapper &other) {
    m_RawLinearVariables = other.m_RawLinearVariables;
  }

  /// Constructor from a raw linear variable.
  inline LinearVariablesWrapper(const RawLinearVariablesType &vec) { m_RawLinearVariables = vec; }

  /// Constructor with a size specification.
  inline LinearVariablesWrapper(unsigned int const &n) { m_RawLinearVariables.resize(n); }

  /// Constructor from a raw vector of variant types.
  template<typename VariantType>
  inline LinearVariablesWrapper(std::vector<VariantType> const &vec) {
    const unsigned long length = vec.size();
    m_RawLinearVariables.resize(length);
    for (unsigned int k = 0; k < length; ++k)
      m_RawLinearVariables[k] = vec[k];
  }

  /// Destructor.
  inline ~LinearVariablesWrapper() {}


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Converts the linear variable back to a raw linear variable.
  // Only for LinearVariableMapWrapper use, to be avoided in the rest of Deformetrica.
  inline RawLinearVariablesType const &GetRawLinearVariables() const { return (m_RawLinearVariables); }
  inline RawLinearVariablesType &GetRawLinearVariables() { return (m_RawLinearVariables); }

  /// Returns the length of the map.
  inline unsigned int size() const { return m_RawLinearVariables.size(); }

  /// Returns an iterator on the first linear variable.
  inline iterator begin() { return m_RawLinearVariables.begin(); }
  inline const_iterator begin() const { return m_RawLinearVariables.begin(); }

  /// Returns an iterator on the end of the linear variables.
  inline iterator end() { return m_RawLinearVariables.end(); }
  inline const_iterator end() const { return m_RawLinearVariables.end(); }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the sum all elements.
  inline ScalarType sum() const {
    ScalarType result = 0;
    for (unsigned long k = 0; k < m_RawLinearVariables.size(); ++k)
      result += m_RawLinearVariables[k].sum();
    return result;
  }

  /// Returns the sum of squares of all elements.
  inline ScalarType sum_of_squares() const {
    ScalarType result = 0;
    for (unsigned long k = 0; k < m_RawLinearVariables.size(); ++k)
      result += m_RawLinearVariables[k].sum_of_squares();
    return result;
  }

  /// Returns the mean squares.
  inline ScalarType mean_squares() const {
    const unsigned long length = m_RawLinearVariables.size();
    ScalarType result = 0;
    for (unsigned long k = 0; k < length; ++k)
      result += m_RawLinearVariables[k].sum_of_squares();
    return result / length;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Access the linear variable corresponding to the key \e str for reading.
  inline LinearVariableType const &operator[](unsigned int const &i) const {
    return m_RawLinearVariables[i];
  }
  /// Access the linear variable corresponding to the key \e str for reading or writing.
  inline LinearVariableType &operator[](unsigned int const &i) {
    return m_RawLinearVariables[i];
  }

  const bool operator==(const LinearVariablesWrapper<ScalarType> &t) const {
    return m_RawLinearVariables.size() == t.m_RawLinearVariables.size() &&
        std::equal(m_RawLinearVariables.begin(), m_RawLinearVariables.end(), t.m_RawLinearVariables.begin());
  }

  virtual const bool compare(const LinearVariablesWrapper<ScalarType>& t, const float tolerance) {
    return m_RawLinearVariables.size() == t.m_RawLinearVariables.size() &&
        this->equal(m_RawLinearVariables.begin(), m_RawLinearVariables.end(), t.m_RawLinearVariables.begin(), tolerance);
  }


 private :
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & m_RawLinearVariables;
  }
  friend class boost::serialization::access;

  /// Raw linear variable.
  RawLinearVariablesType m_RawLinearVariables;

};

template<class ScalarType>
inline LinearVariablesWrapper<ScalarType> operator+(LinearVariablesWrapper<ScalarType> const &left,
                                                    LinearVariablesWrapper<ScalarType> const &right) {
  const unsigned int s = left.size();
  LinearVariablesWrapper<ScalarType> out(s);
  for (unsigned int k = 0; k < s; ++k) {
    out[k] = left[k] + right[k];
  }
  return out;
}

template<class ScalarType>
inline LinearVariablesWrapper<ScalarType> operator-(LinearVariablesWrapper<ScalarType> const &left,
                                                    LinearVariablesWrapper<ScalarType> const &right) {
  const unsigned int s = left.size();
  LinearVariablesWrapper<ScalarType> out(s);
  for (unsigned int k = 0; k < s; ++k) {
    out[k] = left[k] - right[k];
  }
  return out;
}

template<class ScalarType>
inline LinearVariablesWrapper<ScalarType> operator*(LinearVariablesWrapper<ScalarType> const &left,
                                                    ScalarType const &right) {
  const unsigned int s = left.size();
  LinearVariablesWrapper<ScalarType> out(s);
  for (unsigned int k = 0; k < s; ++k) {
    out[k] = left[k] * right;
  }
  return out;
}

template<class ScalarType>
inline LinearVariablesWrapper<ScalarType> operator*(ScalarType const &left,
                                                    LinearVariablesWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariablesWrapper<ScalarType> operator/(LinearVariablesWrapper<ScalarType> const &left,
                                                    ScalarType const &right) {
  return operator*(left, 1 / right);
}

template<class ScalarType>
inline ScalarType dot_product(LinearVariablesWrapper<ScalarType> const &left,
                              LinearVariablesWrapper<ScalarType> const &right) {
  ScalarType result = 0;
  for (unsigned int k = 0; k < left.size(); k++)
    result += dot_product(left[k], right[k]);
  return result;
}

#endif /* _LinearVariablesWrapper_h */
