/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _LinearVariablesMapWrapper_h
#define _LinearVariablesMapWrapper_h

/// Non-core files.
#include "ArmadilloVectorWrapper.h"
#include "ArmadilloMatrixWrapper.h"
#include "MatrixListWrapper.h"
#include "LinearVariablesWrapper.h"

/// Librairies files.
#include "boost/variant.hpp"
#include <map>
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>

template<class ScalarType>
class LinearVariablesMapWrapper;

template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator+(const LinearVariablesMapWrapper<ScalarType> &left,
                                                const LinearVariablesMapWrapper<ScalarType> &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator-(const LinearVariablesMapWrapper<ScalarType> &left,
                                                const LinearVariablesMapWrapper<ScalarType> &right);

template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                ArmadilloVectorWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(ArmadilloVectorWrapper<ScalarType> const &left,
                                                LinearVariablesMapWrapper<ScalarType> const &right);

template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                ScalarType const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(ScalarType const &left,
                                                LinearVariablesMapWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator/(const LinearVariablesMapWrapper<ScalarType> &left,
                                                ScalarType const &right);

template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                int const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(int const &left,
                                                LinearVariablesMapWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator/(const LinearVariablesMapWrapper<ScalarType> &left,
                                                int const &right);

template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                unsigned int const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator*(unsigned int const &left,
                                                LinearVariablesMapWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariablesMapWrapper<ScalarType> operator/(const LinearVariablesMapWrapper<ScalarType> &left,
                                                unsigned int const &right);

/**
 *  \brief      LinearVariablesMapWrapper object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    // TODO
 *
 */
template<class ScalarType>
class LinearVariablesMapWrapper : public Tolerance<LinearVariablesMapWrapper<ScalarType>> {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Friend methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  friend LinearVariablesMapWrapper<ScalarType> operator-<>(LinearVariablesMapWrapper<ScalarType> const &left,
                                                           LinearVariablesMapWrapper<ScalarType> const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator+<>(LinearVariablesMapWrapper<ScalarType> const &left,
                                                           LinearVariablesMapWrapper<ScalarType> const &right);

  friend LinearVariablesMapWrapper<ScalarType> operator*<>(LinearVariablesMapWrapper<ScalarType> const &left,
                                                           ArmadilloVectorWrapper<ScalarType> const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator*<>(ArmadilloVectorWrapper<ScalarType> const &left,
                                                           LinearVariablesMapWrapper<ScalarType> const &right);

  friend LinearVariablesMapWrapper<ScalarType> operator*<>(LinearVariablesMapWrapper<ScalarType> const &left,
                                                           ScalarType const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator*<>(ScalarType const &left,
                                                           LinearVariablesMapWrapper<ScalarType> const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator/<>(const LinearVariablesMapWrapper<ScalarType> &left,
                                                           ScalarType const &right);

  friend LinearVariablesMapWrapper<ScalarType> operator*<>(LinearVariablesMapWrapper<ScalarType> const &left,
                                                           int const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator*<>(int const &left,
                                                           LinearVariablesMapWrapper<ScalarType> const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator/<>(const LinearVariablesMapWrapper<ScalarType> &left,
                                                           int const &right);

  friend LinearVariablesMapWrapper<ScalarType> operator*<>(LinearVariablesMapWrapper<ScalarType> const &left,
                                                           unsigned int const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator*<>(unsigned int const &left,
                                                           LinearVariablesMapWrapper<ScalarType> const &right);
  friend LinearVariablesMapWrapper<ScalarType> operator/<>(const LinearVariablesMapWrapper<ScalarType> &left,
                                                           unsigned int const &right);

 public :

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Vector type.
  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  /// Matrix type.
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  /// List of matrices type.
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Linear variable type.
  typedef LinearVariablesWrapper<ScalarType> LinearVariablesType;

  /// Raw linear variable map type.
  typedef std::map<std::string, LinearVariablesType> RawLinearVariablesMapType;
  /// Iterator on a raw linear variable map type.
  typedef typename RawLinearVariablesMapType::const_iterator const_iterator;
  typedef typename RawLinearVariablesMapType::iterator iterator;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  LinearVariablesMapWrapper() {}

  /// Constructor from a raw linear variable map.
  LinearVariablesMapWrapper(const RawLinearVariablesMapType &map) { m_RawLinearVariablesMap = map; }

  /// Copy contructor.
  LinearVariablesMapWrapper(const LinearVariablesMapWrapper &other) {
    m_RawLinearVariablesMap.clear();
    m_RawLinearVariablesMap = other.m_RawLinearVariablesMap;
  }

  /// Destructor.
  ~LinearVariablesMapWrapper() {}


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Converts the linear variable map back to a raw linear variable map.
  RawLinearVariablesMapType GetRawLinearVariablesMap() const { return m_RawLinearVariablesMap; }

  /// Returns the length of the map.
  unsigned int size() const { return m_RawLinearVariablesMap.size(); }
//    /// Resizes the map to a length \e n.
//    void resize(const unsigned int n) { m_RawLinearVariablesMap.resize(n); }

  /// Fills the linear variables map.
  void fill(ScalarType const &value) {
    for (auto it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        it->second[k].fill(value);
      }
    }
  }

  /// Returns an iterator on the first element of the raw linear variable map.
  const_iterator begin() const noexcept { return m_RawLinearVariablesMap.begin(); }
  iterator begin()       noexcept { return m_RawLinearVariablesMap.begin(); }

  /// Returns an iterator on the end of the raw linear variable map.
  const_iterator end() const noexcept { return m_RawLinearVariablesMap.end(); }
  iterator end()       noexcept { return m_RawLinearVariablesMap.end(); }

  /// Returns an iterator on the searched element, or on the map end if the key \e str does not exist.
  const_iterator find(std::string const &str) const { return m_RawLinearVariablesMap.find(str); }
  iterator find(std::string const &str) { return m_RawLinearVariablesMap.find(str); }

  /// Counts the number of elements with a specific key \e str.
  unsigned int count(std::string const &str) const { return m_RawLinearVariablesMap.count(str); }

  /// Erase a specific entry.
  void erase(const std::string &key) { m_RawLinearVariablesMap.erase(key); }

  /// Clears the map, removing all elements.
  void clear() { m_RawLinearVariablesMap.clear(); }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the sum of all elements.
  ScalarType sum() const {
    ScalarType result(0);
    for (const_iterator it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        result += it->second[k].sum();
      }
    }
    return result;
  }

  /// Returns sum of squares of all elements.
  ScalarType sum_of_squares() const {
    ScalarType result(0);
    for (const_iterator it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        result += it->second[k].sum_of_squares();
      }
    }
    return result;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Access the linear variable corresponding to the key \e str for reading.
  LinearVariablesType const &operator[](std::string const &str) const { return m_RawLinearVariablesMap.at(str); }
  LinearVariablesType const &at(std::string const &str) const { return m_RawLinearVariablesMap.at(str); }
  /// Access the linear variable corresponding to the key \e str for reading or writing.
  LinearVariablesType &operator[](std::string const &str) { return m_RawLinearVariablesMap[str]; }

  /// Add rhs scalar to lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator+=(ScalarType const &rhsScalar) {
    for (auto it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        it->second[k] += rhsScalar;
      }
    }
    return *this;
  }
  /// Add rhs int to lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator+=(int const &rhsInt) {
    return operator+=(ScalarType(rhsInt));
  }
  /// Add rhs unsigned int to lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator+=(unsigned int const &rhsInt) {
    return operator+=(ScalarType(rhsInt));
  }

  /// Subtract rhs scalar from lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator-=(ScalarType const &rhsScalar) {
    return operator+=(-rhsScalar);
  }
  /// Subtract rhs int from lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator-=(int const &rhsInt) {
    return operator+=(-rhsInt);
  }
  /// Subtract rhs unsigned int from lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator-=(unsigned int const &rhsInt) {
    return operator+=(-rhsInt);
  }

  /// Scalar multiplication in situ of lhs linear variables map by \e rhsScalar.
  LinearVariablesMapWrapper<ScalarType> &operator*=(ScalarType const &rhsScalar) {
    for (auto it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        it->second[k] *= rhsScalar;
      }
    }
    return *this;
  }
  /// Int multiplication in situ of lhs linear variables map by \e rhsInt.
  LinearVariablesMapWrapper<ScalarType> &operator*=(int const &rhsInt) {
    return operator*=(ScalarType(rhsInt));
  }
  /// Unsigned int multiplication in situ of lhs linear variables map by \e rhsInt.
  LinearVariablesMapWrapper<ScalarType> &operator*=(unsigned int const &rhsInt) {
    return operator*=(ScalarType(rhsInt));
  }

  /// Scalar division of lhs linear variables map in situ by \e rhsScalar.
  LinearVariablesMapWrapper<ScalarType> &operator/=(ScalarType const &rhsScalar) {
    return operator*=(1 / rhsScalar);
  }
  /// Int division of lhs linear variables map in situ by \e rhsInt.
  LinearVariablesMapWrapper<ScalarType> &operator/=(int const &rhsInt) {
    return operator*=(1.0 / ScalarType(rhsInt));
  }
  /// Unsigned int division of lhs linear variables map in situ by \e rhsInt.
  LinearVariablesMapWrapper<ScalarType> &operator/=(unsigned int const &rhsInt) {
    return operator*=(1.0 / ScalarType(rhsInt));
  }

  /// Add rhs linear variables map to lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator+=(LinearVariablesMapWrapper<ScalarType> const &right) {
    for (auto it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        it->second[k] += right.m_RawLinearVariablesMap.at(it->first)[k];
      }
    }
    return *this;
  }
  /// Subtract rhs linear variable map from lhs linear variables map in situ.
  LinearVariablesMapWrapper<ScalarType> &operator-=(LinearVariablesMapWrapper<ScalarType> const &right) {
    for (auto it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        it->second[k] -= right.m_RawLinearVariablesMap.at(it->first)[k];
      }
    }
    return *this;
  }
  /// Unary minus operator.
  LinearVariablesMapWrapper<ScalarType> operator-() const {
    RawLinearVariablesMapType rawList;
    for (auto it = m_RawLinearVariablesMap.begin(); it != m_RawLinearVariablesMap.end(); ++it) {
      for (unsigned int k = 0; k < it->second.size(); ++k) {
        rawList[it->first] = -it->second[k];
      }
    }
    return LinearVariablesMapWrapper<ScalarType>(rawList);
  }
  /// Equality testing operator.
  const bool operator==(const LinearVariablesMapWrapper<ScalarType> &t) const {
    return m_RawLinearVariablesMap.size() == t.m_RawLinearVariablesMap.size() &&
        std::equal(m_RawLinearVariablesMap.begin(), m_RawLinearVariablesMap.end(), t.m_RawLinearVariablesMap.begin());
  }

  virtual const bool compare(const LinearVariablesMapWrapper<ScalarType>& t, const float tolerance) {
    if (m_RawLinearVariablesMap.size() != t.m_RawLinearVariablesMap.size()) return false;
    auto it1 = m_RawLinearVariablesMap.begin();
    auto it2 = t.m_RawLinearVariablesMap.begin();
    while(it1 != m_RawLinearVariablesMap.end()) {
      if (it1->first != it2->first) return false;
      if (!it1->second.compare(it2->second, tolerance)) return false;
      it1++;
      it2++;
    }

    return true;
  }

 private :
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & m_RawLinearVariablesMap;
  }
  friend class boost::serialization::access;

  /// Raw linear variable map.
  RawLinearVariablesMapType m_RawLinearVariablesMap;

};

#endif /* _LinearVariablesMapWrapper_h */
