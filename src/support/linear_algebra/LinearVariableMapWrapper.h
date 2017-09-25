/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _LinearVariableMapWrapper_h
#define _LinearVariableMapWrapper_h


/// Core files.
#include "LinearVariableVariantVisitors.h"

/// Non-core files.
#include "MatrixListWrapper.h"
#include "LinearVariableWrapper.h"

/// Librairies files.
#include "boost/variant.hpp"
#include <map>
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>

#include "Tolerance.hpp"
using namespace def::algebra::utils;

template<class ScalarType>
class LinearVariableMapWrapper;

template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator+(const LinearVariableMapWrapper<ScalarType> &left,
                                               const LinearVariableMapWrapper<ScalarType> &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator-(const LinearVariableMapWrapper<ScalarType> &left,
                                               const LinearVariableMapWrapper<ScalarType> &right);

template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                               ArmadilloVectorWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(ArmadilloVectorWrapper<ScalarType> const &left,
                                               LinearVariableMapWrapper<ScalarType> const &right);

template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                               ScalarType const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(ScalarType const &left,
                                               LinearVariableMapWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator/(const LinearVariableMapWrapper<ScalarType> &left,
                                               ScalarType const &right);

template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                               int const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(int const &left,
                                               LinearVariableMapWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator/(const LinearVariableMapWrapper<ScalarType> &left,
                                               int const &right);

template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                               unsigned int const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator*(unsigned int const &left,
                                               LinearVariableMapWrapper<ScalarType> const &right);
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> operator/(const LinearVariableMapWrapper<ScalarType> &left,
                                               unsigned int const &right);

/**
 *  \brief      LinearVariableMapWrapper object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    // TODO
 *
 */
template<class ScalarType>
class LinearVariableMapWrapper : public Tolerance<LinearVariableMapWrapper<ScalarType>> {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Friend methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  friend LinearVariableMapWrapper<ScalarType> operator-<>(LinearVariableMapWrapper<ScalarType> const &left,
                                                          LinearVariableMapWrapper<ScalarType> const &right);
  friend LinearVariableMapWrapper<ScalarType> operator+<>(LinearVariableMapWrapper<ScalarType> const &left,
                                                          LinearVariableMapWrapper<ScalarType> const &right);

  friend LinearVariableMapWrapper<ScalarType> operator*<>(LinearVariableMapWrapper<ScalarType> const &left,
                                                          ArmadilloVectorWrapper<ScalarType> const &right);
  friend LinearVariableMapWrapper<ScalarType> operator*<>(ArmadilloVectorWrapper<ScalarType> const &left,
                                                          LinearVariableMapWrapper<ScalarType> const &right);

  friend LinearVariableMapWrapper<ScalarType> operator*<>(LinearVariableMapWrapper<ScalarType> const &left,
                                                          ScalarType const &right);
  friend LinearVariableMapWrapper<ScalarType> operator*<>(ScalarType const &left,
                                                          LinearVariableMapWrapper<ScalarType> const &right);
  friend LinearVariableMapWrapper<ScalarType> operator/<>(const LinearVariableMapWrapper<ScalarType> &left,
                                                          ScalarType const &right);

  friend LinearVariableMapWrapper<ScalarType> operator*<>(LinearVariableMapWrapper<ScalarType> const &left,
                                                          int const &right);
  friend LinearVariableMapWrapper<ScalarType> operator*<>(int const &left,
                                                          LinearVariableMapWrapper<ScalarType> const &right);
  friend LinearVariableMapWrapper<ScalarType> operator/<>(const LinearVariableMapWrapper<ScalarType> &left,
                                                          int const &right);

  friend LinearVariableMapWrapper<ScalarType> operator*<>(LinearVariableMapWrapper<ScalarType> const &left,
                                                          unsigned int const &right);
  friend LinearVariableMapWrapper<ScalarType> operator*<>(unsigned int const &left,
                                                          LinearVariableMapWrapper<ScalarType> const &right);
  friend LinearVariableMapWrapper<ScalarType> operator/<>(const LinearVariableMapWrapper<ScalarType> &left,
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
  typedef LinearVariableWrapper<ScalarType> LinearVariableType;

  /// Raw linear variable map type.
  typedef std::map<std::string, LinearVariableType> RawLinearVariableMapType;
  /// Iterator on a raw linear variable map type.
  typedef typename RawLinearVariableMapType::const_iterator const_iterator;
  typedef typename RawLinearVariableMapType::iterator iterator;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  LinearVariableMapWrapper() {}

  /// Constructor from a raw linear variable map.
  LinearVariableMapWrapper(const RawLinearVariableMapType &map) { m_RawLinearVariableMap = map; }

  /// Copy contructor.
  LinearVariableMapWrapper(const LinearVariableMapWrapper &other) {
    m_RawLinearVariableMap = other.m_RawLinearVariableMap;
  }

  /// Destructor.
  ~LinearVariableMapWrapper() {}


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Converts the linear variable map back to a raw linear variable map.
  RawLinearVariableMapType GetRawLinearVariableMap() const { return (m_RawLinearVariableMap); }

  /// Returns the length of the map.
  unsigned int size() const { return m_RawLinearVariableMap.size(); }
//    /// Resizes the map to a length \e n.
//    void resize(const unsigned int n) { m_RawLinearVariableMap.resize(n); }

  /// Fills the linear variable.
  void fill(ScalarType const &value) {
    fill_visitor<ScalarType> visitor(value);
    for (auto it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
  }

  /// Returns an iterator on the first element of the raw linear variable map.
  const_iterator begin() const noexcept { return m_RawLinearVariableMap.begin(); }
  iterator begin()       noexcept { return m_RawLinearVariableMap.begin(); }

  /// Returns an iterator on the end of the raw linear variable map.
  const_iterator end() const noexcept { return m_RawLinearVariableMap.end(); }
  iterator end()       noexcept { return m_RawLinearVariableMap.end(); }

  /// Returns an iterator on the searched element, or on the map end if the key \e str does not exist.
  const_iterator find(std::string const &str) const { return m_RawLinearVariableMap.find(str); }
  iterator find(std::string const &str) { return m_RawLinearVariableMap.find(str); }

  /// Counts the number of elements with a specific key \e str.
  unsigned int count(std::string const &str) const { return m_RawLinearVariableMap.count(str); }

  /// Erase a specific entry.
  void erase(const std::string &key) { m_RawLinearVariableMap.erase(key); }

  /// Clears the map, removing all elements.
  void clear() { m_RawLinearVariableMap.clear(); }

  /// Vectorizes the linear variable map.
  VectorType vectorize() const {
    VectorType out;
    for (auto it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it) {
      out.push_back(it->second.vectorize());
    }
    return out;
  }

  /// Vectorizes the linear variable map and returns the detected stucture and keys, for later reconstruction.
  VectorType vectorize(std::vector<std::vector<unsigned int>> &structure,
                       std::vector<std::string> &keys) const {
    VectorType out;
    structure.resize(0);
    keys.resize(0);
    for (auto it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it) {
      std::vector<unsigned int> aux;
      keys.push_back(it->first);
      out.push_back(it->second.vectorize(aux));
      structure.push_back(aux);
    }
    return out;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns sum of squares of all elements.
  ScalarType sum_of_squares() const {
    sum_of_squares_visitor<ScalarType> visitor;
    ScalarType result(0);
    for (const_iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      result += boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
    return result;
  }

  /// Returns the mean squares.
  VectorType mean_squares() const {
    sum_of_squares_visitor<ScalarType> ss_visitor;
    number_of_elements_visitor<ScalarType> ne_visitor;
    VectorType ms(m_RawLinearVariableMap.size());
    unsigned int k = 0;
    for (const_iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it, ++k) {
      ms[k] = boost::apply_visitor(ss_visitor)(it->second.GetRawLinearVariable())
          / boost::apply_visitor(ne_visitor)(it->second.GetRawLinearVariable());
    }
    return ms;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Access the linear variable corresponding to the key \e str for reading.
  LinearVariableType const &operator[](std::string const &str) const { return m_RawLinearVariableMap.at(str); }
  LinearVariableType const &at(std::string const &str) const { return m_RawLinearVariableMap.at(str); }
  /// Access the linear variable corresponding to the key \e str for reading or writing.
  LinearVariableType &operator[](std::string const &str) { return m_RawLinearVariableMap[str]; }

  /// Add rhs scalar to lhs linear variable map in situ.
  LinearVariableMapWrapper<ScalarType> &operator+=(ScalarType const &rhsScalar) {
    insitu_scalar_addition_visitor<ScalarType> visitor(rhsScalar);
    for (iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
    return *this;
  }
  /// Add rhs int to lhs linear variable map in situ.
  LinearVariableMapWrapper<ScalarType> &operator+=(int const &rhsInt) {
    return operator+=(ScalarType(rhsInt));
  }
  /// Add rhs unsigned int to lhs linear variable map in situ.
  LinearVariableMapWrapper<ScalarType> &operator+=(unsigned int const &rhsInt) {
    return operator+=(ScalarType(rhsInt));
  }

  /// Subtract rhs scalar from lhs linear variable map in situ.
  LinearVariableMapWrapper<ScalarType> &operator-=(ScalarType const &rhsScalar) {
    return operator+=(-rhsScalar);
  }
  /// Subtract rhs int from lhs linear variables map in situ.
  LinearVariableMapWrapper<ScalarType> &operator-=(int const &rhsInt) {
    return operator+=(-rhsInt);
  }
  /// Subtract rhs unsigned int from lhs linear variables map in situ.
  LinearVariableMapWrapper<ScalarType> &operator-=(unsigned int const &rhsInt) {
    return operator+=(-rhsInt);
  }

  /// Scalar multiplication in situ of lhs linear variable map by \e rhsScalar.
  LinearVariableMapWrapper<ScalarType> &operator*=(ScalarType const &rhsScalar) {
    insitu_scalar_multiplication_visitor<ScalarType> visitor(rhsScalar);
    for (iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
    return *this;
  }
  /// Int multiplication in situ of lhs linear variables map by \e rhsInt.
  LinearVariableMapWrapper<ScalarType> &operator*=(int const &rhsInt) {
    return operator*=(ScalarType(rhsInt));
  }
  /// Unsigned int multiplication in situ of lhs linear variables map by \e rhsInt.
  LinearVariableMapWrapper<ScalarType> &operator*=(unsigned int const &rhsInt) {
    return operator*=(ScalarType(rhsInt));
  }

  /// Scalar division of lhs linear variable map in situ by \e rhsScalar.
  LinearVariableMapWrapper<ScalarType> &operator/=(ScalarType const &rhsScalar) {
    return operator*=(1 / rhsScalar);
  }
  /// Int division of lhs linear variables map in situ by \e rhsInt.
  LinearVariableMapWrapper<ScalarType> &operator/=(int const &rhsInt) {
    return operator*=(1.0 / ScalarType(rhsInt));
  }
  /// Unsigned int division of lhs linear variables map in situ by \e rhsInt.
  LinearVariableMapWrapper<ScalarType> &operator/=(unsigned int const &rhsInt) {
    return operator*=(1.0 / ScalarType(rhsInt));
  }

  /// Add rhs linear variable map to lhs linear variable map in situ.
  LinearVariableMapWrapper<ScalarType> &operator+=(LinearVariableMapWrapper<ScalarType> const &right) {
    insitu_addition_visitor<ScalarType> visitor;
    for (iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      boost::apply_visitor(visitor, it->second.GetRawLinearVariable(),
                           right.m_RawLinearVariableMap.at(it->first).GetRawLinearVariable());
    return *this;
  }
  /// Subtract rhs linear variable map from lhs /LinearVariableMapWrapper.h:1 in situ.
  LinearVariableMapWrapper<ScalarType> &operator-=(LinearVariableMapWrapper<ScalarType> const &right) {
    insitu_subtraction_visitor<ScalarType> visitor;
    for (iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      boost::apply_visitor(visitor, it->second.GetRawLinearVariable(),
                           right.m_RawLinearVariableMap.at(it->first).GetRawLinearVariable());
    return *this;
  }
  /// Unary minus operator.
  LinearVariableMapWrapper<ScalarType> operator-() const {
    unary_minus_visitor<ScalarType> visitor;
    RawLinearVariableMapType rawList;
    for (const_iterator it = m_RawLinearVariableMap.begin(); it != m_RawLinearVariableMap.end(); ++it)
      rawList[it->first] = boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
    return LinearVariableMapWrapper<ScalarType>(rawList);
  }
  /// Equality testing operator.
  const bool operator==(const LinearVariableMapWrapper<ScalarType> &t) const {
    return m_RawLinearVariableMap.size() == t.m_RawLinearVariableMap.size() &&
        std::equal(m_RawLinearVariableMap.begin(), m_RawLinearVariableMap.end(), t.m_RawLinearVariableMap.begin());
  }

  virtual const bool compare(const LinearVariableMapWrapper<ScalarType> &t, const float tolerance) {
    if (m_RawLinearVariableMap.size() != t.m_RawLinearVariableMap.size()) return false;
    auto it1 = m_RawLinearVariableMap.begin();
    auto it2 = t.m_RawLinearVariableMap.begin();
    while (it1 != m_RawLinearVariableMap.end()) {
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
    ar & m_RawLinearVariableMap;
  }
  friend class boost::serialization::access;

  /// Raw linear variable map.
  mutable RawLinearVariableMapType m_RawLinearVariableMap;

};

/// Unvectorizes the input vector, converting it back to a linear variable map based on the input structure and keys.
template<class ScalarType>
LinearVariableMapWrapper<ScalarType> unvectorize(ArmadilloVectorWrapper<ScalarType> const &vec,
                                                 std::vector<std::vector<unsigned int>> const &structure,
                                                 std::vector<std::string> &keys) {
  LinearVariableMapWrapper<ScalarType> out;
  const unsigned long nbOfElements = keys.size();

  unsigned int currentIndex = 0;
  for (unsigned int k = 0; k < nbOfElements; ++k) {
    const unsigned long size = structure[k].size();

    if (size == 0) {
      out[keys[k]] = vec[currentIndex];
      currentIndex += 1;

    } else if (size == 1) {
      out[keys[k]] = vec.subvec(currentIndex, currentIndex + structure[k][0] - 1);
      currentIndex += structure[k][0];

    } else if (size == 2) {
      out[keys[k]] = vec.subvec(currentIndex, currentIndex + structure[k][0] * structure[k][1] - 1)
          .convert_to_matrix_row_wise(structure[k][0], structure[k][1]);
      currentIndex += structure[k][0] * structure[k][1];

    } else {
      const unsigned long nbOfMatrices = (size - 3) / 2;
      MatrixListWrapper<ScalarType> ml(nbOfMatrices);
      for (unsigned int i = 0; i < nbOfMatrices; ++i) {
        const unsigned int ni = structure[k][2 * i] * structure[k][2 * i + 1];
        ml[i] = vec.subvec(currentIndex, currentIndex + ni - 1)
            .convert_to_matrix_row_wise(structure[k][2 * i], structure[k][2 * i + 1]);
        currentIndex += ni;
      }
      out[keys[k]] = ml;
    }
  }
  return out;
}

//#include "LinearVariableMapWrapper_friend.h"

#endif /* _LinearVariableMapWrapper_h */
