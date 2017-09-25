/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _ArmadilloVectorWrapper_h
#define _ArmadilloVectorWrapper_h

#ifndef DEFORMETRICA_CONFIG
#include "DeformetricaConfig.h"
#endif

/// Support files.
#include "ArmadilloMatrixWrapper.h"

/// Libraries files.
#include <vector>
#include <iostream>
#include <armadillo>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/binary_object.hpp>

#include "Tolerance.hpp"

#ifdef USE_DOUBLE_PRECISION
#define ScalarBaseType double
#define ScalarPrecisionType float
#else
#define ScalarBaseType float
#define ScalarPrecisionType double
#endif

template<typename ScalarType>
class ArmadilloMatrixWrapper;

using namespace def::algebra::utils;

/**
 *  \brief      Mathematical vector class (Armadillo).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The ArmadilloVectorWrapper class is a wrapping for mathematical vector class, templated by type of element
 *              and using the Armadillo library.
 */
template<typename ScalarType>
class ArmadilloVectorWrapper : public Tolerance<ArmadilloVectorWrapper<ScalarType>> {

 public :

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Vector type
  typedef typename arma::Col<ScalarType> ArmadilloVectorType;
  /// Matrix type.
  typedef typename arma::Mat<ScalarType> ArmadilloMatrixType;

  /// Iterator on a raw Armadillo vector type.
  typedef typename ArmadilloVectorType::iterator iterator;
  typedef typename ArmadilloVectorType::const_iterator const_iterator;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  ArmadilloVectorWrapper() : m_Vector() {}

  /// Creates vector of \e len elements.
  explicit ArmadilloVectorWrapper(unsigned len) : m_Vector(len) { m_Vector.fill(1e5); }

  /// Creates vector of \e len elements, all set to \e v0.
  explicit ArmadilloVectorWrapper(unsigned len, ScalarType const &v0) : m_Vector(len) { m_Vector.fill(v0); }

  /// Constructor from a std::vector<ScalarType>.
  ArmadilloVectorWrapper(const std::vector<ScalarType> vec) {
    const unsigned int size = vec.size();
    m_Vector.set_size(size);
    for (unsigned int k = 0; k < size; ++k) { m_Vector(k) = vec[k]; }
  }

  /// Constructor from a raw Armadillo vector (do not use it in Deformetrica!).
  ArmadilloVectorWrapper(const ArmadilloVectorType &v) { m_Vector = v; }

  /// Copy constructor.
  ArmadilloVectorWrapper(const ArmadilloVectorWrapper<ScalarType> &other) { m_Vector = other.m_Vector; }

  /// Access directly to the Armadillo vector (do not use it in Deformetrica!).
  ArmadilloVectorType const &toArmadillo() const { return m_Vector; }
  ArmadilloVectorType &toArmadillo() { return m_Vector; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the length, number of elements, dimension of this vector.
  unsigned size() const { return m_Vector.n_rows; }
  /// Resizes to \e n elements.
  bool set_size(unsigned n) {
    m_Vector.set_size(n);
    return true;
  }

  /// Returns the number of elements.
  unsigned n_elem() const { return m_Vector.n_rows; }

  /// Gets value at element \e i.
  ScalarType get(unsigned int i) const { return m_Vector.at(i); }
  /// Gets contiguous values between \e first_index and \e last_index.
  ArmadilloVectorWrapper<ScalarType> subvec(const unsigned int first_index, const unsigned int last_index) const {
    return ArmadilloVectorWrapper<ScalarType>(m_Vector.subvec(first_index, last_index));
  }

  /// Sets all elements of vector to specified value.
  ArmadilloVectorWrapper<ScalarType> &fill(ScalarType const &scalar) {
    m_Vector.fill(scalar);
    return *this;
  }

  /// Sets values of this vector to those of \e V, starting at \e top.
  ArmadilloVectorWrapper<ScalarType> &update(ArmadilloVectorWrapper<ScalarType> const &V, unsigned top = 0) {
    m_Vector(arma::span(top, top + V.size() - 1)) = V.m_Vector;
    return *this;
  }

  /// Converts a (column) vector to a matrix of size \e rows x \e cols, in a column-wise fashion.
  ArmadilloMatrixWrapper<ScalarType> convert_to_matrix_col_wise(unsigned rows, unsigned cols) const {
    ArmadilloMatrixType result(m_Vector);
    result.reshape(rows, cols);
    return ArmadilloMatrixWrapper<ScalarType>(result);
  }
  /// Converts a (column) vector to a matrix of size \e rows x \e cols, in row-wise fashion.
  ArmadilloMatrixWrapper<ScalarType> convert_to_matrix_row_wise(unsigned rows, unsigned cols) const {
    ArmadilloMatrixType result(m_Vector);
    result.reshape(cols, rows);
    return ArmadilloMatrixWrapper<ScalarType>(arma::trans(result));
  }
  /// Same method than convert_to_matrix_row_wise, with another name.
  ArmadilloMatrixWrapper<ScalarType> unvectorize(unsigned rows, unsigned cols) const {
    ArmadilloMatrixType result(m_Vector);
    result.reshape(cols, rows);
    return ArmadilloMatrixWrapper<ScalarType>(arma::trans(result));
  }

  /// Adds a value at the end of the vector.
  void push_back(const ScalarType &s) {
    m_Vector.set_size(m_Vector.n_rows + 1);
    m_Vector(m_Vector.n_rows - 1) = s;
  }
  /// Concatenates vertically.
  void push_back(const ArmadilloVectorWrapper<ScalarType> &v) {
    m_Vector = arma::join_cols(m_Vector, v.m_Vector);
  }

  /// Returns an iterator on the first element of the raw Armadillo vector.
  iterator begin() { return m_Vector.begin(); }
  const_iterator begin() const { return m_Vector.begin(); }

  /// Returns an iterator on the end of the raw Armadillo vector.
  iterator end() { return m_Vector.end(); }
  const_iterator end() const { return m_Vector.end(); }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns magnitude (norm) of vector.
  ScalarType magnitude() const { return arma::norm(m_Vector, 2); } //INSTABLE WITH SCALAR-TYPE=FLOAT
  //ScalarType magnitude() const { return std::sqrt(arma::sum(arma::square(m_Vector))); }

  /// Returns sum of squares of elements.
  ScalarType squared_magnitude() const {
    ScalarType result = arma::norm(m_Vector, 2);
    return result * result;
  } //INSTABLE WITH SCALAR-TYPE=FLOAT
  //ScalarType squared_magnitude() const { return arma::sum(arma::square(m_Vector)); }

  /// Returns sum of elements.
  ScalarType sum() const { return arma::sum(m_Vector); }
  /// Returns sum of squares of elements.
  ScalarType sum_of_squares() const { return arma::sum(arma::square(m_Vector)); }

  /// Returns the mean of the vector.
  ScalarType mean() const { return arma::mean(m_Vector); }

  /// Returns a copy where all elements have been sqrt-ed.
  ArmadilloVectorWrapper<ScalarType> sqrt() const {
    return ArmadilloVectorWrapper<ScalarType>(arma::sqrt(m_Vector));
  }

  /// Raise every element to the square.
  void square() { arma::square(m_Vector); }

  /// Smallest value.
  ScalarType min_value() const { return m_Vector.min(); }
  /// Largest value.
  ScalarType max_value() const { return m_Vector.max(); }

  /// Smallest value index.
//	unsigned int index_min() const { return m_Vector.index_min(); } // REQUIRES ARMADILLO 7.400
  unsigned int index_min() const {
    unsigned int index = 0;
    for (unsigned int k = 1; k < m_Vector.size(); ++k) {
      if (m_Vector(k) < m_Vector(index)) { index = k; }
    }
    return index;
  }
  /// Largest value index.
//	unsigned int index_max() const { return m_Vector.index_max(); } // REQUIRES ARMADILLO 7.400
  unsigned int index_max() const {
    unsigned int index = 0;
    for (unsigned int k = 1; k < m_Vector.size(); ++k) {
      if (m_Vector(k) > m_Vector(index)) { index = k; }
    }
    return index;
  }

  /// Access the contiguous block storing the elements in the vector.
  ScalarType const *memptr() const { return m_Vector.memptr(); }

  /// Access the contiguous block storing the elements in the vector.
  ScalarType *memptr() { return m_Vector.memptr(); }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns reference to the element at index \e i for reading.
  ScalarType const &operator()(unsigned int i) const { return m_Vector(i); }
  /// Returns reference to the element at index \e i for reading or writing.
  ScalarType &operator()(unsigned int i) { return m_Vector(i); }

  /// Returns reference to the element at index \e i.
  /// \warning Please use operator() in Deformetrica instead of operator[].
  ScalarType &operator[](unsigned int i) { return m_Vector(i); }
  /// Returns reference to the element at index \e i.
  /// \warning Please use operator() in Deformetrica instead of operator[].
  ScalarType const &operator[](unsigned int i) const { return m_Vector(i); }

  /// Adds scalar \e rhs to lhs vector in situ.
  ArmadilloVectorWrapper<ScalarType> &operator+=(ScalarType const &rhs) {
    m_Vector += rhs;
    return *this;
  }
  /// Substracts scalar \e rhs to lhs vector in situ.
  ArmadilloVectorWrapper<ScalarType> &operator-=(ScalarType const &rhs) {
    m_Vector -= rhs;
    return *this;
  }
  /// Scalar multiplication of lhs vector in situ by \e rhsScalar.
  ArmadilloVectorWrapper<ScalarType> &operator*=(ScalarType const &rhs) {
    m_Vector *= rhs;
    return *this;
  }
  /// Scalar division of lhs vector in situ by \e rhsScalar.
  ArmadilloVectorWrapper<ScalarType> &operator/=(ScalarType const &rhs) {
    m_Vector /= rhs;
    return *this;
  }

  /// Adds \e rhs to lhs vector in situ.
  ArmadilloVectorWrapper<ScalarType> &operator+=(ArmadilloVectorWrapper<ScalarType> const &rhs) {
    m_Vector += rhs.m_Vector;
    return *this;
  }
  /// Substracts \e rhs to lhs vector in situ.
  ArmadilloVectorWrapper<ScalarType> &operator-=(ArmadilloVectorWrapper<ScalarType> const &rhs) {
    m_Vector -= rhs.m_Vector;
    return *this;
  }

  /// Unary minus operator.
  ArmadilloVectorWrapper<ScalarType> operator-() const { return ArmadilloVectorWrapper<ScalarType>(-m_Vector); }
  /// Unary plus operator.
  ArmadilloVectorWrapper<ScalarType> operator+() const { return ArmadilloVectorWrapper<ScalarType>(+m_Vector); }

  /// Scalar multiplication of lhs vector by \e scalar.
  ArmadilloVectorWrapper<ScalarType> operator*(ScalarType scalar) const {
    return ArmadilloVectorWrapper<ScalarType>(m_Vector * scalar);
  }
  /// Scalar division of lhs vector by \e scalar.
  ArmadilloVectorWrapper<ScalarType> operator/(ScalarType scalar) const {
    return ArmadilloVectorWrapper<ScalarType>(m_Vector / scalar);
  }

  const bool operator==(const ArmadilloVectorWrapper<ScalarType> &t) const {
    if (m_Vector.n_rows != t.m_Vector.n_rows) return false;

    const ScalarType *mem1 = m_Vector.memptr();
    const ScalarType *mem2 = t.m_Vector.memptr();
    return std::equal(mem1, mem1 + m_Vector.n_rows, mem2);
  }

  virtual const bool compare(const ArmadilloVectorWrapper<ScalarType> &Vector, const float tolerance) {
    return (m_Vector.size() == Vector.m_Vector.size())
        && std::equal(m_Vector.begin(), m_Vector.end(), Vector.m_Vector.begin(), ScalarCompare(tolerance));
  }

 private:

  template<class Archive>
  void save(Archive &ar, const unsigned int version) const {
    ScalarType *mem = (ScalarType *) m_Vector.memptr();

    ar & m_Vector.n_rows;
    ar & boost::serialization::make_binary_object(mem, m_Vector.n_rows * sizeof(ScalarType));
  }

  template<class Archive>
  void load(Archive &ar, const unsigned int version) {
    unsigned int rows;
    ar & rows;

    ScalarType *mem = new ScalarType[rows];
    ar & boost::serialization::make_binary_object(mem, rows * sizeof(ScalarType));

    m_Vector = arma::vec(mem, rows);
    delete[] mem;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  friend class boost::serialization::access;

 private :

  ArmadilloVectorType m_Vector;

};

template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator%(const ArmadilloVectorWrapper<ScalarType> &left,
                                             const ArmadilloVectorWrapper<ScalarType> &right);

template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> &left,
                                             const ArmadilloVectorWrapper<ScalarType> &right);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarType &rightScalar);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarPrecisionType &rightScalar);

template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> &left,
                                             const ArmadilloVectorWrapper<ScalarType> &right);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarType &rightScalar);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarPrecisionType &rightScalar);

template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarType &rightScalar);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarPrecisionType &rightScalar);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ScalarType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ScalarPrecisionType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector);

template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(const ScalarType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector);
template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(const ScalarPrecisionType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector);

template<typename ScalarType>
std::ostream &operator<<(std::ostream &os, ArmadilloVectorWrapper<ScalarType> const &rhs);

template<typename ScalarType>
ScalarType dot_product(ArmadilloVectorWrapper<ScalarType> const &v1, ArmadilloVectorWrapper<ScalarType> const &v2);

template<typename ScalarType>
ArmadilloVectorWrapper<ScalarType> cross_3d(ArmadilloVectorWrapper<ScalarType> const &v1,
                                            ArmadilloVectorWrapper<ScalarType> const &v2);

#undef ScalarBaseType
#undef ScalarPrecisionType

#endif /* _ArmadilloVectorWrapper_h */
