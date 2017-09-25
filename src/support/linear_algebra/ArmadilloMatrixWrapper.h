/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _ArmadilloMatrixWrapper_h
#define _ArmadilloMatrixWrapper_h

/// Support files.
#include "ArmadilloVectorWrapper.h"

/// Libraries files.
#include "assert.h"
#include <iostream>
#include <armadillo>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/binary_object.hpp>

#include "Tolerance.hpp"

#ifndef DEFORMETRICA_CONFIG
#include "DeformetricaConfig.h"
#endif

#ifdef USE_DOUBLE_PRECISION
#define ScalarBaseType double
#define ScalarPrecisionType float
#else
#define ScalarBaseType float
#define ScalarPrecisionType double
#endif

template<class ScalarType>
class ArmadilloMatrixWrapper;

template<class ScalarType>
class ArmadilloVectorWrapper;

using namespace def::algebra::utils;

/**
 *  \brief      Mathematical matrix class (Armadillo).
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The ArmadilloMatrixWrapper class is a wrapping for mathematical matrix class, templated by type of element
 *              and using Armadillo library.
 *  \warning    Elements are stored with column-major ordering (ie. column by column).
 *
 */
template<class ScalarType>
class ArmadilloMatrixWrapper : public Tolerance<ArmadilloMatrixWrapper<ScalarType>> {

 public :

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Vector type
  typedef typename arma::Col<ScalarType> ArmadilloVectorType;
  /// Matrix type.
  typedef typename arma::Mat<ScalarType> ArmadilloMatrixType;

  /// Iterator on a raw Armadillo matrix type.
  typedef typename ArmadilloMatrixType::iterator iterator;
  typedef typename ArmadilloMatrixType::const_iterator const_iterator;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  ArmadilloMatrixWrapper() : m_Matrix() {}

  /// Matrix constructor of size \e r rows by \e c columns.
  explicit ArmadilloMatrixWrapper(unsigned const &r, unsigned const &c) : m_Matrix(r, c) {
    m_Matrix.fill(1e5);
  }

  /// Construct a matrix of size \e r rows by \e c columns, and all elements equal to \e v0.
  explicit ArmadilloMatrixWrapper(unsigned const &r, unsigned const &c, ScalarType const &v0) : m_Matrix(r, c) {
    m_Matrix.fill(v0);
  }

//    /// Construct a diagonal squared matrix of size \e n rows by \e n columns, and diagonal elements equal to \e v0.
//    explicit ArmadilloMatrixWrapper(unsigned const& n, ScalarType const& v) { m_Matrix.eye(n, n); m_Matrix *= v; }

  /// Construct a matrix of size \e r rows by \e c columns, initialized by a (column-wise) memory block.
  /// \warning The auxiliary memory \e data_block is not copied !
  explicit ArmadilloMatrixWrapper(ScalarType *data_block,
                                  unsigned const &r,
                                  unsigned const &c) : m_Matrix(data_block, r, c, false, true) {}

  /// Constructor which converts a (column) vector to a matrix.
  ArmadilloMatrixWrapper(const ArmadilloVectorWrapper<ScalarType> &v) : m_Matrix(v.toArmadillo()) {}

  /// Constructor from a std::vector<VectorType>.
  ArmadilloMatrixWrapper(const std::vector<ArmadilloVectorWrapper<ScalarType>> &v) {
    m_Matrix.set_size(v.size(), v[0].size());
    for (unsigned int k = 0; k < v.size(); ++k) { m_Matrix.row(k) = v[k].toArmadillo().t(); }
  }

  /// Constructor from a std::vector<ScalarType>.
  ArmadilloMatrixWrapper(const std::vector<ScalarType> &v) {
    m_Matrix.set_size(v.size(), 1);
    for (unsigned int k = 0; k < v.size(); ++k) { m_Matrix(k, 0) = v[k]; }
  }

  /// Copy constructor.
  ArmadilloMatrixWrapper(const ArmadilloMatrixWrapper<ScalarType> &other) : m_Matrix(other.m_Matrix) {}

  /// Special constructor (do not use it in Deformetrica!).
  ArmadilloMatrixWrapper(const ArmadilloMatrixType &M) : m_Matrix(M) {}

  /// Access directly to the Armadillo matrix (do not use it in Deformetrica!).
  ArmadilloMatrixType const &toArmadillo() const { return m_Matrix; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::Mat<ScalarType> &get_aramadillo_mat() /*const*/ { return (m_Matrix); }
  arma::Mat<ScalarType> get_aramadillo_mat() const { return (m_Matrix); }

  /// Returns the number of rows.
  unsigned rows() const { return m_Matrix.n_rows; }

  /// Returns the number of columns.
  unsigned columns() const { return m_Matrix.n_cols; }

  /// Returns the number of columns.
  unsigned cols() const { return m_Matrix.n_cols; }

  /// Returns the number of elements (This equals to rows() * cols()).
  unsigned size() const { return m_Matrix.n_elem; }

  /// Returns the number of elements (This equals to rows() * cols()).
  unsigned n_elem() const { return m_Matrix.n_elem; }
/*
  ScalarType       * memory_pointer_col_wise() const {  ArmadilloVectorType result = arma::vectorise(m_Matrix, 0); return result.memptr(); }

	/// Access the contiguous block storing the elements in the matrix row-wise.
	ScalarType      * memory_pointer_row_wise() {  ArmadilloVectorRowType result = arma::vectorise(m_Matrix, 1); return result.memptr(); }

	/// Access the contiguous block storing the elements in the matrix row-wise.
	ScalarType const* memory_pointer_row_wise() const { ArmadilloVectorRowType result = arma::vectorise(m_Matrix, 1); return result.memptr(); }
*/
  /// Resize to r rows by c columns. Old data lost.
  bool set_size(unsigned r, unsigned c) {
    m_Matrix.set_size(r, c);
    return true;
  }

  /// Set all elements of matrix to specified value.
  ArmadilloMatrixWrapper<ScalarType> &fill(ScalarType const &scalar) {
    m_Matrix.fill(scalar);
    return *this;
  }

  /// Reshapes the matrix as a diagonal one, diagonal elements being equal to \e value.
  void diagonal_matrix(unsigned const &n, ScalarType const &value) {
    m_Matrix.eye(n, n);
    m_Matrix *= value;
  }

  /// Gets \e n rows beginning at \e rowstart.
  ArmadilloMatrixWrapper<ScalarType> get_n_rows(unsigned rowstart, unsigned n) const {
    return ArmadilloMatrixWrapper<ScalarType>(m_Matrix.rows(rowstart, rowstart + n - 1));
  }

  /// Gets \e n columns beginning at \e colstart.
  ArmadilloMatrixWrapper<ScalarType> get_n_columns(unsigned colstart, unsigned n) const {
    return ArmadilloMatrixWrapper<ScalarType>(m_Matrix.cols(colstart, colstart + n - 1));
  }

  /// Sets columns to those in \e M, starting at \e starting_column, then return *this.
  ArmadilloMatrixWrapper<ScalarType> &set_columns(unsigned starting_column,
                                                  ArmadilloMatrixWrapper<ScalarType> const &M) {
    m_Matrix(arma::span(0, 0 + M.rows() - 1), arma::span(starting_column, starting_column + M.cols() - 1)) = M.m_Matrix;
    return *this;
  }

  /// Sets values of this matrix to those of \e M, starting at (\e top, \e left).
  ArmadilloMatrixWrapper<ScalarType> &update(ArmadilloMatrixWrapper<ScalarType> const &M,
                                             unsigned top = 0,
                                             unsigned left = 0) {
    m_Matrix(arma::span(top, top + M.rows() - 1), arma::span(left, left + M.cols() - 1)) = M.m_Matrix;
    return *this;
  }

  /// Gets a vector equal to the given row.
  ArmadilloVectorWrapper<ScalarType> get_row(unsigned r) const {
    return ArmadilloVectorWrapper<ScalarType>(m_Matrix.row(r).t());
  }

  /// Gets a vector equal to the given column.
  ArmadilloVectorWrapper<ScalarType> get_column(unsigned c) const {
    return ArmadilloVectorWrapper<ScalarType>(m_Matrix.col(c));
  }

  /// Sets the \e i-th row to \e v.
  void set_row(unsigned i, ArmadilloVectorWrapper<ScalarType> const &v) { m_Matrix.row(i) = v.toArmadillo().t(); }

  /// Sets the elements of the i'th row to \e value, then return *this.
  ArmadilloMatrixWrapper<ScalarType> &set_row(unsigned i, ScalarType value) {
    m_Matrix.row(i).fill(value);
    return *this;
  }

  /// Adds \e v to the \e i-th row.
  void increment_row(unsigned i, ArmadilloVectorWrapper<ScalarType> const &v) {
    m_Matrix.row(i) += v.toArmadillo().t();
  }

  /// Multiplies each column d by \e v[d].
  void multiply_cols(ArmadilloVectorWrapper<ScalarType> const &v) {
    for (unsigned int d = 0; d < m_Matrix.n_cols; ++d) { m_Matrix.col(d) *= v[d]; }
  }

  /// Sets \e j-th column to \e v.
  void set_column(unsigned j, ArmadilloVectorWrapper<ScalarType> const &v) { m_Matrix.col(j) = v.toArmadillo(); }

  /// Sets this matrix to an identity matrix.
  void set_identity() { m_Matrix.eye(); }

  /// Converts a matrix to a vector with a column-wise order.
  ArmadilloVectorWrapper<ScalarType> vectorise_col_wise() const {
    return ArmadilloVectorWrapper<ScalarType>(arma::vectorise(m_Matrix, 0));
  }
  /// Converts a matrix to a vector with a row-wise order.
  ArmadilloVectorWrapper<ScalarType> vectorise_row_wise() const {
    return ArmadilloVectorWrapper<ScalarType>(arma::trans(arma::vectorise(m_Matrix, 1)));
  }
  /// Converts a matrix to a vector with a row-wise order (alias).
  ArmadilloVectorWrapper<ScalarType> vectorize() const {
    return ArmadilloVectorWrapper<ScalarType>(arma::trans(arma::vectorise(m_Matrix, 1)));
  }

  const ScalarType *memptr() const { return m_Matrix.memptr(); }

  /// Returns an iterator on the first element of the raw Armadillo matrix.
  iterator begin() { return m_Matrix.begin(); }
  const_iterator begin() const { return m_Matrix.begin(); }

  /// Returns an iterator on the end of the raw Armadillo matrix.
  iterator end() { return m_Matrix.end(); }
  const_iterator end() const { return m_Matrix.end(); }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns Frobenius norm of matrix (sqrt of sum of squares of its elements).
  ScalarType frobenius_norm() const { return arma::norm(m_Matrix, "fro"); }

  //Returns the determinant of a square matrix (error if not square)
  ScalarType determinant() const { return arma::det(m_Matrix); }

  /// Returns transpose.
  ArmadilloMatrixWrapper<ScalarType> transpose() const { return ArmadilloMatrixWrapper<ScalarType>(m_Matrix.t()); }

  /// Returns sum of elements.
  ScalarType sum() const { return arma::accu(m_Matrix); }
  /// Returns sum of squares of all elements.
  ScalarType sum_of_squares() const { return arma::accu(arma::square(m_Matrix)); }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Access an element for reading.
  ScalarType const &operator()(unsigned r, unsigned c) const { return m_Matrix(r, c); }
  /// Access an element for reading or writing.
  ScalarType &operator()(unsigned r, unsigned c) { return m_Matrix(r, c); }

  /// Scalar division of lhs matrix in situ by \e rhsScalar.
  ArmadilloMatrixWrapper<ScalarType> &operator/=(ScalarType rhsScalar) {
    m_Matrix /= rhsScalar;
    return *this;
  }
  /// Scalar multiplication in situ of lhs matrix by \e rhsScalar.
  ArmadilloMatrixWrapper<ScalarType> &operator*=(ScalarType rhsScalar) {
    m_Matrix *= rhsScalar;
    return *this;
  }
  /// Add rhs to lhs matrix in situ.
  ArmadilloMatrixWrapper<ScalarType> &operator+=(ScalarType rhsScalar) {
    m_Matrix += rhsScalar;
    return *this;
  }
  /// Subtract rhs from lhs matrix in situ.
  ArmadilloMatrixWrapper<ScalarType> &operator-=(ScalarType rhsScalar) {
    m_Matrix -= rhsScalar;
    return *this;
  }

  /// Multiply lhs matrix in situ by rhs.
  ArmadilloMatrixWrapper<ScalarType> &operator*=(ArmadilloMatrixWrapper<ScalarType> const &rhsMatrix) {
    m_Matrix *= rhsMatrix.m_Matrix;
    return *this;
  }
  /// Add rhs to lhs matrix in situ.
  ArmadilloMatrixWrapper<ScalarType> &operator+=(ArmadilloMatrixWrapper<ScalarType> const &rhsMatrix) {
    m_Matrix += rhsMatrix.m_Matrix;
    return *this;
  }
  /// Subtract rhs from lhs matrix in situ.
  ArmadilloMatrixWrapper<ScalarType> &operator-=(ArmadilloMatrixWrapper<ScalarType> const &rhsMatrix) {
    m_Matrix -= rhsMatrix.m_Matrix;
    return *this;
  }

  /// Unary minus operator.
  ArmadilloMatrixWrapper<ScalarType> operator-() const { return ArmadilloMatrixWrapper<ScalarType>(-m_Matrix); }
  /// Unary plus operator.
  ArmadilloMatrixWrapper<ScalarType> operator+() const { return ArmadilloMatrixWrapper<ScalarType>(+m_Matrix); }

  const bool operator==(const ArmadilloMatrixWrapper<ScalarType> &t) const {
    if (m_Matrix.n_rows != t.m_Matrix.n_rows) return false;
    if (m_Matrix.n_cols != t.m_Matrix.n_cols) return false;

    const ScalarType *mem1 = m_Matrix.memptr();
    const ScalarType *mem2 = t.m_Matrix.memptr();
    return std::equal(mem1, mem1 + m_Matrix.n_rows * m_Matrix.n_cols, mem2);
  }

  virtual const bool compare(const ArmadilloMatrixWrapper<ScalarType> &Matrix, const float tolerance) {
    if (m_Matrix.n_rows != Matrix.m_Matrix.n_rows) return false;
    if (m_Matrix.n_cols != Matrix.m_Matrix.n_cols) return false;

    return (m_Matrix.size() == Matrix.m_Matrix.size())
        && std::equal(m_Matrix.begin(), m_Matrix.end(), Matrix.m_Matrix.begin(), ScalarCompare(tolerance));
  }

 private:

  template<class Archive>
  void save(Archive &ar, const unsigned int version) const {
    ScalarType *mem = (ScalarType *) m_Matrix.memptr();

    ar & m_Matrix.n_rows;
    ar & m_Matrix.n_cols;
    ar & boost::serialization::make_binary_object(mem, m_Matrix.n_rows * m_Matrix.n_cols * sizeof(ScalarType));
  }

  template<class Archive>
  void load(Archive &ar, const unsigned int version) {
    unsigned int rows;
    unsigned int cols;
    ar & rows;
    ar & cols;

    ScalarType *mem = new ScalarType[rows * cols];
    ar & boost::serialization::make_binary_object(mem, rows * cols * sizeof(ScalarType));

    m_Matrix = arma::mat(mem, rows, cols);
    delete[] mem;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  friend class boost::serialization::access;

 private :

  ArmadilloMatrixType m_Matrix;

};

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator+(const ArmadilloMatrixWrapper<ScalarType> &left,
                                             const ArmadilloMatrixWrapper<ScalarType> &right);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator-(const ArmadilloMatrixWrapper<ScalarType> &left,
                                             const ArmadilloMatrixWrapper<ScalarType> &right);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> &left,
                                             const ArmadilloMatrixWrapper<ScalarType> &right);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarType const &rightScalar);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarPrecisionType const &rightScalar);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(ScalarType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(ScalarPrecisionType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix);
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(ArmadilloMatrixWrapper<ScalarType> const &leftMatrix,
                                             ArmadilloVectorWrapper<ScalarType> const &rightVector);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarType const &rightScalar);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarPrecisionType const &rightScalar);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(ScalarType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(ScalarPrecisionType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(unsigned N, ScalarType const &value);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(unsigned N, ScalarPrecisionType const &value);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(ArmadilloVectorWrapper<ScalarType> const &values);
template<class ScalarType>
ScalarType trace(ArmadilloMatrixWrapper<ScalarType> const &M);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> chol(ArmadilloMatrixWrapper<ScalarType> const &M);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> inverse(ArmadilloMatrixWrapper<ScalarType> const &M);

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> inverse_sympd(ArmadilloMatrixWrapper<ScalarType> const &M);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> eigenvalues_sym(ArmadilloMatrixWrapper<ScalarType> const &M);

/// Computes the eigenvalues and eigenvectors.
template<class ScalarType>
void eigvec_sym(ArmadilloVectorWrapper<ScalarType> &eigVal,
                ArmadilloMatrixWrapper<ScalarType> &eigVec,
                ArmadilloMatrixWrapper<ScalarType> const &symMat) {
  arma::eig_sym(eigVal.toArmadillo(), eigVec.get_aramadillo_mat(), symMat.get_aramadillo_mat());
}

template<class ScalarType>
std::ostream &operator<<(std::ostream &os, ArmadilloMatrixWrapper<ScalarType> const &rhs);

template<class ScalarType>
ScalarType dot_product(ArmadilloMatrixWrapper<ScalarType> const &left, ArmadilloMatrixWrapper<ScalarType> const &right);

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> solve(ArmadilloMatrixWrapper<ScalarType> const &A,
                                         ArmadilloVectorWrapper<ScalarType> const &b);
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> solve(ArmadilloMatrixWrapper<ScalarType> const &A,
                                         ArmadilloMatrixWrapper<ScalarType> const &B);

/// Determinant of the input matrix \e M.
template<class ScalarType>
ScalarType det(ArmadilloMatrixWrapper<ScalarType> const &M);

template<class ScalarType>
ScalarType log_det(ArmadilloMatrixWrapper<ScalarType> const &M);

#undef ScalarBaseType
#undef ScalarPrecisionType

#endif /* _ArmadilloMatrixWrapper_h */
