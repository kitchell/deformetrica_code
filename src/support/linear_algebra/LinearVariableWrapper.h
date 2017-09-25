/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _LinearVariableWrapper_h
#define _LinearVariableWrapper_h


/// Core files.
#include "LinearVariableVariantVisitors.h"

/// Non-core files.
#include "ArmadilloVectorWrapper.h"
#include "ArmadilloMatrixWrapper.h"
#include "MatrixListWrapper.h"

/// Librairies files.
#include "boost/variant.hpp"
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>

#include "Tolerance.hpp"
using namespace def::algebra::utils;

template<class ScalarType>
class LinearVariableWrapper;

template<class ScalarType>
LinearVariableWrapper<ScalarType> operator+(const LinearVariableWrapper<ScalarType> &left,
                                            const LinearVariableWrapper<ScalarType> &right);
template<class ScalarType>
LinearVariableWrapper<ScalarType> operator-(const LinearVariableWrapper<ScalarType> &left,
                                            const LinearVariableWrapper<ScalarType> &right);

template<class ScalarType>
LinearVariableWrapper<ScalarType> operator*(LinearVariableWrapper<ScalarType> const &left, ScalarType const &right);
template<class ScalarType>
LinearVariableWrapper<ScalarType> operator*(ScalarType const &left, LinearVariableWrapper<ScalarType> const &right);

template<class ScalarType>
LinearVariableWrapper<ScalarType> operator/(const LinearVariableWrapper<ScalarType> &left, ScalarType const &right);

template<class ScalarType, typename VariantType>
VariantType recast(LinearVariableWrapper<ScalarType> const &var);
template<class ScalarType, typename VariantType>
VariantType recast(std::vector<LinearVariableWrapper<ScalarType>> const &vec);

/**
 *  \brief      LinearVariableWrapper object class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    // TODO
 *
 */
template<class ScalarType>
class LinearVariableWrapper : public Tolerance<LinearVariableWrapper<ScalarType>>  {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Friend methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  friend LinearVariableWrapper<ScalarType> operator-<>(LinearVariableWrapper<ScalarType> const &left,
                                                       LinearVariableWrapper<ScalarType> const &right);
  friend LinearVariableWrapper<ScalarType> operator+<>(LinearVariableWrapper<ScalarType> const &left,
                                                       LinearVariableWrapper<ScalarType> const &right);

  friend LinearVariableWrapper<ScalarType> operator*<>(LinearVariableWrapper<ScalarType> const &left,
                                                       ScalarType const &right);
  friend LinearVariableWrapper<ScalarType> operator*<>(ScalarType const &left,
                                                       LinearVariableWrapper<ScalarType> const &right);

  friend LinearVariableWrapper<ScalarType> operator/<>(const LinearVariableWrapper<ScalarType> &left,
                                                       ScalarType const &right);

  /// Converts the linear variable back to either a scalar, vector, matrix or matrix list.
  template<typename VariantType>
  friend VariantType recast(LinearVariableWrapper<ScalarType> const &var) {
    return boost::get<VariantType>(var.m_RawLinearVariable);
  }
  /// Converts the vector of linear variable back to a vector of either scalars, vectors, matrices or matrix lists.
  template<typename VariantType>
  friend std::vector<VariantType> recast(std::vector<LinearVariableWrapper<ScalarType>> const &vec) {
    std::vector<VariantType> result;
    for (unsigned int k = 0; k < vec.size(); ++k)
      result.push_back(boost::get<VariantType>(vec[k].m_RawLinearVariable));
    return result;
  }

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

  /// Raw linear variable type.
  typedef boost::variant<ScalarType, VectorType, MatrixType, MatrixListType> RawLinearVariableType;

  /// Raw iterator type on a linear variable.
  typedef boost::variant<std::vector<ScalarType *>, typename MatrixType::iterator,
                         typename MatrixListType::iterator> RawLinearVariableIteratorType;
  /// Raw const iterator type on a linear variable.
  typedef boost::variant<std::vector<ScalarType>, typename MatrixType::const_iterator,
                         typename MatrixListType::const_iterator> RawLinearVariableConstIteratorType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterators :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Iterator.
  class iterator : public std::iterator<
      std::bidirectional_iterator_tag,         // iterator_category
      ScalarType,                              // value_type
      unsigned int,                            // difference_type
      const ScalarType *,                       // pointer
      ScalarType                               // reference
  > {
   public:
    explicit iterator(RawLinearVariableIteratorType const &it) : m_RawLinearVariableIterator(it) {}

    iterator &operator++() {
      iterator_increment_visitor<ScalarType> visitor;
      boost::apply_visitor(visitor)(m_RawLinearVariableIterator);
      return *this;
    }

    bool operator==(iterator const &other) const {
      iterator_equality_visitor<ScalarType> visitor;
      return boost::apply_visitor(visitor, m_RawLinearVariableIterator, other.m_RawLinearVariableIterator);
    }
    bool operator!=(iterator const &other) const { return !(*this == other); }

    ScalarType &operator*() const {
      iterator_dereferencing_visitor<ScalarType> visitor;
      return boost::apply_visitor(visitor)(m_RawLinearVariableIterator);
    }

   private:
    RawLinearVariableIteratorType m_RawLinearVariableIterator;
  };

  /// Const iterator.
  class const_iterator : public std::iterator<
      std::bidirectional_iterator_tag,         // iterator_category
      const ScalarType,                        // value_type
      unsigned int,                            // difference_type
      const ScalarType *const,                 // pointer
      const ScalarType                         // reference
  > {
   public:
    explicit const_iterator(RawLinearVariableConstIteratorType const &it) : m_RawLinearVariableConstIterator(it) {}

    const_iterator &operator++() {
      const_iterator_increment_visitor<ScalarType> visitor;
      boost::apply_visitor(visitor)(m_RawLinearVariableConstIterator);
      return *this;
    }

    bool operator==(const_iterator const &other) const {
      const_iterator_equality_visitor<ScalarType> visitor;
      return boost::apply_visitor(visitor, m_RawLinearVariableConstIterator, other.m_RawLinearVariableConstIterator);
    }
    bool operator!=(const_iterator const &other) const { return !(*this == other); }

    ScalarType const &operator*() const {
      const_iterator_dereferencing_visitor<ScalarType> visitor;
      return boost::apply_visitor(visitor)(m_RawLinearVariableConstIterator);
    }

   private:
    RawLinearVariableConstIteratorType m_RawLinearVariableConstIterator;
  };


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  LinearVariableWrapper() {}

  /// Copy contructor.
  LinearVariableWrapper(const LinearVariableWrapper &other) { m_RawLinearVariable = other.m_RawLinearVariable; }

  /// Constructor from a raw linear variable.
  LinearVariableWrapper(const RawLinearVariableType &var) { m_RawLinearVariable = var; }

  /// Constructor from possible variant types.
  LinearVariableWrapper(const ScalarType &s) { m_RawLinearVariable = s; }
  LinearVariableWrapper(const VectorType &v) { m_RawLinearVariable = v; }
  LinearVariableWrapper(const MatrixType &m) { m_RawLinearVariable = m; }
  LinearVariableWrapper(const MatrixListType &ml) { m_RawLinearVariable = ml; }

  /// Destructor.
  ~LinearVariableWrapper() {}


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Converts the linear variable back to a raw linear variable.
  // Only for LinearVariableMapWrapper use, to be avoided in the rest of Deformetrica.
  RawLinearVariableType const &GetRawLinearVariable() const { return (m_RawLinearVariable); }
  RawLinearVariableType &GetRawLinearVariable() { return (m_RawLinearVariable); }

  /// Converts the linear variable back to either a scalar, vector, matrix or matrix list.
  template<typename VariantType>
  VariantType recast(VariantType &cast) const { return boost::get<VariantType>(m_RawLinearVariable); }

  /// Vectorizes the linear variable, row by row if a matrix is involved.
  VectorType vectorize() const {
    vectorize_visitor<ScalarType> visitor;
    return boost::apply_visitor(visitor)(m_RawLinearVariable);
  }
  /// Vectorizes the linear variable and returns the detected size parameters.
  VectorType vectorize(std::vector<unsigned int> &sizeParameters) const {
    vectorize_with_parameters_visitor<ScalarType> visitor;
    VectorType out = boost::apply_visitor(visitor)(m_RawLinearVariable);
    sizeParameters = visitor.GetParams();
    return out;
  }

  /// Fills the linear variable.
  void fill(ScalarType const &value) {
    fill_visitor<ScalarType> visitor(value);
    boost::apply_visitor(visitor)(m_RawLinearVariable);
  }

  /// Returns an iterator on the first scalar element of the linear variable.
  iterator begin() {
    iterator_begin_visitor<ScalarType> visitor;
    return iterator(boost::apply_visitor(visitor)(m_RawLinearVariable));
  }
  const_iterator begin() const {
    const_iterator_begin_visitor<ScalarType> visitor;
    return const_iterator(boost::apply_visitor(visitor)(m_RawLinearVariable));
  }

  /// Returns an iterator on the end of the linear variable.
  iterator end() {
    iterator_end_visitor<ScalarType> visitor;
    return iterator(boost::apply_visitor(visitor)(m_RawLinearVariable));
  }
  const_iterator end() const {
    const_iterator_end_visitor<ScalarType> visitor;
    return const_iterator(boost::apply_visitor(visitor)(m_RawLinearVariable));
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the number of elements.
  unsigned int n_elem() const {
    number_of_elements_visitor<ScalarType> visitor;
    return boost::apply_visitor(visitor)(m_RawLinearVariable);
  }

  /// Returns the sum of all elements.
  ScalarType sum() const {
    sum_visitor<ScalarType> visitor;
    return boost::apply_visitor(visitor)(m_RawLinearVariable);
  }

  /// Returns the sum of squares of all elements.
  ScalarType sum_of_squares() const {
    sum_of_squares_visitor<ScalarType> visitor;
    return boost::apply_visitor(visitor)(m_RawLinearVariable);
  }

  /// Returns the mean squares.
  ScalarType mean_squares() const {
    sum_of_squares_visitor<ScalarType> ss_visitor;
    number_of_elements_visitor<ScalarType> ne_visitor;
    return boost::apply_visitor(ss_visitor)(m_RawLinearVariable)
        / boost::apply_visitor(ne_visitor)(m_RawLinearVariable);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Add rhs scalar to lhs linear variable in situ.
  LinearVariableWrapper<ScalarType> &operator+=(ScalarType const &rhsScalar) {
    insitu_scalar_addition_visitor<ScalarType> visitor(rhsScalar);
    boost::apply_visitor(visitor)(m_RawLinearVariable);
    return *this;
  }

  /// Subtract rhs scalar from lhs linear variable in situ.
  LinearVariableWrapper<ScalarType> &operator-=(ScalarType const &rhsScalar) {
    return operator+=(-rhsScalar);
  }

  /// Scalar multiplication in situ of lhs linear variable by \e rhsScalar.
  LinearVariableWrapper<ScalarType> &operator*=(ScalarType const &rhsScalar) {
    insitu_scalar_multiplication_visitor<ScalarType> visitor(rhsScalar);
    boost::apply_visitor(visitor)(m_RawLinearVariable);
    return *this;
  }

  /// Scalar division of lhs linear variable in situ by \e rhsScalar.
  LinearVariableWrapper<ScalarType> &operator/=(ScalarType const &rhsScalar) {
    return operator*=(1 / rhsScalar);
  }

  /// Add rhs linear variable to lhs linear variable in situ.
  LinearVariableWrapper<ScalarType> &operator+=(LinearVariableWrapper<ScalarType> const &right) {
    insitu_addition_visitor<ScalarType> visitor;
    boost::apply_visitor(visitor, m_RawLinearVariable, right.m_RawLinearVariable);
    return *this;
  }

  /// Subtract rhs linear variable from lhs linear variable in situ.
  LinearVariableWrapper<ScalarType> &operator-=(LinearVariableWrapper<ScalarType> const &right) {
    insitu_subtraction_visitor<ScalarType> visitor;
    boost::apply_visitor(visitor, m_RawLinearVariable, right.m_RawLinearVariable);
    return *this;
  }

  /// Unary minus operator.
  LinearVariableWrapper<ScalarType> operator-() const {
    unary_minus_visitor<ScalarType> visitor;
    return LinearVariableWrapper<ScalarType>(boost::apply_visitor(visitor)(m_RawLinearVariable));
  }

  const bool operator==(const LinearVariableWrapper<ScalarType> &t) const {
    if (m_RawLinearVariable.empty() != t.m_RawLinearVariable.empty()) return false;
    if (m_RawLinearVariable.which() != t.m_RawLinearVariable.which()) return false;

    switch (m_RawLinearVariable.which()) {
      case 0: return (boost::get<ScalarType>(m_RawLinearVariable) == boost::get<ScalarType>(t.m_RawLinearVariable));
      case 1: return (boost::get<VectorType>(m_RawLinearVariable) == boost::get<VectorType>(t.m_RawLinearVariable));
      case 2: return (boost::get<MatrixType>(m_RawLinearVariable) == boost::get<MatrixType>(t.m_RawLinearVariable));
      case 3:
        return (boost::get<MatrixListType>(m_RawLinearVariable) == boost::get<MatrixListType>(t.m_RawLinearVariable));
    }

    return false;
  }

  virtual const bool compare(const LinearVariableWrapper<ScalarType>& t, const float tolerance) {
    if (m_RawLinearVariable.empty() != t.m_RawLinearVariable.empty()) return false;
    if (m_RawLinearVariable.which() != t.m_RawLinearVariable.which()) return false;

    switch (m_RawLinearVariable.which()) {
      case 0: return ScalarCompare(tolerance).operator()(boost::get<ScalarType>(m_RawLinearVariable),boost::get<ScalarType>(t.m_RawLinearVariable));
      case 1: return boost::get<VectorType>(m_RawLinearVariable).compare(boost::get<VectorType>(t.m_RawLinearVariable), tolerance);
      case 2: return boost::get<MatrixType>(m_RawLinearVariable).compare(boost::get<MatrixType>(t.m_RawLinearVariable), tolerance);
      case 3: return boost::get<MatrixListType>(m_RawLinearVariable).compare(boost::get<MatrixListType>(t.m_RawLinearVariable), tolerance);
    }

    return false;
  }

 private :

  template<class Archive>
  void save(Archive &ar, const unsigned int version) const {
    auto which = m_RawLinearVariable.which();
    ar & which;

    switch (which) {
      case 0:ar & boost::get<ScalarType>(m_RawLinearVariable);
        break;
      case 1:ar & boost::get<VectorType>(m_RawLinearVariable);
        break;
      case 2:ar & boost::get<MatrixType>(m_RawLinearVariable);
        break;
      case 3:ar & boost::get<MatrixListType>(m_RawLinearVariable);
        break;
    }
  }

  template<class Archive>
  void load(Archive &ar, const unsigned int version) {
    int which;
    ar & which;

    switch (which) {
      case 0: {
        ScalarType t;
        ar & t;
        m_RawLinearVariable = t;
        break;
      }
      case 1: {
        VectorType t;
        ar & t;
        m_RawLinearVariable = t;
        break;
      }
      case 2: {
        MatrixType t;
        ar & t;
        m_RawLinearVariable = t;
        break;
      }
      case 3: {
        MatrixListType t;
        ar & t;
        m_RawLinearVariable = t;
        break;
      }
    }
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  friend class boost::serialization::access;

  /// Raw linear variable.
  RawLinearVariableType m_RawLinearVariable;

};

template<class ScalarType>
inline ScalarType dot_product(LinearVariableWrapper<ScalarType> const &left,
                              LinearVariableWrapper<ScalarType> const &right) {
  dot_product_visitor<ScalarType> visitor;
  return boost::apply_visitor(visitor)(left, right);
}

/// Vectorizes the input linear variable.
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> vectorize(LinearVariableWrapper<ScalarType> const &var) {
  vectorize_visitor<ScalarType> visitor;
  return boost::apply_visitor(visitor)(var.GetRawLinearVariable());
}

/// Vectorizes the input linear variable, along with the detected size parameters.
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> vectorize(LinearVariableWrapper<ScalarType> const &var,
                                             std::vector<unsigned int> &sizeParameters) {
  vectorize_with_parameters_visitor<ScalarType> visitor;
  ArmadilloVectorWrapper<ScalarType> out = boost::apply_visitor(visitor)(var.GetRawLinearVariable());
  sizeParameters = visitor.GetParams();
  return out;
}

/// Unvectorizes the input vector, converting it back to a linear variable based on the input size parameters.
template<class ScalarType>
LinearVariableWrapper<ScalarType> unvectorize(ArmadilloVectorWrapper<ScalarType> const &vec,
                                              std::vector<unsigned int> const &sizeParams) {
  LinearVariableWrapper<ScalarType> out;
  const unsigned long size = sizeParams.size();

  if (size == 0) { out = vec[0]; }
  else if (size == 1) { out = vec; }
  else if (size == 2) { out = vec.convert_to_matrix_row_wise(sizeParams[0], sizeParams[1]); }
  else {
    const unsigned int nbOfMatrices = (size - 3) / 2;
    unsigned int currentIndex = 0;
    MatrixListWrapper<ScalarType> ml(nbOfMatrices);
    for (unsigned int k = 0; k < nbOfMatrices; ++k) {
      const unsigned int nk = sizeParams[2 * k] * sizeParams[2 * k + 1];
      ml[k] = vec.subvec(currentIndex, currentIndex + nk - 1)
          .convert_to_matrix_row_wise(sizeParams[2 * k], sizeParams[2 * k + 1]);
      currentIndex += nk;
    }
    out = ml;
  }
  return out;
}

//#include "LinearVariableWrapper_friend.h"

#endif /* _LinearVariableWrapper_h */