/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#ifndef _MatrixListWrapper_h
#define _MatrixListWrapper_h

/// Support files.
#include "ArmadilloVectorWrapper.h"
#include "ArmadilloMatrixWrapper.h"

/// Librairies files.
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>

#include "Tolerance.hpp"

using namespace def::algebra::utils;

template<class ScalarType>
class MatrixListWrapper;

template<class ScalarType>
MatrixListWrapper<ScalarType> operator+(const MatrixListWrapper<ScalarType> &left,
                                        const MatrixListWrapper<ScalarType> &right);
template<class ScalarType>
MatrixListWrapper<ScalarType> operator-(const MatrixListWrapper<ScalarType> &left,
                                        const MatrixListWrapper<ScalarType> &right);

template<class ScalarType>
MatrixListWrapper<ScalarType> operator*(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                        ScalarType const &rightScalar);
template<class ScalarType>
MatrixListWrapper<ScalarType> operator*(ScalarType const &leftScalar,
                                        const MatrixListWrapper<ScalarType> &rightMatrixListType);
template<class ScalarType>
MatrixListWrapper<ScalarType> operator/(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                        ScalarType const &rightScalar);

template<class ScalarType>
MatrixListWrapper<ScalarType> operator*(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                        int const &rightScalar);
template<class ScalarType>
MatrixListWrapper<ScalarType> operator*(int const &leftScalar,
                                        const MatrixListWrapper<ScalarType> &rightMatrixListType);
template<class ScalarType>
MatrixListWrapper<ScalarType> operator/(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                        int const &rightScalar);

/**
 *  \brief      Mathematical matrix list class.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica X.X
 *
 *  \details    The MatrixListWrapper class is a wrapping for matrix lists, adding basic operations to the
 *  			std::vector<MatrixType> structure.
 *
 */
template<class ScalarType>
class MatrixListWrapper : public Tolerance<MatrixListWrapper<ScalarType>> {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Friend methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  friend MatrixListWrapper<ScalarType> operator-<>(const MatrixListWrapper<ScalarType> &left,
                                                   const MatrixListWrapper<ScalarType> &right);
  friend MatrixListWrapper<ScalarType> operator+<>(const MatrixListWrapper<ScalarType> &left,
                                                   const MatrixListWrapper<ScalarType> &right);

  friend MatrixListWrapper<ScalarType> operator*<>(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                                   ScalarType const &rightScalar);
  friend MatrixListWrapper<ScalarType> operator*<>(ScalarType const &leftScalar,
                                                   const MatrixListWrapper<ScalarType> &rightMatrixListType);
  friend MatrixListWrapper<ScalarType> operator/<>(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                                   ScalarType const &rightScalar);

  friend MatrixListWrapper<ScalarType> operator*<>(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                                   int const &rightScalar);
  friend MatrixListWrapper<ScalarType> operator*<>(int const &leftScalar,
                                                   const MatrixListWrapper<ScalarType> &rightMatrixListType);
  friend MatrixListWrapper<ScalarType> operator/<>(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                                   int const &rightScalar);

 public :

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Typedefs :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Vector type.
  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  /// Matrix type.
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  /// Raw matrix list type.
  typedef std::vector<MatrixType> RawMatrixListType;


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
    typedef typename MatrixType::iterator MatrixIteratorType;

    explicit iterator(RawMatrixListType &ml, bool const &endInitialization) : m_RawMatrixList(ml) {
      if (endInitialization) {
        m_ListIterator = ml.size();
        m_MatrixIterator = ml[ml.size() - 1].end();
      } else {
        m_ListIterator = 0;
        m_MatrixIterator = ml[0].begin();
      }
    }

    iterator &operator++() {
      ++m_MatrixIterator;
      if (m_MatrixIterator == m_RawMatrixList[m_ListIterator].end()) {
        ++m_ListIterator;
        if (m_ListIterator != m_RawMatrixList.size())
          m_MatrixIterator = m_RawMatrixList[m_ListIterator].begin();
      }
      return *this;
    }

    bool operator==(iterator const &other) const {
      return ((m_MatrixIterator == other.m_MatrixIterator) && (m_ListIterator == other.m_ListIterator));
    }
    bool operator!=(iterator const &other) const { return !(*this == other); }

    ScalarType &operator*() const { return *m_MatrixIterator; }

   private:
    RawMatrixListType &m_RawMatrixList;
    MatrixIteratorType m_MatrixIterator;
    unsigned int m_ListIterator;
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
    explicit const_iterator(RawMatrixListType const &ml, bool const &endInitialization = 0) : m_RawMatrixList(ml) {
      if (endInitialization) {
        const unsigned int length = ml.size();
        m_ListIterator = length;
        m_ColIterator = ml[length - 1].cols();
        m_RowIterator = ml[length - 1].rows();
      } else {
        m_ListIterator = 0;
        m_ColIterator = 0;
        m_RowIterator = 0;
      }
    }

    const_iterator &operator++() {
      ++m_RowIterator;
      if (m_RowIterator == m_RawMatrixList[m_ListIterator].rows()) {
        ++m_ColIterator;
        if (m_ColIterator != m_RawMatrixList[m_ListIterator].cols()) {
          m_RowIterator = 0;
        } else {
          ++m_ListIterator;
          if (m_ListIterator != m_RawMatrixList.size()) {
            m_RowIterator = 0;
            m_ColIterator = 0;
          }
        }
      }
      return *this;
    }

    bool operator==(const_iterator const &other) const {
      return ((m_RowIterator == other.m_RowIterator) && (m_ColIterator == other.m_ColIterator)
          && (m_ListIterator == other.m_ListIterator));
    }
    bool operator!=(const_iterator const &other) const { return !(*this == other); }

    ScalarType const &operator*() const { return m_RawMatrixList[m_ListIterator](m_RowIterator, m_ColIterator); }

   private:
    RawMatrixListType m_RawMatrixList;

    unsigned int m_RowIterator;
    unsigned int m_ColIterator;
    unsigned int m_ListIterator;
  };


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Default constructor.
  MatrixListWrapper() {}

  /// Constructor from a raw matrix list.
  MatrixListWrapper(RawMatrixListType const &ml) {
    m_RawMatrixList.resize(ml.size());
    for (unsigned int k = 0; k < ml.size(); k++)
      m_RawMatrixList[k] = ml[k];
  }

  /// Constructor with an input matrix list size.
  MatrixListWrapper(const unsigned int &n) { m_RawMatrixList.resize(n); }

  /// Copy constructor.
  MatrixListWrapper(const MatrixListWrapper<ScalarType> &other) { m_RawMatrixList = other.m_RawMatrixList; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Methods :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Converts the matrix list back to a raw matrix list.
  RawMatrixListType GetRawMatrixList() const { return (m_RawMatrixList); }

  /// Returns the length of the list.
  unsigned int size() const { return m_RawMatrixList.size(); }
  /// Resizes the list to a length \e n.
  void resize(const unsigned int n) { m_RawMatrixList.resize(n); }

  /// Returns a reversed copy of the matrix list.
  MatrixListWrapper<ScalarType> reverse() const {
    RawMatrixListType out = m_RawMatrixList;
    std::reverse(out.begin(), out.end());
    return MatrixListWrapper<ScalarType>(out);
  }

  /// Returns a sub-MatrixList starting at index i1 and ending at i2.
  MatrixListWrapper<ScalarType> sublist(const unsigned int i1, const unsigned int i2) const {
    if (m_RawMatrixList.size() == 0) { return MatrixListWrapper<ScalarType>(m_RawMatrixList); }
    else if (i2 >= i1) {
      RawMatrixListType out(i2 - i1 + 1);
      for (unsigned int k = 0; k < i2 - i1 + 1; ++k) { out[k] = m_RawMatrixList[i1 + k]; }
      return MatrixListWrapper<ScalarType>(out);
    } else { return this->sublist(i2, i1).reverse(); }
  }

  /// Vectorizes the matrix list, matrix by matrix and row by row.
  VectorType vectorize() const {
    VectorType out = m_RawMatrixList[0].vectorise_row_wise();
    for (unsigned k = 1; k < m_RawMatrixList.size(); ++k)
      out.push_back(m_RawMatrixList[k].vectorise_row_wise());
    return out;
  }
  /// Returns the size parameters of the matrix list.
  std::vector<unsigned int> GetSizeParameters() const {
    const unsigned int nbOfMatrices = m_RawMatrixList.size();
    std::vector<unsigned int> out(2 * nbOfMatrices + 3); // Offset by 3 to discriminate the cases nbOfMatrices = 0 or 1.
    for (unsigned int k = 0; k < nbOfMatrices; ++k) {
      out[2 * k] = m_RawMatrixList[k].rows();
      out[2 * k + 1] = m_RawMatrixList[k].cols();
    }
    return out;
  }

  /// Returns the number of elements.
  unsigned int n_elem() const {
    unsigned int result(0);
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      result += m_RawMatrixList[k].n_elem();
    return result;
  }

  /// Fills the matrix list.
  void fill(ScalarType const &value) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k].fill(value);
  }

  /// Returns an iterator on the first scalar element of the raw matrix list.
  iterator begin() { return iterator(m_RawMatrixList, 0); }
  const_iterator begin() const { return const_iterator(m_RawMatrixList, 0); }

  /// Returns an iterator on the end of the raw matrix list.
  iterator end() { return iterator(m_RawMatrixList, 1); }
  const_iterator end() const { return const_iterator(m_RawMatrixList, 1); }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Arithmetic operations :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the sum of all elements.
  ScalarType sum() const {
    ScalarType result(0);
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      result += m_RawMatrixList[k].sum();
    return result;
  }

  /// Returns the sum of squares of all elements.
  ScalarType sum_of_squares() const {
    ScalarType result(0);
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      result += m_RawMatrixList[k].sum_of_squares();
    return result;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Operator overloading :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Access the \e i th matrix of the list for reading.
  MatrixType const &operator[](const unsigned int i) const { return m_RawMatrixList[i]; }
  /// Access the \e i th matrix of the list for reading or writing.
  MatrixType &operator[](const unsigned int i) { return m_RawMatrixList[i]; }

  /// Access the \e ith matrix of the list for reading. Performs boundary checks.
  MatrixType const &at(const unsigned int i) const { return m_RawMatrixList.at(i); }
  /// Access the \e ith matrix of the list for reading or writing. Performs boundary checks.
  MatrixType &at(const unsigned int i) { return m_RawMatrixList.at(i); }

  /// Add rhs scalar to lhs matrix list in situ.
  MatrixListWrapper<ScalarType> &operator+=(ScalarType rhsScalar) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k] += rhsScalar;
    return *this;
  }

  /// Subtract rhs scalar from lhs matrix list in situ.
  MatrixListWrapper<ScalarType> &operator-=(ScalarType rhsScalar) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k] -= rhsScalar;
    return *this;
  }

  /// Scalar multiplication in situ of lhs matrix list by \e rhsScalar.
  MatrixListWrapper<ScalarType> &operator*=(ScalarType rhsScalar) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k] *= rhsScalar;
    return *this;
  }

  /// Scalar division of lhs matrix list in situ by \e rhsScalar.
  MatrixListWrapper<ScalarType> &operator/=(ScalarType rhsScalar) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k] /= rhsScalar;
    return *this;
  }

  /// Add rhs matrix list to lhs matrix list in situ.
  MatrixListWrapper<ScalarType> &operator+=(MatrixListWrapper<ScalarType> const &rhsMatrixListType) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k] += rhsMatrixListType.m_RawMatrixList[k];
    return *this;
  }

  /// Subtract rhs matrix list from lhs matrix list in situ.
  MatrixListWrapper<ScalarType> &operator-=(MatrixListWrapper<ScalarType> const &rhsMatrixListType) {
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      m_RawMatrixList[k] -= rhsMatrixListType.m_RawMatrixList[k];
    return *this;
  }

  /// Unary minus operator.
  MatrixListWrapper<ScalarType> operator-() const {
    RawMatrixListType rawMatrixList(m_RawMatrixList.size());
    for (unsigned int k = 0; k < m_RawMatrixList.size(); k++)
      rawMatrixList[k] = -m_RawMatrixList[k];
    return MatrixListWrapper<ScalarType>(rawMatrixList);
  }

  const bool operator==(const MatrixListWrapper<ScalarType> &t) const {
    if (m_RawMatrixList.size() != t.m_RawMatrixList.size()) return false;

    return std::equal(m_RawMatrixList.begin(), m_RawMatrixList.end(), t.m_RawMatrixList.begin());
  }

  virtual const bool compare(const MatrixListWrapper<ScalarType>& t, const float tolerance) {
    if (m_RawMatrixList.size() != t.m_RawMatrixList.size()) return false;

    return this->equal(m_RawMatrixList.begin(), m_RawMatrixList.end(), t.m_RawMatrixList.begin(), tolerance);
  }

 private:

  template<class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & m_RawMatrixList;
  }
  friend class boost::serialization::access;

  RawMatrixListType m_RawMatrixList;

};

template<class ScalarType>
inline ScalarType dot_product(MatrixListWrapper<ScalarType> const &left,
                              MatrixListWrapper<ScalarType> const &right) {
  ScalarType result = 0;
  for (unsigned int k = 0; k < left.size(); k++)
    result += dot_product(left[k], right[k]);
  return result;
}

template<class ScalarType>
inline MatrixListWrapper<ScalarType> concatenate(MatrixListWrapper<ScalarType> const &left,
                                                 MatrixListWrapper<ScalarType> const &right) {
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  const std::vector<MatrixType> rml_left = left.GetRawMatrixList();
  const std::vector<MatrixType> rml_right = right.GetRawMatrixList();
  const unsigned int rml_left_size = rml_left.size();
  const unsigned int rml_right_size = rml_right.size();

  std::vector<MatrixType> out(rml_left_size + rml_right_size);
  for (unsigned int k = 0; k < rml_left_size; ++k) { out[k] = rml_left[k]; }
  for (unsigned int k = 0; k < rml_right_size; ++k) { out[rml_left_size + k] = rml_right[k]; }
  return MatrixListType(out);
}

//#include "MatrixListWrapper_friend.h"

#endif /* _MatrixListWrapper_h */
