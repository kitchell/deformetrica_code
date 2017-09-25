/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah. All rights reserved. This file is     *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _LinearVariableVariantVisitors_h
#define _LinearVariableVariantVisitors_h

/// Non-core files.
#include "ArmadilloVectorWrapper.h"
#include "ArmadilloMatrixWrapper.h"
#include "MatrixListWrapper.h"

/// Libraries files.
#include "boost/variant.hpp"

/**
 *	\brief      LinearVariableVariantVisitors objects
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    Needed tool classes to exploit the boost::variant data structure.
 */


template<class ScalarType>
class vectorize_visitor : public boost::static_visitor<ArmadilloVectorWrapper<ScalarType>> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  VectorType operator()(ScalarType const &var) const { return VectorType(1, var); }
  VectorType operator()(VectorType const &var) const { return var; }
  VectorType operator()(MatrixType const &var) const { return var.vectorize(); }
  VectorType operator()(MatrixListType const &var) const { return var.vectorize(); }

}; /* class vectorize_visitor */


template<class ScalarType>
class vectorize_with_parameters_visitor
    : public boost::static_visitor<ArmadilloVectorWrapper<ScalarType>> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Getter.
  std::vector<unsigned int> GetParams() const { return m_Params; }

  VectorType operator()(ScalarType const &var) {
    m_Params.resize(0);
    return VectorType(1, var);
  }
  VectorType operator()(VectorType const &var) {
    m_Params.resize(1);
    m_Params[0] = var.size();
    return var;
  }
  VectorType operator()(MatrixType const &var) {
    m_Params.resize(2);
    m_Params[0] = var.rows();
    m_Params[1] = var.cols();
    return var.vectorize();
  }
  VectorType operator()(MatrixListType const &var) {
    const unsigned int nbOfMatrices = var.size();
    m_Params.resize(2 * nbOfMatrices + 3); // Offset by 3 to discriminate the cases nbOfMatrices = 0 or 1.
    for (unsigned int k = 0; k < nbOfMatrices; ++k) {
      m_Params[2 * k] = var[k].rows();
      m_Params[2 * k + 1] = var[k].cols();
    }
    return var.vectorize();
  }

 private:
  std::vector<unsigned int> m_Params;

}; /* class vectorize_with_parameters_visitor */


template<class ScalarType>
class fill_visitor : public boost::static_visitor<> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Constructor.
  fill_visitor(ScalarType const &value) : m_FillValue(value) {}

  void operator()(ScalarType var) const { var = m_FillValue; }
  void operator()(VectorType var) const { var.fill(m_FillValue); }
  void operator()(MatrixType var) const { var.fill(m_FillValue); }
  void operator()(MatrixListType var) const { var.fill(m_FillValue); }

 private:
  ScalarType m_FillValue;

}; /* class fill_visitor */


template<class ScalarType>
class sum_visitor : public boost::static_visitor<ScalarType> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  ScalarType operator()(ScalarType const &var) const { return var; }
  ScalarType operator()(VectorType const &var) const { return var.sum(); }
  ScalarType operator()(MatrixType const &var) const { return var.sum(); }
  ScalarType operator()(MatrixListType const &var) const { return var.sum(); }

}; /* class sum_visitor */

template<class ScalarType>
class sum_of_squares_visitor : public boost::static_visitor<ScalarType> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  ScalarType operator()(ScalarType const &var) const { return var * var; }
  ScalarType operator()(VectorType const &var) const { return var.sum_of_squares(); }
  ScalarType operator()(MatrixType const &var) const { return var.sum_of_squares(); }
  ScalarType operator()(MatrixListType const &var) const { return var.sum_of_squares(); }

}; /* class sum_of_squares_visitor */


template<class ScalarType>
class dot_product_visitor : public boost::static_visitor<ScalarType> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  ScalarType operator()(ScalarType const &left, ScalarType const &right) const { return left * right; }
  ScalarType operator()(VectorType const &left, VectorType const &right) const { return dot_product(left, right); }
  ScalarType operator()(MatrixType const &left, MatrixType const &right) const { return dot_product(left, right); }
  ScalarType operator()(MatrixListType const &left, MatrixListType const &right) const {
    return dot_product(left, right);
  }

}; /* class dot_product_visitor */


template<class ScalarType>
class number_of_elements_visitor : public boost::static_visitor<ScalarType> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  ScalarType operator()(ScalarType const &var) const { return 1; }
  ScalarType operator()(VectorType const &var) const { return var.n_elem(); }
  ScalarType operator()(MatrixType const &var) const { return var.n_elem(); }
  ScalarType operator()(MatrixListType const &var) const { return var.n_elem(); }

}; /* class number_of_elements_visitor */


template<class ScalarType>
class addition_visitor : public boost::static_visitor<boost::variant<ScalarType,
                                                                     ArmadilloVectorWrapper<ScalarType>,
                                                                     ArmadilloMatrixWrapper<ScalarType>,
                                                                     MatrixListWrapper<ScalarType>>> {
 public:

  /// Raw linear variable type.
  typedef boost::variant<ScalarType, ArmadilloVectorWrapper<ScalarType>,
                         ArmadilloMatrixWrapper<ScalarType>, MatrixListWrapper<ScalarType>>
      RawLinearVariableType;

  template<typename VariantType1, typename VariantType2>
  RawLinearVariableType operator()(const VariantType1 &, const VariantType2 &) const {
    std::cerr << "Exception: attempt to add incompatible data structures during the manipulation of "
        "boost::variant<ScalarType, VectorType, MatrixType, MatrixListType> type." << std::endl;
    return RawLinearVariableType();
  }

  template<typename VariantType>
  RawLinearVariableType operator()(const VariantType &left, const VariantType &right) const {
    RawLinearVariableType result = left + right;
    return result;
  }

}; /* class addition_visitor */


template<class ScalarType>
class subtraction_visitor : public boost::static_visitor<boost::variant<ScalarType,
                                                                        ArmadilloVectorWrapper<ScalarType>,
                                                                        ArmadilloMatrixWrapper<ScalarType>,
                                                                        MatrixListWrapper<ScalarType>>> {
 public:

  /// Raw linear variable type.
  typedef boost::variant<ScalarType, ArmadilloVectorWrapper<ScalarType>,
                         ArmadilloMatrixWrapper<ScalarType>, MatrixListWrapper<ScalarType>>
      RawLinearVariableType;

  template<typename VariantType1, typename VariantType2>
  RawLinearVariableType operator()(const VariantType1 &, const VariantType2 &) const {
    std::cerr << "Exception: attempt to subtract incompatible data structures during the manipulation of "
        "boost::variant<ScalarType, VectorType, MatrixType, MatrixListType> type." << std::endl;
    return RawLinearVariableType();
  }

  template<typename VariantType>
  RawLinearVariableType operator()(const VariantType &left, const VariantType &right) const {
    RawLinearVariableType result = left - right;
    return result;
  }

}; /* class subtraction_visitor */


template<class ScalarType>
class right_scalar_multiplication_visitor : public boost::static_visitor<
    boost::variant<ScalarType,
                   ArmadilloVectorWrapper<ScalarType>,
                   ArmadilloMatrixWrapper<ScalarType>,
                   MatrixListWrapper<ScalarType>>> {
 public:

  /// Raw linear variable type.
  typedef boost::variant<ScalarType, ArmadilloVectorWrapper<ScalarType>,
                         ArmadilloMatrixWrapper<ScalarType>, MatrixListWrapper<ScalarType>>
      RawLinearVariableType;

  /// Constructor.
  right_scalar_multiplication_visitor(ScalarType const &scalar) : m_Right(scalar) {}

  template<typename VariantType>
  RawLinearVariableType operator()(const VariantType &left) const {
    RawLinearVariableType result = left * m_Right;
    return result;
  }

 private:
  ScalarType m_Right;

}; /* class right_scalar_multiplication_visitor */


template<class ScalarType>
class insitu_scalar_addition_visitor : public boost::static_visitor<> {
 public:

  /// Constructor.
  insitu_scalar_addition_visitor(ScalarType const &scalar) : m_Right(scalar) {}

  template<typename VariantType>
  void operator()(VariantType &left) const { left += m_Right; }

 private:
  ScalarType m_Right;

}; /* class insitu_scalar_addition_visitor */


template<class ScalarType>
class insitu_scalar_multiplication_visitor : public boost::static_visitor<> {
 public:

  /// Constructor.
  insitu_scalar_multiplication_visitor(ScalarType const &scalar) : m_Right(scalar) {}

  template<typename VariantType>
  void operator()(VariantType &left) const { left *= m_Right; }

 private:
  ScalarType m_Right;

}; /* class insitu_scalar_multiplication_visitor */


template<class ScalarType>
class insitu_addition_visitor : public boost::static_visitor<> {
 public:

  template<typename VariantType>
  void operator()(VariantType &left, const VariantType &right) const {
    left += right;
  }

  template<typename VariantType1, typename VariantType2>
  void operator()(VariantType1 &var1, const VariantType2 &var2) const {
    std::cerr << "Exception: attempt to add in situ incompatible data structures during the manipulation"
        " of boost::variant<ScalarType, VectorType, MatrixType, MatrixListType> type." << std::endl;
  }

}; /* class insitu_addition_visitor */


template<class ScalarType>
class insitu_subtraction_visitor : public boost::static_visitor<> {
 public:

  template<typename VariantType1, typename VariantType2>
  void operator()(VariantType1 &, const VariantType2 &) const {
    std::cerr << "Exception: attempt to subtract in situ incompatible data structures during the manipulation"
        " of boost::variant<ScalarType, VectorType, MatrixType, MatrixListType> type." << std::endl;
  }

  template<typename VariantType>
  void operator()(VariantType &left, const VariantType &right) const {
    left -= right;
  }

}; /* class insitu_subtraction_visitor */


template<class ScalarType>
class unary_minus_visitor : public boost::static_visitor<boost::variant<ScalarType,
                                                                        ArmadilloVectorWrapper<ScalarType>,
                                                                        ArmadilloMatrixWrapper<ScalarType>,
                                                                        MatrixListWrapper<ScalarType>>> {
 public:

  /// Raw linear variable type.
  typedef boost::variant<ScalarType, ArmadilloVectorWrapper<ScalarType>,
                         ArmadilloMatrixWrapper<ScalarType>,
                         MatrixListWrapper<ScalarType>> RawLinearVariableType;

  template<typename VariantType>
  RawLinearVariableType operator()(const VariantType &var) const {
    RawLinearVariableType result = -var;
    return result;
  }

}; /* class unary_minus_visitor */


////////////////////////////////////////////////////////////////////////////////////////////////////
// Iterator :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
class iterator_increment_visitor : public boost::static_visitor<> {
 public:

  typedef std::vector<ScalarType *> ScalarIteratorType;
  typedef typename ArmadilloMatrixWrapper<ScalarType>::iterator ArmadilloIteratorType;
  typedef typename MatrixListWrapper<ScalarType>::iterator MatrixListIteratorType;

  void operator()(ScalarIteratorType &it) const { it.push_back(0); }
  void operator()(ArmadilloIteratorType &it) const { ++it; }
  void operator()(MatrixListIteratorType &it) const { ++it; }

}; /* class iterator_increment_visitor */


template<class ScalarType>
class iterator_equality_visitor : public boost::static_visitor<bool> {
 public:

  typedef std::vector<ScalarType *> ScalarIteratorType;
  typedef typename ArmadilloMatrixWrapper<ScalarType>::iterator ArmadilloIteratorType;
  typedef typename MatrixListWrapper<ScalarType>::iterator MatrixListIteratorType;

  bool operator()(ScalarIteratorType const &it1, ScalarIteratorType const &it2) const {
    return it1.size() == it2.size();
  }
  bool operator()(ArmadilloIteratorType const &it1, ArmadilloIteratorType const &it2) const { return it1 == it2; }
  bool operator()(MatrixListIteratorType const &it1, MatrixListIteratorType const &it2) const { return it1 == it2; }

  template<typename VariantType1, typename VariantType2>
  bool operator()(VariantType1 const &it1, VariantType2 const &it2) const {
    std::cout << "Warning: attempt to compare different iterator types during the manipulation of "
        "LinearVariableWrapper::iterator." << std::endl;
    return 0;
  }

}; /* class iterator_equality_visitor */


template<class ScalarType>
class iterator_dereferencing_visitor : public boost::static_visitor<ScalarType &> {
 public:

  typedef std::vector<ScalarType *> ScalarIteratorType;
  typedef typename ArmadilloMatrixWrapper<ScalarType>::iterator ArmadilloIteratorType;
  typedef typename MatrixListWrapper<ScalarType>::iterator MatrixListIteratorType;

  ScalarType &operator()(ScalarIteratorType const &it) const { return *it[0]; }
  ScalarType &operator()(ArmadilloIteratorType const &it) const { return *it; }
  ScalarType &operator()(MatrixListIteratorType const &it) const { return *it; }

}; /* class iterator_dereferencing_visitor */


template<class ScalarType>
class iterator_begin_visitor : public boost::static_visitor<boost::variant<std::vector<ScalarType *>,
                                                                           typename ArmadilloMatrixWrapper<ScalarType>::iterator,
                                                                           typename MatrixListWrapper<ScalarType>::iterator>> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Raw linear variable iterator type.
  typedef boost::variant<std::vector<ScalarType *>, typename MatrixType::iterator,
                         typename MatrixListType::iterator> RawLinearVariableIteratorType;

  RawLinearVariableIteratorType operator()(ScalarType &var) const {
    std::vector<ScalarType *> it;
    it.push_back(&var);
    return it;
  }
  RawLinearVariableIteratorType operator()(VectorType &var) const { return var.begin(); }
  RawLinearVariableIteratorType operator()(MatrixType &var) const { return var.begin(); }
  RawLinearVariableIteratorType operator()(MatrixListType &var) const { return var.begin(); }

}; /* class iterator_begin_visitor */


template<class ScalarType>
class iterator_end_visitor : public boost::static_visitor<boost::variant<std::vector<ScalarType *>,
                                                                         typename ArmadilloMatrixWrapper<ScalarType>::iterator,
                                                                         typename MatrixListWrapper<ScalarType>::iterator>> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Raw linear variable iterator type.
  typedef boost::variant<std::vector<ScalarType *>, typename MatrixType::iterator,
                         typename MatrixListType::iterator> RawLinearVariableIteratorType;

  RawLinearVariableIteratorType operator()(ScalarType &var) const {
    std::vector<ScalarType *> it;
    it.push_back(&var);
    it.push_back(0);
    return it;
  }
  RawLinearVariableIteratorType operator()(VectorType &var) const { return var.end(); }
  RawLinearVariableIteratorType operator()(MatrixType &var) const { return var.end(); }
  RawLinearVariableIteratorType operator()(MatrixListType &var) const { return var.end(); }

}; /* class iterator_end_visitor */


////////////////////////////////////////////////////////////////////////////////////////////////////
// Const iterator :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
class const_iterator_increment_visitor : public boost::static_visitor<> {
 public:

  typedef std::vector<ScalarType> ScalarIteratorType;
  typedef typename ArmadilloMatrixWrapper<ScalarType>::const_iterator ArmadilloIteratorType;
  typedef typename MatrixListWrapper<ScalarType>::const_iterator MatrixListIteratorType;

  void operator()(ScalarIteratorType &it) const { it.push_back(0); }
  void operator()(ArmadilloIteratorType &it) const { ++it; }
  void operator()(MatrixListIteratorType &it) const { ++it; }

}; /* class const_iterator_increment_visitor */


template<class ScalarType>
class const_iterator_equality_visitor : public boost::static_visitor<bool> {
 public:

  typedef std::vector<ScalarType> ScalarIteratorType;
  typedef typename ArmadilloMatrixWrapper<ScalarType>::const_iterator ArmadilloIteratorType;
  typedef typename MatrixListWrapper<ScalarType>::const_iterator MatrixListIteratorType;

  bool operator()(ScalarIteratorType const &it1, ScalarIteratorType const &it2) const {
    return it1.size() == it2.size();
  }
  bool operator()(ArmadilloIteratorType const &it1, ArmadilloIteratorType const &it2) const { return it1 == it2; }
  bool operator()(MatrixListIteratorType const &it1, MatrixListIteratorType const &it2) const { return it1 == it2; }

  template<typename VariantType1, typename VariantType2>
  bool operator()(VariantType1 const &it1, VariantType2 const &it2) const {
    std::cout << "Warning: attempt to compare different iterator types during the manipulation of "
        "LinearVariableWrapper::const_iterator." << std::endl;
    return 0;
  }

}; /* class const_iterator_equality_visitor */


template<class ScalarType>
class const_iterator_dereferencing_visitor : public boost::static_visitor<ScalarType const &> {
 public:

  typedef std::vector<ScalarType> ScalarIteratorType;
  typedef typename ArmadilloMatrixWrapper<ScalarType>::const_iterator ArmadilloIteratorType;
  typedef typename MatrixListWrapper<ScalarType>::const_iterator MatrixListIteratorType;

  ScalarType const &operator()(ScalarIteratorType const &it) const { return it[0]; }
  ScalarType const &operator()(ArmadilloIteratorType const &it) const { return *it; }
  ScalarType const &operator()(MatrixListIteratorType const &it) const { return *it; }

}; /* class const_iterator_dereferencing_visitor */


template<class ScalarType>
class const_iterator_begin_visitor : public boost::static_visitor<
    boost::variant<std::vector<ScalarType>,
                   typename ArmadilloMatrixWrapper<ScalarType>::const_iterator,
                   typename MatrixListWrapper<ScalarType>::const_iterator>> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Raw linear variable const iterator type.
  typedef boost::variant<std::vector<ScalarType>, typename MatrixType::const_iterator,
                         typename MatrixListType::const_iterator> RawLinearVariableIteratorType;

  RawLinearVariableIteratorType operator()(ScalarType const &var) const {
    std::vector<ScalarType> it;
    it.push_back(var);
    return it;
  }
  RawLinearVariableIteratorType operator()(VectorType const &var) const { return var.begin(); }
  RawLinearVariableIteratorType operator()(MatrixType const &var) const { return var.begin(); }
  RawLinearVariableIteratorType operator()(MatrixListType const &var) const { return var.begin(); }

}; /* class const_iterator_begin_visitor */


template<class ScalarType>
class const_iterator_end_visitor : public boost::static_visitor<
    boost::variant<std::vector<ScalarType>,
                   typename ArmadilloMatrixWrapper<ScalarType>::const_iterator,
                   typename MatrixListWrapper<ScalarType>::const_iterator>> {
 public:

  typedef ArmadilloVectorWrapper<ScalarType> VectorType;
  typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;
  typedef MatrixListWrapper<ScalarType> MatrixListType;

  /// Raw linear variable const iterator type.
  typedef boost::variant<std::vector<ScalarType>, typename MatrixType::const_iterator,
                         typename MatrixListType::const_iterator> RawLinearVariableIteratorType;

  RawLinearVariableIteratorType operator()(ScalarType const &var) const {
    std::vector<ScalarType> it;
    it.push_back(var);
    it.push_back(0);
    return it;
  }
  RawLinearVariableIteratorType operator()(VectorType const &var) const { return var.end(); }
  RawLinearVariableIteratorType operator()(MatrixType const &var) const { return var.end(); }
  RawLinearVariableIteratorType operator()(MatrixListType const &var) const { return var.end(); }

}; /* class const_iterator_end_visitor */


#endif /* _LinearVariableVariantVisitors_h */
