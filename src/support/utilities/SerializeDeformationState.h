/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#pragma once

#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/version.hpp>
#include <src/support/linear_algebra/LinearAlgebra.h>
#include <src/support/linear_algebra/Tolerance.hpp>
#include <deque>
#include <vector>

#define OPERATOR_MACRO(type, name) \
  DeformationState& operator<<(const type &t) { \
    name.push_back(t); \
    return *this; \
  } \
  DeformationState& operator>>(type &t) { \
    t = name.front(); \
    name.pop_front(); \
    return *this; \
  }

using namespace def::algebra;
using namespace def::algebra::utils;

namespace def {
namespace utils {


class DeformationState : public Tolerance<DeformationState> {
 public:
  typedef std::vector<std::string> VecStringType;
  typedef std::vector<int> VecIntType;
  typedef std::vector<unsigned int> VecUIntType;
  typedef std::vector<float> VecFloatType;
  typedef std::vector<double> VecDoubleType;
  typedef std::vector<VectorType> VecVectorType;
  typedef std::vector<MatrixType> VecMatrixType;
  typedef std::vector<MatrixListType> VecMatrixListType;

  DeformationState(){}

  void reset() {
    DeformationState *df = new DeformationState;
    std::swap(*this, *df);
    delete df;
  }

  void save(const std::string &file) {
    std::ofstream ofs(file);

    boost::archive::text_oarchive oa(ofs);
    oa << *this;
  }

  void save_and_reset(const std::string &file) {
    std::ofstream ofs(file);

    boost::archive::text_oarchive oa(ofs);
    oa << *this;
    reset();
  }

  void load(const std::string &file) {
    std::ifstream ifs(file);

    DeformationState *df = new DeformationState;
    boost::archive::text_iarchive ia(ifs);
    ia >> *df;
    std::swap(*this, *df);
    delete df;
  }

  void deformation(DeformationState&& t) {
    *this = std::move(t);
  }

  OPERATOR_MACRO(bool,bool_type_);
  OPERATOR_MACRO(int,int_type_);
  OPERATOR_MACRO(unsigned int,uint_type_);
  OPERATOR_MACRO(float,float_type_);
  OPERATOR_MACRO(double,double_type_);
  OPERATOR_MACRO(std::string,string_type_);
  OPERATOR_MACRO(VectorType,vector_type_);
  OPERATOR_MACRO(MatrixType,matrix_type_);
  OPERATOR_MACRO(MatrixListType,matrix_list_type_);
  OPERATOR_MACRO(LinearVariableType,linear_variable_type_);
  OPERATOR_MACRO(LinearVariableMapType,linear_variable_map_type_);
  OPERATOR_MACRO(LinearVariablesType,linear_variables_type_);
  OPERATOR_MACRO(LinearVariablesMapType,linear_variables_map_type_);
  OPERATOR_MACRO(VecStringType, vec_string_type_);
  OPERATOR_MACRO(VecIntType, vec_int_type_);
  OPERATOR_MACRO(VecUIntType, vec_uint_type_);
  OPERATOR_MACRO(VecFloatType, vec_float_type_);
  OPERATOR_MACRO(VecDoubleType, vec_double_type_);
  OPERATOR_MACRO(VecVectorType, vec_vector_type_);
  OPERATOR_MACRO(VecMatrixType, vec_matrix_type_);
  OPERATOR_MACRO(VecMatrixListType, vec_matrix_list_type_);

  const bool operator==(const DeformationState& ds) const {
#define EQ(t) (t.size() == ds.t.size() && std::equal(t.begin(), t.end(), ds.t.begin()))
    return  EQ(bool_type_)  &&
            EQ(int_type_)   &&
            EQ(uint_type_)   &&
            EQ(float_type_)   &&
            EQ(double_type_)   &&
            EQ(string_type_)   &&
            EQ(vector_type_)   &&
            EQ(matrix_type_)   &&
            EQ(matrix_list_type_)   &&
            EQ(linear_variable_type_)   &&
            EQ(linear_variable_map_type_)   &&
            EQ(linear_variables_type_)   &&
            EQ(linear_variables_map_type_)   &&
            EQ(vec_string_type_)   &&
            EQ(vec_int_type_)   &&
            EQ(vec_uint_type_)   &&
            EQ(vec_float_type_)   &&
            EQ(vec_double_type_)   &&
            EQ(vec_vector_type_)   &&
            EQ(vec_matrix_type_)   &&
            EQ(vec_matrix_list_type_);
#undef EQ
  }


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

  template<class InputIt1, class InputIt2>
  bool vequal(InputIt1 first1, InputIt1 last1,
             InputIt2 first2, const float tolerance) const
  {
    for (; first1 != last1; ++first1, ++first2) {
      if (!equal(first1->begin(),first1->end(),first2->begin(),tolerance)) {
        return false;
      }
    }
    return true;
  }

  virtual const bool compare(const DeformationState& ds, const float tolerance) {
    ScalarCompare compare(tolerance);
#define EQ(t) (t.size() == ds.t.size() && std::equal(t.begin(), t.end(), ds.t.begin()))
#define TO(t) (t.size() == ds.t.size() && std::equal(t.begin(), t.end(), ds.t.begin(), compare))
#define CO(t) (t.size() == ds.t.size() && equal(t.begin(), t.end(), ds.t.begin(), tolerance))
#define COV(t) (t.size() == ds.t.size() && vequal(t.begin(), t.end(), ds.t.begin(), tolerance))

    auto CHK = [](const bool b, const std::string& msg) {
#ifndef NDEBUG
      if (!b) std::cerr << std::endl << msg << " type mismatch" << std::endl;
#endif
      return b;
    };

    bool ret = true;
    ret = CHK(EQ(bool_type_),"bool") && ret;
    ret = CHK(EQ(int_type_),"int") && ret;
    ret = CHK(EQ(uint_type_),"uint") && ret;
    ret = CHK(TO(float_type_),"float") && ret;
    ret = CHK(TO(double_type_),"double") && ret;
    ret = CHK(EQ(string_type_),"string") && ret;
    ret = CHK(CO(vector_type_),"vector") && ret;
    ret = CHK(CO(matrix_type_),"matrix") && ret;
    ret = CHK(CO(matrix_list_type_),"matrix list") && ret;
    ret = CHK(CO(linear_variable_type_),"linear variable") && ret;
    ret = CHK(CO(linear_variable_map_type_),"linear map variable") && ret;
    ret = CHK(CO(linear_variables_type_),"linear variables") && ret;
    ret = CHK(CO(linear_variables_map_type_),"linear map variables") && ret;
    ret = CHK(EQ(vec_string_type_),"vec string") && ret;
    ret = CHK(EQ(vec_int_type_),"vec int") && ret;
    ret = CHK(TO(vec_float_type_),"vec float") && ret;
    ret = CHK(TO(vec_double_type_),"vec double") && ret;
    ret = CHK(COV(vec_vector_type_),"vec vector") && ret;
    ret = CHK(COV(vec_matrix_type_),"vec matrix") && ret;
    ret = CHK(COV(vec_matrix_list_type_),"vec matrix list") && ret;

    return ret;
#undef EQ
#undef TO
#undef CO
  }

  friend std::ostream& operator<<(std::ostream& os, const DeformationState& ser);

    private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & bool_type_;
    ar & int_type_;
    ar & uint_type_;
    ar & float_type_;
    ar & double_type_;
    ar & string_type_;
    ar & vector_type_;
    ar & matrix_type_;
    ar & matrix_list_type_;
    ar & linear_variable_type_;
    ar & linear_variable_map_type_;
    ar & linear_variables_type_;
    ar & linear_variables_map_type_;
    ar & vec_string_type_;
    ar & vec_int_type_;
    ar & vec_uint_type_;
    ar & vec_float_type_;
    ar & vec_double_type_;
    ar & vec_vector_type_;
    ar & vec_matrix_type_;
    ar & vec_matrix_list_type_;
  }

 private:
  std::deque<bool> bool_type_;
  std::deque<int> int_type_;
  std::deque<unsigned int> uint_type_;
  std::deque<float> float_type_;
  std::deque<double> double_type_;
  std::deque<std::string> string_type_;
  std::deque<VectorType> vector_type_;
  std::deque<MatrixType> matrix_type_;
  std::deque<MatrixListType> matrix_list_type_;
  std::deque<LinearVariableType> linear_variable_type_;
  std::deque<LinearVariableMapType> linear_variable_map_type_;
  std::deque<LinearVariablesType> linear_variables_type_;
  std::deque<LinearVariablesMapType> linear_variables_map_type_;
  std::deque<VecStringType> vec_string_type_;
  std::deque<VecIntType> vec_int_type_;
  std::deque<VecUIntType> vec_uint_type_;
  std::deque<VecFloatType> vec_float_type_;
  std::deque<VecDoubleType> vec_double_type_;
  std::deque<VecVectorType> vec_vector_type_;
  std::deque<VecMatrixType> vec_matrix_type_;
  std::deque<VecMatrixListType> vec_matrix_list_type_;
};


class SerializeDeformationState {
 public:
  SerializeDeformationState(const SerializeDeformationState&) = delete;
  SerializeDeformationState(SerializeDeformationState&&) = delete;

  SerializeDeformationState& operator=(const SerializeDeformationState&) = delete;
  SerializeDeformationState& operator=(SerializeDeformationState&&) = delete;

  static SerializeDeformationState* instance() {
    static SerializeDeformationState* instance = nullptr;
    if (instance) return instance;
    return (instance = new SerializeDeformationState());
  }

  DeformationState& state() {
    return deformation_state_;
  }

 protected:
  SerializeDeformationState() {}
  ~SerializeDeformationState() {}

  DeformationState deformation_state_;
};

extern DeformationState serialize;

}
}

BOOST_CLASS_VERSION(def::utils::DeformationState, 1)
