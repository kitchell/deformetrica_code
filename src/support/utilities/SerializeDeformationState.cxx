/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "SerializeDeformationState.h"
namespace def {
namespace utils {

DeformationState serialize;

std::ostream& operator<<(std::ostream& os, const DeformationState& ser) {
#define PR(t,n) os << std::endl << "##n" << std::endl; std::copy(n.begin(),n.end(),std::ostream_iterator<t>(os," "));
  PR(bool, ser.bool_type_)
  PR(int, ser.int_type_)
  PR(unsigned int, ser.uint_type_)
  return os;
}



}
}