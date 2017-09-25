/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "AutomaticRelevanceDeterminationDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
AutomaticRelevanceDeterminationDistribution<ScalarType>
::AutomaticRelevanceDeterminationDistribution() {
}

template<class ScalarType>
AutomaticRelevanceDeterminationDistribution<ScalarType>
::~AutomaticRelevanceDeterminationDistribution() {}

template<class ScalarType>
AutomaticRelevanceDeterminationDistribution<ScalarType>
::AutomaticRelevanceDeterminationDistribution(const AutomaticRelevanceDeterminationDistribution &other) : Superclass(other) {
  m_Gamma = other.m_Gamma;
}

template<class ScalarType>
ScalarType
AutomaticRelevanceDeterminationDistribution<ScalarType>
::ComputeLogLikelihood(MatrixType const &obs, ScalarType const &temperature) const{

  ScalarType out;
  unsigned int q = obs.cols();///Number of pca components
  assert(m_Gamma.size()==q);
  unsigned int d = obs.cols();///Size of the components i.e. dim of momenta space
  for (unsigned int i=0;i<q;++i)
    out += 0.5*d*std::log(m_Gamma[i]/(2*M_PI));


  ///Because we want to access the columns of the matrix, could be cleaner.
  for (unsigned int i;i<q;++i)
    out -= 0.5 * m_Gamma[i] * obs.get_column(i).sum_of_squares();

  return out;

}


template class AutomaticRelevanceDeterminationDistribution<double>;