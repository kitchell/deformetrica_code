/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

namespace def {
namespace io {
class XmlModel;
}
};

class SparseDiffeoParameters;

#include <memory>
#include "Regression.h"

/**
 estimateGeodesicRegression()
 Manages the estimation of sparse geodesic regression.  Reads and organizes input and handles output of results.

 @param[in] paramDiffeos An object containing diffeo parameters
 @param[in] numObservations The number of observations per object
 @param[in] numObjects The number of unique objects (eg 2 for multiple observations of left/right hemisphere)
 @param[in] paramObjectsList A vector of objects containing parameter values for each object
 @param[in] templatefnList A vector of vectors of string paths to data defining the inital template for each observation and object
 @param[in] observationfnList A vector of vectors of string paths to data defining the observations for each observation and object
 @param[in] observationTimesList A vector of vectors of times for each observation and object
 */

template<unsigned int Dimension>
void estimateGeodesicRegression(std::shared_ptr<const def::io::XmlModel> xml_model);
