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

#include <memory>

template<unsigned int Dimension>
void estimateLongitudinalRegistration(std::shared_ptr<const def::io::XmlModel> xml_model);