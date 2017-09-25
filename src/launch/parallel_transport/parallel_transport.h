#pragma once

namespace def {
namespace io {
class XmlModel;
}
};

#include <memory>


///Launcher for parallel transport computations. The standard use is the transport of a matching along
///a regression. The xml then contains paths to the momenta and control points of both of them
///and the number of time points is for the regression. the matching is assumed to have 10 time
///points.
///In both cases, the template object should be the origin point of the regression.
template<unsigned int Dimension>
void parallel_transport(std::shared_ptr<const def::io::XmlModel> xml_model);