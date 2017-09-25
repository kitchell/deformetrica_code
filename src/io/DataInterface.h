#ifndef DEFORMETRICA_DATAINTERFACE_H
#define DEFORMETRICA_DATAINTERFACE_H

#include <cstddef>

namespace def {
namespace interface {

template<typename ScalarType>
struct mat {
  ScalarType *ptr;
  std::size_t rows;
  std::size_t cols;

  mat(ScalarType **p, std::size_t rows, std::size_t cols)
      : ptr(*p), rows(rows), cols(cols)
  {
    ScalarType e0 = ptr[0];
    ScalarType e1 = ptr[1];
    ScalarType e2 = ptr[2];
    ScalarType e3 = ptr[3];
  }

  mat(const mat& x)
      : ptr(x.ptr), rows(x.rows), cols(x.cols)
  {
    ScalarType e = *ptr;
  }

};



}
}

#endif //DEFORMETRICA_DATAINTERFACE_H
