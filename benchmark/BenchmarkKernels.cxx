#include <iostream>
#include "benchmark/benchmark_api.h"
#include "src/support/kernels/ExactKernel.h"
#include "src/support/kernels/P3MKernel.h"
#include "src/support/kernels/Compact.h"

#ifdef USE_CUDA
#include "src/support/kernels/CUDAExactKernel.h"
#endif

#include <time.h>
#include <map>

#include "LinearAlgebra.h"

using namespace def::algebra;

/*
#ifdef USE_DOUBLE_PRECISION
typedef double ScalarType;
#else
typedef float ScalarType;
#endif
*/

enum {
  RUN_EXACT,
#ifdef USE_CUDA
  RUN_CUDA,
#endif
  RUN_P3M,
  RUN_COMPACT
};

namespace tear_up {

static VectorType generate_random_vector(const size_t size) {
  static std::map<size_t, VectorType> rand_vector;

  if (rand_vector.find(size) != rand_vector.end())
    return rand_vector[size];

  srand(time(0));

  VectorType vect(size);
  for (int i = 0; i < size; ++i) {
    vect[i] = (ScalarType) std::rand() / (ScalarType) (RAND_MAX);
  }

  rand_vector[size] = std::move(vect);

  return rand_vector[size];
}

static MatrixType generate_random_matrix(const size_t rows, const size_t cols) {
  typedef std::tuple<size_t, size_t> index;
  static std::map<index, MatrixType> rand_matrix;

  index i = std::make_tuple(rows, cols);

  if (rand_matrix.find(i) != rand_matrix.end())
    return rand_matrix[i];

  srand(time(0));

  MatrixType tmp(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      tmp(i, j) = (ScalarType) std::rand() / (ScalarType) (RAND_MAX);
    }
  }

  rand_matrix[i] = std::move(tmp);

  return rand_matrix[i];
}

static ExactKernel<ScalarType, 3> *get_implementation(const size_t type) {
  switch (type) {
    case RUN_EXACT:return new ExactKernel<ScalarType, 3>();
#ifdef USE_CUDA
      case RUN_CUDA:return new CUDAExactKernel<ScalarType, 3>();
#endif
    case RUN_P3M:return new P3MKernel<ScalarType, 3>();
    case RUN_COMPACT:return new Compact<ScalarType, 3>();
    default:assert(0);
  }
}

}


void Convolve_kernel(benchmark::State &state) {
  const ScalarType kernel_width = 3;
  auto X = tear_up::generate_random_matrix(state.range(0), 3);
  auto Y = tear_up::generate_random_matrix(state.range(0) * 5, 3);
  auto W = tear_up::generate_random_matrix(state.range(0) * 5, 6);

  auto kernel = tear_up::get_implementation(state.range(1));
  kernel->SetSources(Y);
  kernel->SetWeights(W);
  kernel->SetKernelWidth(kernel_width);

  // only for P3MKernel
  if (auto p3m = dynamic_cast<P3MKernel<ScalarType, 3> *>(kernel)) {
    // Determine the bounding box and other quantities for P3M kernels
    MatrixType DataDomain(3, 2, 0.0f);
    for (int k = 0; k < X.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (X(k, d) < DataDomain(d, 0) ? X(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (X(k, d) > DataDomain(d, 1) ? X(k, d) : DataDomain(d, 1));
      }
    }

    for (int k = 0; k < Y.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (Y(k, d) < DataDomain(d, 0) ? Y(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (Y(k, d) > DataDomain(d, 1) ? Y(k, d) : DataDomain(d, 1));
      }
    }

    p3m->SetDataDomain(DataDomain);
    p3m->SetWorkingSpacingRatio(0.2);
    p3m->SetPaddingFactor(3.0 * kernel_width);
  }

  while (state.KeepRunning()) {
    kernel->Convolve(X);
  }
}

void ConvolveGradient_kernel(benchmark::State &state) {
  const ScalarType kernel_width = 3;
  auto X = tear_up::generate_random_matrix(state.range(0), 3);
  auto Y = tear_up::generate_random_matrix(state.range(0) * 5, 3);
  auto W = tear_up::generate_random_matrix(state.range(0) * 5, 3);

  auto kernel = tear_up::get_implementation(state.range(1));
  kernel->SetSources(Y);
  kernel->SetWeights(W);
  kernel->SetKernelWidth(kernel_width);

  // only for P3MKernel
  if (auto p3m = dynamic_cast<P3MKernel<ScalarType, 3> *>(kernel)) {
    // Determine the bounding box and other quantities for P3M kernels
    MatrixType DataDomain(3, 2, 0.0f);
    for (int k = 0; k < X.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (X(k, d) < DataDomain(d, 0) ? X(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (X(k, d) > DataDomain(d, 1) ? X(k, d) : DataDomain(d, 1));
      }
    }

    for (int k = 0; k < Y.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (Y(k, d) < DataDomain(d, 0) ? Y(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (Y(k, d) > DataDomain(d, 1) ? Y(k, d) : DataDomain(d, 1));
      }
    }

    p3m->SetDataDomain(DataDomain);
    p3m->SetWorkingSpacingRatio(0.2);
    p3m->SetPaddingFactor(3.0 * kernel_width);
  }

  while (state.KeepRunning()) {
    kernel->ConvolveGradient(X, X);
  }
}

void ConvolveHessian_kernel(benchmark::State &state) {
  const ScalarType kernel_width = 3;
  auto X = tear_up::generate_random_matrix(state.range(0), 3);
  auto Y = tear_up::generate_random_matrix(state.range(0) * 5, 3);
  auto W = tear_up::generate_random_matrix(state.range(0) * 5, 3);

  auto kernel = tear_up::get_implementation(state.range(1));
  kernel->SetSources(Y);
  kernel->SetWeights(W);
  kernel->SetKernelWidth(kernel_width);

  // only for P3MKernel
  if (auto p3m = dynamic_cast<P3MKernel<ScalarType, 3> *>(kernel)) {
    // Determine the bounding box and other quantities for P3M kernels
    MatrixType DataDomain(3, 2, 0.0f);
    for (int k = 0; k < X.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (X(k, d) < DataDomain(d, 0) ? X(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (X(k, d) > DataDomain(d, 1) ? X(k, d) : DataDomain(d, 1));
      }
    }

    for (int k = 0; k < Y.rows(); k++) {
      for (int d = 0; d < 3; d++) {
        DataDomain(d, 0) = (Y(k, d) < DataDomain(d, 0) ? Y(k, d) : DataDomain(d, 0));
        DataDomain(d, 1) = (Y(k, d) > DataDomain(d, 1) ? Y(k, d) : DataDomain(d, 1));
      }
    }

    p3m->SetDataDomain(DataDomain);
    p3m->SetWorkingSpacingRatio(0.2);
    p3m->SetPaddingFactor(3.0 * kernel_width);
  }

  while (state.KeepRunning()) {
    kernel->ConvolveSpecialHessian(W);
  }
}


void Evaluate_Kernel_Vector(benchmark::State &state) {
  const ScalarType kernel_width = 3;
  auto X = tear_up::generate_random_vector(state.range(0));
  auto Y = X * M_PI;

  auto kernel = tear_up::get_implementation(state.range(1));
  kernel->SetKernelWidth(kernel_width);

  while (state.KeepRunning()) {
    kernel->EvaluateKernel(X,Y);
  }
}

void Evaluate_Kernel_Matrix(benchmark::State &state) {
  const ScalarType kernel_width = 3;
  auto X = tear_up::generate_random_matrix(state.range(0), 3);
  auto Y = X * M_PI;

  auto kernel = tear_up::get_implementation(state.range(1));
  kernel->SetKernelWidth(kernel_width);

  const unsigned int R = X.rows();
  const unsigned int C = X.cols();

  while (state.KeepRunning()) {
    for (int i = 0; i < R; i++) {
      for (int j = 0; j < C; j++) {
        kernel->EvaluateKernel(X, Y, i, j);
      }
    }
  }
}


// Declare ranges
#define BASIC_BENCHMARK_EVALUATE_KERNEL(x, y) \
    BENCHMARK(x)->ArgPair(1<<16, y)->ArgPair(1<<18, y)->ArgPair(1<<20, y)->ArgPair(1<<22, y)->Unit(benchmark::kMillisecond);
#define BASIC_BENCHMARK_TEST(x, y) \
    BENCHMARK(x)->ArgPair(1<<12, y)->ArgPair(1<<14, y)->ArgPair(1<<16, y)->Unit(benchmark::kMillisecond);
#define BASIC_BENCHMARK_TEST_SMALL(x, y) \
    BENCHMARK(x)->ArgPair(1<<8, y)->ArgPair(1<<10, y)->ArgPair(1<<12, y)->Unit(benchmark::kMillisecond);

// Run tests
BASIC_BENCHMARK_TEST_SMALL(ConvolveGradient_kernel, RUN_EXACT);
BASIC_BENCHMARK_TEST_SMALL(ConvolveGradient_kernel, RUN_COMPACT);
BASIC_BENCHMARK_TEST(ConvolveGradient_kernel, RUN_P3M);
#ifdef USE_CUDA
BASIC_BENCHMARK_TEST(ConvolveGradient_kernel, RUN_CUDA);
#endif


BASIC_BENCHMARK_TEST_SMALL(ConvolveGradient_kernel, RUN_EXACT);
BASIC_BENCHMARK_TEST_SMALL(ConvolveGradient_kernel, RUN_COMPACT);

BASIC_BENCHMARK_EVALUATE_KERNEL(Evaluate_Kernel_Matrix, RUN_EXACT);
BASIC_BENCHMARK_EVALUATE_KERNEL(Evaluate_Kernel_Matrix, RUN_COMPACT);

BASIC_BENCHMARK_EVALUATE_KERNEL(Evaluate_Kernel_Vector, RUN_EXACT);
BASIC_BENCHMARK_EVALUATE_KERNEL(Evaluate_Kernel_Vector, RUN_COMPACT);


BASIC_BENCHMARK_TEST_SMALL(Convolve_kernel, RUN_EXACT);
BASIC_BENCHMARK_TEST(Convolve_kernel, RUN_P3M);
#ifdef USE_CUDA
BASIC_BENCHMARK_TEST(Convolve_kernel, RUN_CUDA);
#endif

BASIC_BENCHMARK_TEST_SMALL(ConvolveHessian_kernel, RUN_EXACT);
BASIC_BENCHMARK_TEST(ConvolveHessian_kernel, RUN_P3M);
#ifdef USE_CUDA
BASIC_BENCHMARK_TEST(ConvolveHessian_kernel, RUN_CUDA);
#endif

BENCHMARK_MAIN();
