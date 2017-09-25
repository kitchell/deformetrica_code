find_package(Threads REQUIRED)

macro(compile_benchmark_test name)
  add_executable(${name} "${name}.cxx")
  target_link_libraries(${name} benchmark deformetrica-lib ${CMAKE_THREAD_LIBS_INIT})
endmacro(compile_benchmark_test)

compile_benchmark_test(BenchmarkKernels)

