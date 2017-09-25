
#Define variable to use for tests
add_definitions(-DFUNCTIONAL_TESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/functional_tests")
add_definitions(-DUNIT_TESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/unit_tests")

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB basic_test_files unit_tests/io/TestReadConfiguration.cxx unit_tests/io/TestReadConfiguration.h)
file(GLOB basic_test_files unit_tests/io/TestReadParametersXML.cxx unit_tests/io/TestReadParametersXML.h ${basic_test_files})
file(GLOB basic_test_files unit_tests/parallel-transport/TestParallelTransport.cxx unit_tests/parallel-transport/TestParallelTransport.h ${basic_test_files})
file(GLOB basic_test_files unit_tests/geometries/TestVolumeGradient.cxx unit_tests/geometries/TestVolumeGradient.h ${basic_test_files})
file(GLOB basic_test_files unit_tests/serialize/TestSerialization.cxx unit_tests/serialize/TestSerialization.h ${basic_test_files})
file(GLOB basic_test_files unit_tests/linear_algebra/TestBoostWrappers.cxx unit_tests/linear_algebra/TestBoostWrappers.h ${basic_test_files})

file(GLOB basic_test_files unit_tests/kernels/AbstractTestKernelPrecision.cxx unit_tests/kernels/AbstractTestKernelPrecision.h ${basic_test_files})
file(GLOB basic_test_files unit_tests/kernels/TestKernelPrecisionP3M.cxx unit_tests/kernels/TestKernelPrecisionP3M.h ${basic_test_files})
if(USE_CUDA)
    file(GLOB cuda_test_files unit_tests/kernels/TestKernelPrecisionCUDA.cxx unit_tests/kernels/TestKernelPrecisionCUDA.h)
endif()

file(GLOB functional_test_files functional_tests/TestFunctional.cxx functional_tests/TestFunctional.h)

add_executable(RunUnitTests RunTests.cxx ${basic_test_files} ${cuda_test_files})

target_link_libraries(RunUnitTests
    gtest
    gtest_main
    deformetrica-lib
    )

add_test(run-unit-test RunUnitTests)

if (UNIX OR APPLE)
    add_executable(RunFunctionalTests RunTests.cxx ${functional_test_files})
    target_link_libraries(RunFunctionalTests
        gtest
        gtest_main
        deformetrica-lib
        )
endif ()
