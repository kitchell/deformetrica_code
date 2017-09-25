cmake_Minimum_Required(VERSION 3.2)

project(Deformetrica CXX)
Set(DEFORMETRICA_VERSION_MAJOR 3)
Set(DEFORMETRICA_VERSION_MINOR 0)
Set(DEFORMETRICA_VERSION_PATCH 0)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif ()

function(PRINT_ENV print)
    if (${print})
        get_cmake_property(_variableNames VARIABLES)
        message("==============================================")
        message("       LIST OF CURRENT ENVIRONMENT VARIABLES")
        foreach (_variableName ${_variableNames})
            message(STATUS "${_variableName}=${${_variableName}}")
        endforeach ()
        message("==============================================")
    endif ()
endfunction()

option(BUILD_TESTS "Build tests ?" ON)
option(BUILD_BENCHMARKS "Build benchmarks ?" ON)
option(BUILD_UTILS "Build Utils ?" ON)
Option(CMAKE_DEBUG "Print the current environment variables ?" OFF)
option(FORCE_INSTALL "Ignore the version number of the dependencies ?" OFF)

##Deformetrica 3.0.0 us eonly double and deep precision
#Option(USE_DOUBLE_PRECISION "Deformetrica will make computations in double precision ?" ON)
#Option(USE_FAST_MATH "Deformetrica will use approximation maths function ?" OFF)

##This is a C++11 project!
include(CheckCXXCompilerFlag)
Set(ENABLE_CXXFLAGS_TO_CHECK
    -std=gnu++1z
    -std=c++1z
    -std=gnu++1y
    -std=c++1y
    -std=gnu++11
    -std=c++11
    -std=gnu++0x
    -std=c++0x)

Set(CXX_COMPILER_ENABLED FALSE)

foreach (flag ${ENABLE_CXXFLAGS_TO_CHECK})
    string(REPLACE "-std=" "_" flag_var ${flag})
    string(REPLACE "+" "x" flag_var ${flag_var})
    CHECK_CXX_COMPILER_FLAG("${flag}" COMPILER_HAS_CXX_FLAG${flag_var})
    if (COMPILER_HAS_CXX_FLAG${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        Set(CXX_COMPILER_ENABLED TRUE)
        break()
    endif ()
endforeach ()

if (NOT CXX_COMPILER_ENABLED)
    PRINT_ENV(${CMAKE_DEBUG})
    Mesage(FATAL_ERROR "The current CPP compiler does not support C++11: Please provide a good one.")
endif ()

#do not disable enable_testing(): the lib/cuda_convolution will not compile!
enable_testing()

if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif ()

# RELEASE version by default :
if (NOT CMAKE_BUILD_TYPE)
    Set(CMAKE_BUILD_TYPE "RELEASE")
endif ()

# Extra flags for the C++ compiler :
if (CMAKE_BUILD_TYPE MATCHES "RELEASE")
    Message("Deformetrica will compile in RELEASE mode")
    Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DARMA_NO_DEBUG -DNDEBUG -msse4.2")
endif ()

# Extra flags for the C++ compiler if version is Gprof :
if (CMAKE_BUILD_TYPE MATCHES "GPROF")
    Message("Deformetrica will compile in GPROF mode")
    Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
    List(APPEND CUDA_NVCC_FLAGS "-g")
endif ()

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
    Message("Deformetrica will compile in DEBUG mode")
    Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    List(APPEND CUDA_NVCC_FLAGS "-g")
endif ()

## Find packages:

# -Linear Algebra: ARMADILLO
find_Package(Armadillo)
if (NOT ARMADILLO_FOUND)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "ARMADILLO library is not found")
endif ()

if (NOT FORCE_INSTALL AND (NOT ARMADILLO_VERSION_MAJOR OR "${ARMADILLO_VERSION_MAJOR}" LESS 6))
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "ARMADILLO version must be great or equal to version 6.0")
endif ()

# -Image Processing, Images I/O: VTK
Find_Package(VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkIOLegacy vtkFiltersCore vtkFiltersGeneral NO_MODULE)
if (NOT VTK_FOUND)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "VTK library is not found")
endif ()

if (NOT VTK_MAJOR_VERSION OR ${VTK_MAJOR_VERSION} LESS 6)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "VTK version must be great or equal to version 6.0")
endif ()

# -Image Processing, Segmentation, Registration: ITK
find_Package(ITK)
if (NOT ITK_FOUND)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "ITK library is not found")
endif ()

if (NOT ITK_VERSION_MAJOR OR ${ITK_VERSION_MAJOR} LESS 4)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "ITK version must be great or equal to version 4.0")
endif ()

# -General purpose programming: Boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost COMPONENTS system filesystem serialization)
if (NOT Boost_FOUND)
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "Boost library is not found")
endif ()

if ((NOT Boost_MAJOR_VERSION) OR (NOT "${Boost_MAJOR_VERSION}" STREQUAL "1") OR (${Boost_MINOR_VERSION} LESS 45))
    PRINT_ENV(${CMAKE_DEBUG})
    Message(FATAL_ERROR "Boost version must be great or equal to version 1.45")
endif ()

# -GPU Programming: CUDA
Find_Package(CUDA)

if (CUDA_FOUND)
    Option(USE_CUDA "Use cuda implementation of convolutions ?" ON)

    if (USE_CUDA)
        Message("Using CUDA ${CUDA_VERSION}")
        Include_Directories(${CUDA_INCLUDE_DIRS})

        Message("------------------------------------------------------------------------------------------------------------")
        Message("CUDA: Please check the CUDA_NVCC_FLAGS within CMakeLists.txt in order to be update with your current nvidia card")
        Message("CUDA: Compiling with an old architecture parameters may produce wrong results!")
        Message("------------------------------------------------------------------------------------------------------------")
        List(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30")
        List(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_35,code=sm_35")
        List(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_52,code=sm_52")
        List(APPEND CUDA_NVCC_FLAGS "--compiler-options=-fPIC")

        if (USE_FAST_MATH)
            List(APPEND CUDA_NVCC_FLAGS "--use_fast_math")
        endif ()
        Set(CUDA_VERBOSE_BUILD ON)
        Set(CUDA_PROPAGATE_HOST_FLAGS OFF)

        Message("Deformetrica will be compiled with Cuda")

        Cuda_Add_Library(
            CudaConvolution STATIC
            ${CMAKE_CURRENT_SOURCE_DIR}/lib/cuda_convolutions/RadialFunction.h
            ${CMAKE_CURRENT_SOURCE_DIR}/lib/cuda_convolutions/ScalarRadialKernel.h
            ${CMAKE_CURRENT_SOURCE_DIR}/lib/cuda_convolutions/GpuConv1D.cu
        )
        Set(CUDA_LIBS ${CUDA_LIBS} CudaConvolution)

    else ()
        Message("CUDA is found but it will not be used")
    endif (USE_CUDA)
else ()
    Message("CUDA was not found.")
    SET(USE_CUDA OFF)
endif ()

set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}

    ${CMAKE_CURRENT_SOURCE_DIR}/src

    ${CMAKE_CURRENT_SOURCE_DIR}/src/launch
    ${CMAKE_CURRENT_SOURCE_DIR}/src/launch/atlas
    ${CMAKE_CURRENT_SOURCE_DIR}/src/launch/matching
    ${CMAKE_CURRENT_SOURCE_DIR}/src/launch/regression


    ${CMAKE_CURRENT_SOURCE_DIR}/src/core
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/estimators
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/estimators_tools
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/estimators_tools/samplers
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/models
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/models/atlases
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/models/longitudinal_models
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/models/regressions
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/model_tools
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/model_tools/deformations
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/observations
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/observations/data_sets
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/observations/deformable_objects
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/observations/deformable_objects/geometries
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/observations/deformable_objects/geometries/images
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/observations/deformable_objects/geometries/landmarks

    ${CMAKE_CURRENT_SOURCE_DIR}/src/io

    ${CMAKE_CURRENT_SOURCE_DIR}/src/support
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support/fast_math
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support/linear_algebra
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support/linear_algebra/armadillo
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support/probability_distributions
    ${CMAKE_CURRENT_SOURCE_DIR}/src/support/utilities

    ${PROJECT_BINARY_DIR}
    )

Include_Directories(${SOURCE_FILES})
## Print the version of libraries and inlcude headers
Message("Using ARMADILLO Version ${ARMADILLO_VERSION_STRING} '${ARMADILLO_VERSION_NAME}'")
Include_Directories(${ARMADILLO_INCLUDE_DIRS})

Message("Using ITK Version ${ITK_VERSION}")
Include(${ITK_USE_FILE})

Message("Using VTK Version ${VTK_VERSION}")
Include(${VTK_USE_FILE})

Message("Using BOOST Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
include_directories(${Boost_INCLUDE_DIRS})

Add_Library(deformetrica-lib
    STATIC
    src/launch/deformetrica.cxx
    src/launch/atlas/estimate_atlas.cxx
    src/launch/longitudinal_atlas/estimate_longitudinal_atlas.cxx
    src/launch/longitudinal_registration/estimate_longitudinal_registration.cxx
    src/launch/parallel_transport/parallel_transport.cxx
    src/launch/regression/estimate_geodesic_regression.cxx
    src/launch/shooting/deform.cxx
    src/io/SparseDiffeoParametersXMLFile.cxx
    src/io/DeformableObjectParametersXMLFile.cxx
    src/io/SparseDiffeoParameters.cxx
    src/io/DeformableObjectParameters.cxx
    src/io/XmlConfigurationConverter.cxx
    src/io/DeformableObjectReader.cxx
    src/io/DeformationFieldIO.cxx
    src/io/MatrixDLM.cxx
    src/io/SparseDiffeoWriter.cxx
    src/support/utilities/SimpleTimer.cxx
    src/support/utilities/myvtkPolyDataNormals.cxx
    src/support/utilities/SerializeDeformationState.h
    src/support/utilities/GeneralSettings.cxx
    src/support/utilities/SerializeDeformationState.cxx
    src/support/utilities/AnatomicalCoordinateSystem.cxx
    src/support/utilities/GridFunctions.cxx
    src/support/probability_distributions/AbstractNormalDistribution.cxx
    src/support/probability_distributions/AbstractProbabilityDistribution.cxx
    src/support/probability_distributions/ConditionedNormalDistribution.cxx
    src/support/probability_distributions/DirichletDistribution.cxx
    src/support/probability_distributions/DisplacementFieldNormalDistribution.cxx
    src/support/probability_distributions/InverseWishartDistribution.cxx
    src/support/probability_distributions/MultiScalarInverseWishartDistribution.cxx
    src/support/probability_distributions/MultiScalarNormalDistribution.cxx
    src/support/probability_distributions/NormalDistribution.cxx
    src/support/probability_distributions/UniformDistribution.cxx
    src/support/probability_distributions/AutomaticRelevanceDeterminationDistribution.cxx
    src/support/kernels/AbstractKernel.cxx
    src/support/kernels/CUDAExactKernel.cxx
    src/support/kernels/ExactKernel.cxx
    src/support/kernels/P3MKernel.cxx
    src/support/kernels/Compact.cxx
    src/support/kernels/KernelFactory.cxx
    src/support/linear_algebra/ArmadilloMatrixWrapper.cxx
    src/support/linear_algebra/ArmadilloVectorWrapper.cxx
    src/support/linear_algebra/LinearVariableMapWrapper.cxx
    src/support/linear_algebra/LinearVariablesMapWrapper.cxx
    src/support/linear_algebra/LinearVariableWrapper.cxx
    src/support/linear_algebra/MatrixListWrapper.cxx
    src/core/models/atlases/AbstractAtlas.cxx
    src/core/models/atlases/BayesianAtlas.cxx
    src/core/models/atlases/BayesianAtlasMixture.cxx
    src/core/models/atlases/DeterministicAtlas.cxx
    src/core/models/atlases/LdaAtlas.cxx
    src/core/models/longitudinal_models/LongitudinalAtlas.cxx
    src/core/models/longitudinal_models/LongitudinalRegistration.cxx
    src/core/models/regressions/Regression.cxx
    src/core/models/AbstractStatisticalModel.cxx
    src/core/observations/data_sets/CrossSectionalDataSet.cxx
    src/core/observations/data_sets/LongitudinalDataSet.cxx
    src/core/observations/data_sets/TimeSeriesDataSet.cxx
    src/core/observations/deformable_objects/DeformableMultiObject.cxx
    src/core/observations/deformable_objects/geometries/AbstractGeometry.cxx
    src/core/observations/deformable_objects/geometries/images/EQLAImage.cxx
    src/core/observations/deformable_objects/geometries/images/LCCImage.cxx
    src/core/observations/deformable_objects/geometries/images/LinearInterpImage.cxx
    src/core/observations/deformable_objects/geometries/images/MutualInformationImage.cxx
    src/core/observations/deformable_objects/geometries/images/ParametricImage.cxx
    src/core/observations/deformable_objects/geometries/images/SSDImage.cxx
    src/core/observations/deformable_objects/geometries/landmarks/Landmark.cxx
    src/core/observations/deformable_objects/geometries/landmarks/PointCloud.cxx
    src/core/observations/deformable_objects/geometries/landmarks/NonOrientedPolyLine.cxx
    src/core/observations/deformable_objects/geometries/landmarks/NonOrientedSurfaceMesh.cxx
    src/core/observations/deformable_objects/geometries/landmarks/OrientedPolyLine.cxx
    src/core/observations/deformable_objects/geometries/landmarks/OrientedSurfaceMesh.cxx
    src/core/observations/deformable_objects/geometries/landmarks/OrientedVolumeMesh.cxx
    src/core/estimators/AbstractEstimator.cxx
    src/core/estimators/FastGradientAscent.cxx
    src/core/estimators/GradientAscent.cxx
    src/core/estimators/McmcSaem.cxx
    src/core/estimators/PowellsMethod.cxx
    src/core/estimators_tools/samplers/AbstractSampler.cxx
    src/core/estimators_tools/samplers/AmalaSampler.cxx
    src/core/estimators_tools/samplers/MalaSampler.cxx
    src/core/estimators_tools/samplers/SrwMhwgSampler.cxx
    src/core/model_tools/deformations/AbstractDeformations.cxx
    src/core/model_tools/deformations/AdjointEquationsIntegrator.cxx
    src/core/model_tools/deformations/Diffeos.cxx
    src/io/XmlDataSet.hpp
    src/support/utilities/Utils.hpp
    )


ADD_EXECUTABLE(deformetrica
  src/launch/deformetrica.cxx)

Add_Executable(CreateDeformationField
    src/launch/deformation/CreateDeformationField.cxx)

#
# Target link :
#
List(APPEND VTK_DEPENDANCIES ${VTK_LIBRARIES})

Target_Link_Libraries(deformetrica-lib
    ${ARMADILLO_LIBRARIES}
    ${ITK_LIBRARIES}
    ${VTK_DEPENDANCIES}
    ${Boost_LIBRARIES}
    ${CUDA_LIBS})

Target_Link_Libraries(deformetrica deformetrica-lib)
Target_Link_Libraries(CreateDeformationField deformetrica-lib)

#Add Tests
if (BUILD_TESTS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest)
    add_subdirectory(tests)
endif ()

#Adding Benchmarks
if (BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/benchmark)
    add_subdirectory(benchmark)
endif ()

if (BUILD_UTILS)
    add_subdirectory(utilities/c++)
endif ()

# Configure a header file to pass some of the CMake Settings to the source code :
Configure_File(
    "${PROJECT_SOURCE_DIR}/DeformetricaConfig.h.in"
    "${PROJECT_BINARY_DIR}/DeformetricaConfig.h")

## INSTALLATION
install(TARGETS deformetrica DESTINATION bin)

PRINT_ENV(${CMAKE_DEBUG})
