CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(DebugMethods)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

FIND_PACKAGE(ITK)
IF (ITK_FOUND)
  INCLUDE(${USE_ITK_FILE})
ELSE (ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

FIND_PACKAGE(VTK)
IF (VTK_FOUND)
  INCLUDE(${USE_VTK_FILE})
ELSE (VTK_FOUND)
  MESSAGE(FATAL_ERROR "VTK not found. Please set VTK_DIR.")
ENDIF (VTK_FOUND)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/algo/deformable_objects
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/algo/deformations
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/algo/energy_minimization
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/algo/integrators
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/algo/kernels
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/algo/utilities
  ${CMAKE_CURRENT_SOURCE_DIR}/../../codes/src/io
)

#ADD_EXECUTABLE(surfcentering surfcentering.cxx surfio.cxx)
#ADD_EXECUTABLE(surfconv surfconv.cxx surfio.cxx)
#ADD_EXECUTABLE(surfdisp SurfaceViewer.cxx surfdisp.cxx surfio.cxx)
#ADD_EXECUTABLE(surfextract surfextract.cxx)
#ADD_EXECUTABLE(surfmerge surfmerge.cxx surfio.cxx)
#ADD_EXECUTABLE(surfsmooth surfsmooth.cxx surfio.cxx)
#ADD_EXECUTABLE(surfroi surfroi.cxx surfio.cxx)
#ADD_EXECUTABLE(templatedisp SurfaceViewer.cxx templatedisp.cxx surfio.cxx)
#ADD_EXECUTABLE(templateshoot SurfaceViewer.cxx templateshoot.cxx surfio.cxx ../app/myvtkPolyDataNormals.cxx)

ADD_EXECUTABLE(TestVectorizeImage TestVectorizeImage.cxx)

TARGET_LINK_LIBRARIES(TestVectorizeImage
 ${ITK_LIBRARIES}
 ${VTK_LIBRARIES}
)

ADD_EXECUTABLE(TestAnatomicalCoordinateSystem TestAnatomicalCoordinateSystem.cxx)

TARGET_LINK_LIBRARIES(TestAnatomicalCoordinateSystem
 ${ITK_LIBRARIES}
# ${VTK_LIBRARIES}
)

