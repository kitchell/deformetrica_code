#ADD_EXECUTABLE(deformetrica-maps-ellipsoid MapsEllipsoid.cxx)
ADD_EXECUTABLE(deformetrica-maps-ellipsoid-with-source MapsEllipsoidWithSource.cxx)
ADD_EXECUTABLE(deformetrica-compute-deformation-distance ComputeDeformationDistance.cxx)
ADD_EXECUTABLE(deformetrica-compute-varifold-distance ComputeVarifoldDistance.cpp)
ADD_EXECUTABLE(deformetrica-compute-current-distance ComputeCurrentDistance.cpp)
ADD_EXECUTABLE(deformetrica-compute-image-distance ComputeImageDistance.cpp)
ADD_EXECUTABLE(deformetrica-simulate-longitudinal-data SimulateLongitudinalData.cpp)
ADD_EXECUTABLE(deformetrica-simulate-lda-data SimulateLdaAtlas.cpp)
ADD_EXECUTABLE(deformetrica-smooth-mesh SmoothMesh.cpp)
ADD_EXECUTABLE(deformetrica-decimate-mesh DecimateMesh.cpp)
ADD_EXECUTABLE(deformetrica-reorient-mesh-normals ReorientMeshNormals.cpp)
ADD_EXECUTABLE(deformetrica-inspect-serialization InspectSerialization.cpp)

#TARGET_LINK_LIBRARIES(deformetrica-maps-ellipsoid
#        ${ITK_LIBRARIES}
#        ${VTK_DEPENDANCIES}
#        )
#
#message ("vtk " ${VTK_DEPENDANCIES})

TARGET_LINK_LIBRARIES(deformetrica-maps-ellipsoid-with-source
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-compute-deformation-distance
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-compute-varifold-distance
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-compute-current-distance
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-compute-image-distance
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-simulate-longitudinal-data
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-simulate-lda-data
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-smooth-mesh
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-decimate-mesh
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-reorient-mesh-normals
    deformetrica-lib
    )

TARGET_LINK_LIBRARIES(deformetrica-inspect-serialization
    deformetrica-lib
    )


#if ( "${VTK_MAJOR_VERSION}" LESS 6 )
#
##    Add_Dependencies (deformetrica-maps-ellipsoid
##            ${ITK_LIBRARIES}
##            ${VTK_DEPENDANCIES})
#
#    Add_Dependencies (deformetrica-maps-ellipsoid-with-source
#            ${ITK_LIBRARIES}
#            ${VTK_DEPENDANCIES})
#
#endif ( "${VTK_MAJOR_VERSION}" LESS 6 )


#install (TARGETS deformetrica-maps-ellipsoid DESTINATION bin)
install(TARGETS deformetrica-maps-ellipsoid-with-source DESTINATION bin)
install(TARGETS deformetrica-compute-deformation-distance DESTINATION bin)
install(TARGETS deformetrica-compute-image-distance DESTINATION bin)
install(TARGETS deformetrica-compute-varifold-distance DESTINATION bin)
install(TARGETS deformetrica-compute-current-distance DESTINATION bin)
install(TARGETS deformetrica-simulate-longitudinal-data DESTINATION bin)
install(TARGETS deformetrica-simulate-lda-data DESTINATION bin)
install(TARGETS deformetrica-smooth-mesh DESTINATION bin)
install(TARGETS deformetrica-decimate-mesh DESTINATION bin)
install(TARGETS deformetrica-reorient-mesh-normals DESTINATION bin)
install(TARGETS deformetrica-inspect-serialization DESTINATION bin)


