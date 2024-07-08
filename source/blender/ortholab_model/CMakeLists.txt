set(ORTHOLAB_LIB_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/extern/ortholab_lib)


set(INC
)

set(INC_SYS
  ${BOOST_INCLUDE_DIR}
  ${EIGEN3_INCLUDE_DIRS}
  ${ORTHOLAB_LIB_INCLUDE_DIRS}
)

set(LIB
)

set(SRC
  base_helpers.cpp
  extrusion.cpp
  file_services.cpp
  kvstore.cpp
  meshing.cpp
  pad.cpp
  object.cpp 
  project.cpp
  save.cpp
  scene_object.cpp
  smooth_mesh.cpp
  teeth.cpp
  tooth.cpp
  vedis.cpp
  visurep_base.cpp
  wire.cpp
  function.cpp
)

blender_add_lib(ortholab_model "${SRC}" "${INC}" "${INC_SYS}" "${LIB}")

target_link_libraries(ortholab_model PRIVATE ${PLATFORM_LINKLIBS})

set_target_properties(ortholab_model
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/${_LIB_SUB_FOLDER}"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/${_LIB_SUB_FOLDER}"
)
