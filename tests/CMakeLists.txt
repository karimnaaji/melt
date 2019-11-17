file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories(..)
set(CMAKE_CXX_FLAGS "-g -O0 -std=c++14")

set(MODEL_FILES
  models/suzanne.obj
  models/column.obj
  models/cube.obj
  models/sphere.obj
  models/teapot.obj
  models/bunny.obj
)

foreach(model ${MODEL_FILES})
  get_filename_component(model_name "${model}" NAME_WE)
  get_filename_component(model_ext "${model}" EXT)
  string(REGEX REPLACE "\\." "" model_ext "${model_ext}")
  set(model_name "${model_name}_${model_ext}")
  set(MODEL_OUTPUT "${PROJECT_SOURCE_DIR}/generated/${model_name}.h")
  add_custom_command(
    OUTPUT "${MODEL_OUTPUT}"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    MAIN_DEPENDENCY "${model}"
    COMMAND ${CMAKE_COMMAND} -Din="${model}" -Dout="${MODEL_OUTPUT}" -Dname="${model_name}" -P ${PROJECT_SOURCE_DIR}/models.cmake)
  list(APPEND MODEL_OUTPUT_FILES "${MODEL_OUTPUT}")
endforeach()

add_definitions(-DMTR_ENABLED)
add_custom_target(model-headers DEPENDS ${MODEL_OUTPUT_FILES})

macro(add_resources TARGET RESOURCE_DIR DEST_DIR)
  add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${RESOURCE_DIR} ${DEST_DIR})
endmacro(add_resources)

foreach(src_file ${TEST_SOURCES})
  message(STATUS ${src_file})
  string(REPLACE ".cpp" "" test_case ${src_file})
  string(REGEX MATCH "([^/]*)$" test_name ${test_case})

  set(EXECUTABLE_NAME "${test_name}.out")
  add_executable(${EXECUTABLE_NAME} ${src_file})
  add_resources(${EXECUTABLE_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/models models)
endforeach()
