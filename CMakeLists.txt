cmake_minimum_required(VERSION 2.8.12)

project(melt)

if(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 /O2")
else()
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O0 -Wall -pedantic")
endif()

include_directories(${PROJECT_SOURCE_DIR})
include_directories(deps)
include_directories(deps/glm)

set(MODEL_FILES
  tests/models/suzanne.obj
  tests/models/column.obj
  tests/models/cube.obj
  tests/models/sphere.obj
  tests/models/teapot.obj
  tests/models/bunny.obj
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

add_executable(compile-check melt.cpp)

add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
