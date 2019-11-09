file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories(..)
set(CMAKE_CXX_FLAGS "-g -O0 -std=c++14")

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