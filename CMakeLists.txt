cmake_minimum_required(VERSION 2.8.12)

project(melt)

if(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 /O2")
else()
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O0 -Wall -pedantic")
endif()
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "Build the GLFW example programs")
set(GLFW_BUILD_TESTS OFF CACHE BOOL "Build the GLFW test programs")
set(GLFW_BUILD_DOCS OFF CACHE BOOL "Build the GLFW documentation")
set(GLFW_INSTALL OFF CACHE BOOL "Generate installation target")

add_subdirectory(deps/glfw)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(deps)
include_directories(deps/glfw/deps)
include_directories(deps/glfw/include)
include_directories(deps/glm)

add_definitions(-DMTR_ENABLED)

add_executable(melt-gui melt-gui.cpp deps/minitrace.c deps/glfw/deps/glad.c)
target_link_libraries(melt-gui glfw ${GLFW_LIBRARIES})

add_executable(compile-check melt.cpp)

add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
