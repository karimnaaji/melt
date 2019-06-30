cmake_minimum_required(VERSION 2.8.12)

project(melt-gui)

if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 /O2")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O0 -Wall")
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

add_executable(${PROJECT_NAME} melt-gui.cpp deps/minitrace.c deps/glfw/deps/glad.c)

target_link_libraries(${PROJECT_NAME} glfw ${GLFW_LIBRARIES})
