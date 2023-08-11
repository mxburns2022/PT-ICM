cmake_minimum_required(VERSION 3.11.0)

set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
project(PT_ICM)
# add files
file( GLOB SOURCES src/*.cc )
file( GLOB ESTIM_SOURCES estimator/*.cc )
file( GLOB TEST_SOURCES src/unittests/*.cc )
set ( execname pt_icm )
add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} BEFORE PRIVATE include)
target_link_libraries(${PROJECT_NAME} PUBLIC boost_program_options)

# Build file for cae standalone
add_library(replica STATIC src/replica.cc src/problem.cc src/utils.cc)
target_include_directories(replica BEFORE PRIVATE include)
add_compile_definitions(BOOST)
add_compile_definitions(WHOLE)
# Add source files
