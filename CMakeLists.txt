cmake_minimum_required(VERSION 3.11.0)

set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()
project(pt_icm)
# add files
file( GLOB SOURCES src/*.cc )
file( GLOB ESTIM_SOURCES estimator/*.cc )
file( GLOB TEST_SOURCES src/unittests/*.cc )
add_executable(${PROJECT_NAME} ${SOURCES})
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  target_compile_options(${PROJECT_NAME} BEFORE PUBLIC -O3 -fopenmp -funroll-loops)
endif()
if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  target_compile_options(${PROJECT_NAME} BEFORE PUBLIC -g -ggdb)
endif()
add_compile_options(-Wall -pedantic)
target_include_directories(${PROJECT_NAME} BEFORE PRIVATE include)
target_link_libraries(${PROJECT_NAME} PUBLIC boost_program_options)
add_library(replica STATIC src/replica.cc src/problem.cc src/utils.cc)
target_include_directories(replica BEFORE PRIVATE include)
add_compile_definitions(BOOST)
target_compile_definitions(${PROJECT_NAME} PUBLIC WHOLE)

add_executable(cae src/cae_estimator.cc)
target_link_libraries(cae replica)
target_include_directories(cae BEFORE PRIVATE include)
target_compile_options(cae BEFORE PUBLIC -O3 -Wall -pedantic)

# Build file for cae standalone

# Add source files
