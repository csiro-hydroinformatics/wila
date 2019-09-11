cmake_minimum_required(VERSION 3.0)

## start setting
SET (this_target mhcpptest)
PROJECT(${this_target})

# Use the FindBoost macro to find boost.
# Note that if boost is installed in a non-standard location you may
# need to help it by setting some variables. For example, on
# CENTOS 5 or 6:
#     cmake -DBOOST_INCLUDEDIR=/usr/include/boost148 \
#           -DBOOST_LIBRARYDIR=/usr/lib64/boost148 .
find_package(Boost 1.48.0 COMPONENTS system date_time thread REQUIRED)

IF(DEFINED ENV{INTEL_CC_HOME})
  message("Found env var INTEL_CC_HOME, using TBB_ROOT_DIR=$ENV{INTEL_CC_HOME}/tbb as a hint for TBB package discovery")
  set(TBB_ROOT_DIR $ENV{INTEL_CC_HOME}/tbb)
ENDIF()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/..")

if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake_modules")
    INCLUDE(CodeCoverage)
endif()

find_package(TBB REQUIRED)
find_path(THREADPOOL_INCLUDE_DIRS boost/threadpool.hpp PATHS threadpool ../threadpool ~/src/csiro/stash/per202/swift/Externals/threadpool)

INCLUDE_DIRECTORIES(
    ../include
    ${Boost_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${THREADPOOL_INCLUDE_DIRS}
)


## section: source files
# Add your source files here (one file per line), please SORT in alphabetical order for future maintenance
SET (${this_target}_SOURCE_FILES
    common.cpp
    main.cpp
)

## section: header files
# Add your header files here(one file per line), please SORT in alphabetical order for future maintenance!
SET(${this_target}_HEADER_FILES
    )

add_definitions(-DBOOST_THREAD_USE_LIB)

## section: add definitions
#   add prefix -D. example> -DSHP
#  - DO NOT add  the following definitions(already defined in ${OSP_DEFINITIONS}:
#   -DSHP, -DWIN32, -D_WINDOWS, -D_DEBUG, -D_USRDLL, -D_CRT_SECURE_NO_DEPRECATE
ADD_DEFINITIONS(
    # -DCMAKE_CXX_FLAGS:STRING=' -std=c++11 '    
    )

# JM needed to use this as of 2015-02, after an upgrade to latest Debian setup. The Add_definitions macro was not working, somehow.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

## section: add target

ADD_EXECUTABLE(${this_target} ${${this_target}_SOURCE_FILES} )



## section: set link libraries
TARGET_LINK_LIBRARIES( ${this_target}
	${Boost_LIBRARIES}
	${TBB_LIBRARIES}
)

