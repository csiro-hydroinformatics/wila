cmake_minimum_required(VERSION 2.8)

SET (VERSION 0.5.5)

## section: Macro
MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)
  IF(MSVC)
    GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
    SET(PrecompiledBinary "${CMAKE_CURRENT_BINARY_DIR}/${PrecompiledBasename}.pch")
    SET(Sources ${${SourcesVar}})

    SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_OUTPUTS "${PrecompiledBinary}")
    SET_SOURCE_FILES_PROPERTIES(${Sources}
                                PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledBinary}\" /FI\"${PrecompiledBinary}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_DEPENDS "${PrecompiledBinary}")  
    # Add precompiled header to SourcesVar
    LIST(APPEND ${SourcesVar} ${PrecompiledSource})
  ENDIF(MSVC)
ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)

## start setting
SET (this_target wila)
PROJECT(${this_target})

## section: source files
# Add your source files here (one file per line), please SORT in alphabetical order for future maintenance
SET (${this_target}_SOURCE_FILES
    src/sce.cpp
)

## section: header files
# Add your header files here(one file per line), please SORT in alphabetical order for future maintenance!
SET(${this_target}_HEADER_FILES
    include/wila/core.hpp
    include/wila/evaluations.hpp
    include/wila/random.hpp
    include/wila/sce.hpp
    include/wila/sce.h
    include/wila/utils.hpp
)

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

find_package(TBB REQUIRED)
#find_path(TBB_INCLUDE_DIRS concurrent_vector.h PATH_SUFFIXES tbb)
#message ("The path to the TBB headers is ${TBB_INCLUDE_DIRS}")
#find_library (TBB_LIBRARIES NAMES tbb)
message ("The path to the Thread Building Blocks (tbb) library is ${TBB_LIBRARIES}")

find_path(THREADPOOL_INCLUDE_DIRS boost/threadpool.hpp PATHS threadpool ../threadpool)

INCLUDE_DIRECTORIES(
    ./include
    ./include/wila
    ../include
    ${Boost_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${THREADPOOL_INCLUDE_DIRS}
)

## section: precompiled header	
SET_SOURCE_FILES_PROPERTIES(${this_target}_HEADER_FILES
                            PROPERTIES HEADER_FILE_ONLY TRUE)
LIST(APPEND ${this_target}_SOURCE_FILES ${${this_target}_HEADER_FILES})

# JM needed to use this as of 2015-02, after an upgrade to latest Debian setup. The Add_definitions macro was not working, somehow.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

## section: add target
ADD_LIBRARY(${this_target} SHARED ${${this_target}_SOURCE_FILES} )
ADD_EXECUTABLE(wila_tests tests/common.cpp tests/main.cpp)
TARGET_LINK_LIBRARIES(wila_tests 
    pthread
    ${this_target}
    ${Boost_LIBRARIES}
    ${TBB_LIBRARIES}
)

# Files to be installed - used by packaging stuff below
INSTALL(FILES ${${this_target}_HEADER_FILES} DESTINATION include/wila)
INSTALL(TARGETS ${this_target} 
  LIBRARY DESTINATION lib
)

# Packages - currently just an RPM
set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "RPM;TGZ")
set(CPACK_PACKAGE_NAME "wila")
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_RELEASE ${CPACK_RPM_PACKAGE_RELEASE})
set(CPACK_PACKAGE_CONTACT "J-M <jean-michel.perraud@csiro.au>")
set(CPACK_PACKAGE_VENDOR "https://github.com/jmp75/wila")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
include(CPack)
