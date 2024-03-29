cmake_minimum_required(VERSION 3.10)

IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Coverage Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Coverage Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

set(TARGET_LIB_NAME wila)
SET (this_target ${TARGET_LIB_NAME})

set(TARGET_MAJOR 0)
set(TARGET_MINOR 7)
set(TARGET_RELEASE 2)
# set(TARGET_NAME libmoirai${TARGET_MAJOR}${TARGET_MINOR})
set(TARGET_VERSION ${TARGET_MAJOR}.${TARGET_MINOR}.${TARGET_RELEASE})
SET (VERSION ${TARGET_VERSION})

PROJECT(${this_target} VERSION ${TARGET_VERSION} DESCRIPTION "A simple metaheuristics framework to couple optimization and simulation models")

include(GNUInstallDirs)

## section: header files
# Add your header files here(one file per line), please SORT in alphabetical order for future maintenance!
SET(${this_target}_HEADER_FILES
    include/wila/common.h
    include/wila/constructs.hpp
    include/wila/core.h
    include/wila/core.hpp
    include/wila/evaluations.hpp
    include/wila/interop_c_cpp.hpp
    include/wila/interop_c_structs.h
    include/wila/interop_rcpp.hpp
    include/wila/logging.hpp
    include/wila/multithreading.hpp
    include/wila/random.hpp
    include/wila/sce.h
    include/wila/sce.hpp
    include/wila/urs.hpp
    include/wila/utils.hpp
    include/wila/wrappers.hpp
)

## section: source files
# Add your source files here (one file per line), please SORT in alphabetical order for future maintenance
# 2020-10 this is for legacy reasons; added perhaps for needs by the australian BoM
SET (${this_target}_SOURCE_FILES
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

message ("Adding to CMAKE_MODULE_PATH: ${PROJECT_SOURCE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake_modules")
    INCLUDE(CodeCoverage)
endif()

IF(DEFINED ENV{WILA_CMAKE_VERBOSE})
    message ("before testing: TBB_FOUND=${TBB_FOUND}")
ENDIF()

if(NOT TBB_FOUND) # Trying to fix appveyor build, which seems not to find FindTBB.cmake. 
find_package(TBB REQUIRED)
ENDIF()
#find_path(TBB_INCLUDE_DIRS concurrent_vector.h PATH_SUFFIXES tbb)
#message ("The path to the TBB headers is ${TBB_INCLUDE_DIRS}")
#find_library (TBB_LIBRARIES NAMES tbb)

find_path (THREADPOOL_INCLUDE_DIRS boost/threadpool.hpp)
find_path (CINTEROP_INCLUDE_DIRS cinterop/common_c_interop.h)
find_path (CATCH_INCLUDE_DIRS NAMES catch/catch.hpp )

IF(DEFINED ENV{WILA_CMAKE_VERBOSE})
	message ("The path to the boost headers is ${Boost_INCLUDE_DIRS}")
	message ("The path to the boost libraries is ${Boost_LIBRARIES}")
	message ("TBB_FOUND=${TBB_FOUND}")
	# message ("The path to the TBB headers is ${TBB_INCLUDE_DIRS}")
	# message ("The path to the TBB libraries is ${TBB_LIBRARIES}")
	message ("The path to the threadpool headers is ${THREADPOOL_INCLUDE_DIRS}")
	message ("The path to the cinterop headers is ${CINTEROP_INCLUDE_DIRS}")
	message ("The path to the catch headers is ${CATCH_INCLUDE_DIRS}")
ENDIF()


INCLUDE_DIRECTORIES(
    ./include
    ./include/wila
    ${Boost_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${THREADPOOL_INCLUDE_DIRS}
    ${CATCH_INCLUDE_DIRS}
    ${CINTEROP_INCLUDE_DIRS}
)

## section: precompiled header
SET_SOURCE_FILES_PROPERTIES(${this_target}_HEADER_FILES
                            PROPERTIES HEADER_FILE_ONLY TRUE)
LIST(APPEND ${this_target}_SOURCE_FILES ${${this_target}_HEADER_FILES})

# JM needed to use this as of 2015-02, after an upgrade to latest Debian setup. The Add_definitions macro was not working, somehow.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

## section: add target
# 2020-10 AS_SHARED_LIB is deprecated but kept for legacy reasons; may address needs of the australian BoM.
IF (AS_SHARED_LIB)
    ADD_LIBRARY(${this_target} SHARED ${${this_target}_SOURCE_FILES} )
    set_target_properties(${this_target} PROPERTIES LINKER_LANGUAGE CXX)
ELSE()
    add_library(${this_target} INTERFACE)
ENDIF()

ADD_EXECUTABLE(wila_tests tests/common.cpp tests/main.cpp)
TARGET_LINK_LIBRARIES(wila_tests
    pthread
    ${this_target}
    ${Boost_LIBRARIES}
    TBB::tbb
)

# Files to be installed
INSTALL(FILES ${${this_target}_HEADER_FILES} DESTINATION include/wila)
configure_file(${TARGET_LIB_NAME}.pc.in ${TARGET_LIB_NAME}.pc @ONLY)
# INSTALL(TARGETS wila_tests RUNTIME DESTINATION bin)
install(FILES ${CMAKE_BINARY_DIR}/${TARGET_LIB_NAME}.pc
    # DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig) from examplar, but this appears not prevalent, most I see are under /usr/lib/x86_64-linux-gnu/pkgconfig, so:
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
    # The above also prevents a lintian error

########### Add uninstall target ###############
# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
########### End uninstall target ###############


# Packages - currently just an RPM
set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "RPM;TGZ")
set(CPACK_PACKAGE_NAME "wila")
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_RELEASE ${CPACK_RPM_PACKAGE_RELEASE})
set(CPACK_PACKAGE_CONTACT "J-M <jean-michel.perraud@csiro.au>")
set(CPACK_PACKAGE_VENDOR "https://github.com/csiro-hydroinformatics/wila")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
include(CPack)

# Add new build types
message("* Adding build types...")
SET(CMAKE_CXX_FLAGS_COVERAGE
    "${GCC_DEBUG_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_C_FLAGS_COVERAGE
    "${GCC_DEBUG_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C compiler during coverage builds."
    FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used for linking binaries during coverage builds."
    FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used by the shared libraries linker during coverage builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )
