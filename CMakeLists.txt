cmake_minimum_required(VERSION 3.8)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
include(FetchContent)

# Set name of project and library
set(PROJECT_NAME hephaestus)
set(PROJECT_LIB_NAME ${PROJECT_NAME}-lib)

project(${PROJECT_NAME} VERSION 0.0.0 LANGUAGES CXX)

############################################################
# Find and import MFEM 
############################################################
# The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for MFEM ...")

set(MFEM_DIR "${PROJECT_SOURCE_DIR}/../mfem/build" CACHE PATH "Path to the MFEM build or install prefix.")
set(MFEM_COMMON_INCLUDES "${MFEM_DIR}/../miniapps/common" CACHE PATH "Path to the MFEM common miniapp headers.")

if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()

############################################################
# Check C++ compiler is the same as used in MFEM 
############################################################
message(STATUS "Found MFEM config in: ${mfem_DIR} (version ${MFEM_VERSION})")

# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
message(STATUS "MFEM compiler is: ${MFEM_CXX_COMPILER}")

if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

message(STATUS "MFEM_INCLUDE_DIRS is set to: ${MFEM_INCLUDE_DIRS}")
message(STATUS "MFEM_COMMON_INCLUDES is set to: ${MFEM_COMMON_INCLUDES}")

find_library(MFEM_COMMON_LIBRARY mfem-common HINTS ${mfem_DIR}/miniapps/common)
message(STATUS "MFEM_COMMON_LIBRARY is found in the folder: ${MFEM_COMMON_LIBRARY}")

enable_testing()

## put all targets in bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
## put all libraries in lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
## export compile commands for clang-format
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

############################################################
# Fetch spdlog for logging purposes
############################################################
set(SPDLOG_BUILD_PIC ON)

FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.x
)

FetchContent_MakeAvailable(spdlog)

############################################################
# Set the build type environment variable
############################################################
if(BUILD_TYPE STREQUAL "TEST")
  set(BUILD_TYPE_COMPILE_FLAGS "-g;-O0;--coverage")
  set(TEST_LIBRARIES "gcov")
elseif(BUILD_TYPE STREQUAL "PROFILE")
  set(BUILD_TYPE_COMPILE_FLAGS "-g;-O0;-pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg") 
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg") 
  set(TEST_LIBRARIES "gcov")
elseif(BUILD_TYPE STREQUAL "DEBUG")
  set(BUILD_TYPE_COMPILE_FLAGS "-g;-O0")
  find_program(CLANG_TIDY NAMES "clang-tidy" "clang-tidy-10" "clang-tidy-11" REQUIRED)
  set (CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY} -checks=* CACHE STRING "Specify the clang-tidy command line")
  set(TEST_LIBRARIES "")
elseif(BUILD_TYPE STREQUAL "RELEASE")
  set(BUILD_TYPE_COMPILE_FLAGS "-O2")
  set(TEST_LIBRARIES "")
endif()

# Check clang is being used if build type is set to "DEBUG-MEMORY"
if(BUILD_TYPE STREQUAL "DEBUG-MEMORY" AND NOT ${CMAKE_CXX_COMPILER} STREQUAL "/usr/bin/clang++")
    message(FATAL_ERROR "DEBUG-MEMORY option requires the clang++ compiler.")
endif()

############################################################
# Add directories with source files
############################################################
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)

############################################################
# all configurable options
############################################################
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

############################################################
# Create doxygen documentation `make doc`
############################################################
if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
