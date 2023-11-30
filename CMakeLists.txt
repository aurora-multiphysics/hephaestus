cmake_minimum_required(VERSION 3.8)

# no space in these names! used as shared object name and python module name
set(PROJECT_NAME hephaestus)
set(PROJECT_LIB_NAME ${PROJECT_NAME}-lib)

project(${PROJECT_NAME} VERSION 0.0.0 LANGUAGES CXX)
find_package(Boost REQUIRED COMPONENTS ${Boost_PACKAGES})

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")

set(MFEM_DIR "${PROJECT_SOURCE_DIR}/../mfem_build" CACHE PATH "Path to the MFEM build or install prefix.")
set(MFEM_COMMON_INCLUDES "${PROJECT_SOURCE_DIR}/../mfem/miniapps/common" CACHE PATH "Path to the MFEM common miniapp headers.")

if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
message(STATUS "Mfem compiler is: ${MFEM_CXX_COMPILER}")
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

message(STATUS "MFEM_INCLUDE_DIRS is set to: ${MFEM_INCLUDE_DIRS}")
message(STATUS "MFEM_LIBRARIES is set to: ${MFEM_LIBRARIES}")

find_library(MFEM_COMMON_LIBRARY mfem-common HINTS ${mfem_DIR}/miniapps/common)
message(STATUS "MFEM_COMMON_LIBRARY is found in the folder: ${MFEM_COMMON_LIBRARY}")

enable_testing()

## put all targets in bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
## put all libraries in lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

############################################################
# Set the build type environment variable
############################################################
if(BUILD_TYPE STREQUAL "TEST")
  set(BUILD_TYPE_COMPILE_FLAGS "-g;-O0;--coverage")
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
