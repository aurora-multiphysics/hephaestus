Find_package(Threads REQUIRED)

include(ExternalProject)

ExternalProject_Add(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

#ExternalProject_Get_Property(Catch2 source_dir)
#set(CATCH2_INCLUDE_DIRS ${source_dir}/Catch2/include)
#
#ExternalProject_Get_Property(Catch2 binary_dir)
#set(CATCH2_LIBRARY_PATH  ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}libCatch2.a)
#set(CATCH2_MAIN_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}libCatch2Main.a)
#set(CATCH2_LIBRARY CATCH2)
#set(CATCH2_MAIN_LIBRARY CATCH2_main)
#add_library(${CATCH2_LIBRARY} UNKNOWN IMPORTED)
#add_library(${CATCH2_MAIN_LIBRARY} UNKNOWN IMPORTED)
#set_target_properties(${CATCH2_LIBRARY} PROPERTIES IMPORTED_LOCATION ${CATCH2_LIBRARY_PATH})
#set_target_properties(${CATCH2_MAIN_LIBRARY} PROPERTIES IMPORTED_LOCATION ${CATCH2_MAIN_LIBRARY_PATH})
#add_dependencies(${CATCH2_LIBRARY} Catch2)
#add_dependencies(${CATCH2_MAIN_LIBRARY} Catch2)

# Set up to generate coverage reports
find_package(PythonInterp REQUIRED)
find_program(GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/scripts/test)

IF(BUILD_TYPE STREQUAL "TEST")
  add_custom_target(gcovr
    # Create folder
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/gcovr

    # Running gcovr
    COMMAND ${Python_EXECUTABLE} ${GCOVR_PATH} --html --html-details -s
            -r ${PROJECT_SOURCE_DIR} --object-directory=${PROJECT_BINARY_DIR}
            -o gcovr/index.html
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Running gcovr to produce HTML code coverage report."
  )
  # Show info where to find the report
  add_custom_command(TARGET gcovr POST_BUILD
    COMMAND ;
    COMMENT "Open ./gcovr/index.html in your browser to view the full coverage report."
  )
ENDIF()

add_subdirectory(unit)
add_subdirectory(integration)
add_subdirectory(regression)

