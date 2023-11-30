Find_package(Threads REQUIRED)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)


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

