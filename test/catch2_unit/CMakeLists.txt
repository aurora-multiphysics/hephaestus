# get hephaestsus src headers for use in tests
file(GLOB_RECURSE test_src_files "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cpp")

set (${PROJECT_NAME}_INCLUDE_DIRS "")
foreach (_srcFile ${test_src_files})
    get_filename_component(_dir ${_srcFile} PATH)
    list (APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${_dir})
endforeach()
list (REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)

file(GLOB_RECURSE src_unit_catch2 "*.cpp") # get all files in catch2 test directory

# add_executable(src_unit_catch2 "${PROJECT_SOURCE_DIR}/test/catch2_unit/test_inputs_catch2.cpp") # for single file executable
add_executable(unit_tests_catch2 ${src_unit_catch2})
set_property(TARGET unit_tests_catch2 PROPERTY CXX_STANDARD 17)
target_link_libraries(unit_tests_catch2 Catch2::Catch2WithMain)
add_compile_options(unit_tests_catch2 ${BUILD_TYPE_COMPILER_FLAGS})


target_include_directories(unit_tests_catch2 PUBLIC ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_include_directories(unit_tests_catch2 PUBLIC ${${PROJECT_NAME}_INCLUDE_DIRS})
target_link_libraries(unit_tests_catch2 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(unit_tests_catch2 ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)