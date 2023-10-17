file(GLOB_RECURSE example_src "*.cpp")

file(GLOB_RECURSE test_src_files "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cpp")

set (${PROJECT_NAME}_INCLUDE_DIRS "")
foreach (_srcFile ${test_src_files})
    get_filename_component(_dir ${_srcFile} PATH)
    list (APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${_dir})
endforeach()
list (REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)


add_executable(magnetostatic ${example_src})
add_compile_options(magnetostatic ${BUILD_TYPE_COMPILER_FLAGS})
target_include_directories(magnetostatic PUBLIC ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_include_directories(magnetostatic PUBLIC ${${PROJECT_NAME}_INCLUDE_DIRS})
target_include_directories(magnetostatic PUBLIC ${PROJECT_SOURCE_DIR}/data/)

set_property(TARGET magnetostatic PROPERTY CXX_STANDARD 17)

target_link_libraries(magnetostatic ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} pthread)
target_link_libraries(magnetostatic ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(magnetostatic ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)
