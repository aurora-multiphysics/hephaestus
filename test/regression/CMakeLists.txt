file(GLOB_RECURSE src_regression "*.cpp")

add_executable(regression_tests ${src_regression})
add_compile_options(regression_tests ${BUILD_TYPE_COMPILER_FLAGS})
target_include_directories(regression_tests PRIVATE ${GTEST_INCLUDE_DIRS})
target_include_directories(regression_tests PRIVATE ${GMOCK_INCLUDE_DIRS})
target_include_directories(regression_tests PUBLIC ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_include_directories(regression_tests PUBLIC ${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}_lib)

target_link_libraries(regression_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} pthread)
target_link_libraries(regression_tests ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(regression_tests ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)

include(GoogleTest)
gtest_discover_tests(regression_tests)
