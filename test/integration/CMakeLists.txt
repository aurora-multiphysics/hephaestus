file(GLOB_RECURSE src_integration "*.cpp")

add_executable(integration_tests ${src_integration})
add_compile_options(integration_tests ${BUILD_TYPE_COMPILER_FLAGS})
target_include_directories(integration_tests PRIVATE ${GTEST_INCLUDE_DIRS})
target_include_directories(integration_tests PRIVATE ${GMOCK_INCLUDE_DIRS})
target_include_directories(integration_tests PUBLIC ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_include_directories(integration_tests PUBLIC ${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}_lib)
set_property(TARGET integration_tests PROPERTY CXX_STANDARD 17)

target_link_libraries(integration_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} pthread)
target_link_libraries(integration_tests ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(integration_tests ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)

include(GoogleTest)
gtest_discover_tests(integration_tests EXTRA_ARGS --data_directory ${PROJECT_SOURCE_DIR}/data/)
