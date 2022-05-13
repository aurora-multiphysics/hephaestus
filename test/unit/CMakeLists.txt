file(GLOB_RECURSE src_unit "*.cpp")

add_executable(unit_tests ${src_unit})
add_compile_options(unit_tests ${BUILD_TYPE_COMPILER_FLAGS})
target_include_directories(unit_tests PRIVATE ${GTEST_INCLUDE_DIRS})
target_include_directories(unit_tests PRIVATE ${GMOCK_INCLUDE_DIRS})
target_include_directories(unit_tests PUBLIC ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_include_directories(unit_tests PUBLIC ${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}_lib)

target_link_libraries(unit_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} pthread)
target_link_libraries(unit_tests ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(unit_tests ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)

include(GoogleTest)
gtest_discover_tests(unit_tests)
