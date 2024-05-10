file(GLOB_RECURSE src_unit "*.cpp")
file(GLOB_RECURSE test_src_files "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cpp")

set (${PROJECT_NAME}_INCLUDE_DIRS "")
foreach (_srcFile ${test_src_files})
    get_filename_component(_dir ${_srcFile} PATH)
    list (APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${_dir})
endforeach()
list (REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)

add_executable(integration_tests ${src_unit} Main.cpp)
add_compile_options(integration_tests ${BUILD_TYPE_COMPILER_FLAGS})
target_include_directories(integration_tests PUBLIC ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_include_directories(integration_tests PUBLIC ${${PROJECT_NAME}_INCLUDE_DIRS})
set_property(TARGET integration_tests PROPERTY CXX_STANDARD 17)

target_link_libraries(integration_tests spdlog::spdlog)
target_link_libraries(integration_tests Catch2::Catch2)
target_link_libraries(integration_tests ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(integration_tests ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)

include(CTest)
include(Catch)
catch_discover_tests(integration_tests WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/bin" DISCOVERY_MODE PRE_TEST)
