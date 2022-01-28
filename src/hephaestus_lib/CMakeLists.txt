file(GLOB_RECURSE src_files "*.h" "*.hpp" "*.cpp")

message(STATUS "src files:")
message(STATUS ${src_files})

# make a shared object library to link against for testing
add_library(${PROJECT_LIB_NAME} SHARED ${src_files})
target_compile_options(${PROJECT_LIB_NAME} PRIVATE ${BUILD_TYPE_COMPILE_FLAGS})
set_target_properties(${PROJECT_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(${PROJECT_LIB_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_include_directories(${PROJECT_LIB_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${MFEM_COMMON_INCLUDES} ${MFEM_INCLUDE_DIRS})
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC ${Boost_LIBRARIES} ${TEST_LIBRARIES} ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)

