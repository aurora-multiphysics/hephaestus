file(GLOB_RECURSE src_files "*.h" "*.hpp" "*.cpp")

message(STATUS "src files:")
message(STATUS ${src_files})

include_directories(${Boost_INCLUDE_DIRS})
add_executable(${PROJECT_NAME} ${src_files})

add_compile_options(${PROJECT_NAME} ${BUILD_TYPE_COMPILE_FLAGS})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}-lib ${MFEM_LIBRARIES} ${MFEM_COMMON_LIBRARY} -lrt)
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES} ${TEST_LIBRARIES} )



