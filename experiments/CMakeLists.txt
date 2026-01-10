# Experiments (only adds new executables; does not modify src/test/thirdparty/include)

add_executable(msb5_fail_trace msb5_fail_trace.cpp)
target_link_libraries(msb5_fail_trace HEDB)
target_include_directories(msb5_fail_trace PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_definitions(msb5_fail_trace PRIVATE USE_RANDEN)

add_executable(msb_wraparound_experiment msb_wraparound_experiment.cpp)
target_link_libraries(msb_wraparound_experiment tfhe++)
target_include_directories(msb_wraparound_experiment PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_definitions(msb_wraparound_experiment PRIVATE USE_RANDEN)

add_executable(msb6_test msb6_test.cpp)
target_link_libraries(msb6_test HEDB)
target_include_directories(msb6_test PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_definitions(msb6_test PRIVATE USE_RANDEN)

add_executable(msb_timing_bench msb_timing_bench.cpp)
target_link_libraries(msb_timing_bench HEDB)
target_include_directories(msb_timing_bench PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_definitions(msb_timing_bench PRIVATE USE_RANDEN)

add_executable(msb_delta_half_fail_trace msb_delta_half_fail_trace.cpp)
target_link_libraries(msb_delta_half_fail_trace HEDB)
target_include_directories(msb_delta_half_fail_trace PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_definitions(msb_delta_half_fail_trace PRIVATE USE_RANDEN)
