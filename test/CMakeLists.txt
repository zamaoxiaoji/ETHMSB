add_executable(comparison_test comparison_test.cpp)
target_link_libraries(comparison_test HEDB)

add_executable(conversion_test conversion_test.cpp)
target_link_libraries(conversion_test HEDB)

add_executable(tpch_q6 tpch_q6.cpp)
target_link_libraries(tpch_q6 HEDB)

add_executable(tpch_q6_copy tpch_q6_copy.cpp)
target_link_libraries(tpch_q6_copy HEDB)

add_executable(comparison_test_copy comparison_test_copy.cpp)
target_link_libraries(comparison_test_copy HEDB)

add_executable(tpch_q6_mytest1 tpch_q6_mytest1.cpp)
target_link_libraries(tpch_q6_mytest1 HEDB)

add_executable(add_test add_test.cpp)
target_link_libraries(add_test HEDB)

add_executable(add_replace_and_test add_replace_and_test.cpp)
target_link_libraries(add_replace_and_test HEDB)

add_executable(tpch_q6_mytest_v1 tpch_q6_mytest_v1.cpp)
target_link_libraries(tpch_q6_mytest_v1 HEDB)

add_executable(tpch_q6_mytest_v2 tpch_q6_mytest_v2.cpp)
target_link_libraries(tpch_q6_mytest_v2 HEDB)

add_executable(tpch_q6_mytest_v3 tpch_q6_mytest_v3.cpp)
target_link_libraries(tpch_q6_mytest_v3 HEDB)

add_executable(tpch_q6_mytest_v4 tpch_q6_mytest_v4.cpp)
target_link_libraries(tpch_q6_mytest_v4 HEDB)

add_executable(tpch_q6_mytest_v5 tpch_q6_mytest_v5.cpp)
target_link_libraries(tpch_q6_mytest_v5 HEDB)

add_executable(lwe_mul lwe_mul.cpp)
target_link_libraries(lwe_mul HEDB)

add_executable(tpch_q14 tpch_q14.cpp)
target_link_libraries(tpch_q14 HEDB)

add_executable(tpch_q14_v1 tpch_q14_v1.cpp)
target_link_libraries(tpch_q14_v1 HEDB)

add_executable(tpch_q14_v0 tpch_q14_v0.cpp)
target_link_libraries(tpch_q14_v0 HEDB)

add_executable(tpch_q14_v0_copy tpch_q14_v0_copy.cpp)
target_link_libraries(tpch_q14_v0_copy HEDB)

add_executable(add_test_copy add_test_copy.cpp)
target_link_libraries(add_test_copy HEDB)

add_executable(nested_cmux_from_cb nested_cmux_from_cb.cpp)
target_link_libraries(nested_cmux_from_cb HEDB)

add_executable(tpch_q19 tpch_q19.cpp)
target_link_libraries(tpch_q19 HEDB)

add_executable(tpch_q19_v1 tpch_q19_v1.cpp)
target_link_libraries(tpch_q19_v1 HEDB)

add_executable(distinct_test distinct_test.cpp)
target_link_libraries(distinct_test HEDB)

add_executable(distinct_test_v2 distinct_test_v2.cpp)
target_link_libraries(distinct_test_v2 HEDB)

add_executable(distinct_test_v3 distinct_test_v3.cpp)
target_link_libraries(distinct_test_v3 HEDB)

add_executable(new_test new_test.cpp)
target_link_libraries(new_test HEDB)

add_executable(new_test2 new_test2.cpp)
target_link_libraries(new_test2 HEDB)

add_executable(com_test com_test.cpp)
target_link_libraries(com_test HEDB)
