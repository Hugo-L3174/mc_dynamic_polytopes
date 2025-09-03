find_package(Boost REQUIRED COMPONENTS unit_test_framework)


enable_testing()

add_executable(test_dynamic_polytope
    test_dynamic_polytope.cpp
)

target_link_libraries(test_dynamic_polytope PRIVATE DynamicPolytope)
target_link_libraries(
  test_dynamic_polytope PUBLIC Boost::unit_test_framework Boost::disable_autolinking Boost::dynamic_linking
)
target_compile_definitions(test_dynamic_polytope PRIVATE -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

add_test(NAME DynamicPolytopeTest COMMAND test_dynamic_polytope)
