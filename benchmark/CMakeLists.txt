find_package(benchmark REQUIRED)

# Configure Google Benchmark to skip its own tests and docs.
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

add_executable(polytope_bench bench.cpp)

target_link_libraries(polytope_bench PRIVATE
  benchmark::benchmark
  benchmark::benchmark_main
  DynamicPolytope
)
