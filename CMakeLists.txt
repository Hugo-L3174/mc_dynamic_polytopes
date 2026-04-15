cmake_minimum_required(VERSION 3.22)

include(cmake/get-jrl-cmakemodules.cmake)

project(
  mc_dynamic_polytopes
  VERSION 1.0.1
  DESCRIPTION "mc_rtc library for dynamic balance polytopes computations using politopix"
  HOMEPAGE_URL "https://github.com/Hugo-L3174/mc_dynamic_polytopes"
  LANGUAGES CXX
)

jrl_configure_defaults()

# cmake macro for project dependencies
jrl_find_package(polito REQUIRED NO_MODULE)
jrl_find_package(Qhull REQUIRED)

# Check if the project is built inside mc_rtc
if(NOT TARGET mc_rtc::mc_control)
  jrl_find_package(mc_rtc REQUIRED)
endif()


set(HEADERS
  include/mc_dynamic_polytopes/AsyncJob.h
  include/mc_dynamic_polytopes/ContactPolytope.h
  include/mc_dynamic_polytopes/DynamicPolytope.h
  include/mc_dynamic_polytopes/GUIComputations.h
  include/mc_dynamic_polytopes/PolytopeFuncs.h
  include/mc_dynamic_polytopes/VRPRegionMinkSum.h
  include/mc_dynamic_polytopes/WrenchCones.h
  include/mc_dynamic_polytopes/ZMPRegion.h
)

add_library(DynamicPolytope SHARED
  ${PROJECT_SOURCE_DIR}/src/WrenchCones.cpp
  ${PROJECT_SOURCE_DIR}/src/GUIComputations.cpp
  ${PROJECT_SOURCE_DIR}/src/DynamicPolytope.cpp
  ${PROJECT_SOURCE_DIR}/src/ContactPolytope.cpp
  ${PROJECT_SOURCE_DIR}/src/PolytopeFuncs.cpp
  ${PROJECT_SOURCE_DIR}/src/ZMPRegion.cpp
  ${PROJECT_SOURCE_DIR}/src/VRPRegionMinkSum.cpp
)

target_include_directories(DynamicPolytope
  PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

jrl_target_headers(DynamicPolytope PUBLIC HEADERS ${HEADERS} BASE_DIRS include)
set_target_properties(DynamicPolytope PROPERTIES VERSION ${PROJECT_VERSION})


target_link_libraries(DynamicPolytope
  PUBLIC
  mc_rtc::mc_control_fsm
  Politopix::polito
  Qhull::qhull_r
)

add_subdirectory(tests)

#### Exporting targets and installing ####
jrl_add_export_component(NAME DynamicPolytope TARGETS DynamicPolytope)
jrl_export_package()


# jrl_print_dependencies_summary()
jrl_print_options_summary()
