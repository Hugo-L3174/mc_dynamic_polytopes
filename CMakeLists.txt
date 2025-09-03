cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)

set(PROJECT_NAME mc_dynamic_polytopes)
set(PROJECT_VERSION 0.1.0)
set(PROJECT_DESCRIPTION "")
set(PROJECT_URL "")

# jrl-cmakemodules minimal settings
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CXX_DISABLE_WERROR 1)

include(cmake/base.cmake)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

# cmake macro for project dependencies
add_project_dependency(polito REQUIRED NO_MODULE)
add_project_dependency(Qhull REQUIRED)

# Check if the project is built inside mc_rtc
if(NOT TARGET mc_rtc::mc_control)
  find_package(mc_rtc REQUIRED)
endif()


# add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})
# target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
# set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DEIGEN_CDDLIB_EXPORTS")
# target_link_libraries(${PROJECT_NAME} PUBLIC cdd::cdd)
# target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")

# install(TARGETS ${PROJECT_NAME}
#     EXPORT "${TARGETS_EXPORT_NAME}"
#     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
# )
# install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})

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


target_link_libraries(DynamicPolytope
  PUBLIC
  mc_rtc::mc_control_fsm
  Politopix::polito
  Qhull::qhull_r
)

add_subdirectory(tests)

#### Exporting targets and installing ####

install(
  TARGETS DynamicPolytope
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(DIRECTORY include/mc_dynamic_polytopes
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
