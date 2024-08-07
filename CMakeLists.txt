cmake_minimum_required(VERSION 3.5)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)

set(PROJECT_NAME mc_dynamic_polytopes)
set(PROJECT_DESCRIPTION "mc_dynamic_polytopes")
set(PROJECT_URL "")

project(${PROJECT_NAME} CXX)

# Check if the project is built inside mc_rtc
if(NOT TARGET mc_rtc::mc_control)
  find_package(mc_rtc REQUIRED)
endif()

find_package(polito REQUIRED)


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
  ${PROJECT_SOURCE_DIR}/src/DynamicPolytope.cpp
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
  polito
)

#### Exporting targets and installing ####

include(GNUInstallDirs)
install(TARGETS DynamicPolytope
    EXPORT DynamicPolytope-export
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/
)

install(EXPORT DynamicPolytope-export
  FILE mc_dynamic_polytopesConfig.cmake
  NAMESPACE mc_dynamic_polytopes::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mc_dynamic_polytopes
)

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mc_dynamic_polytopes
)
