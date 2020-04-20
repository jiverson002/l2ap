# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.0)

# ...
project(l2ap
  VERSION      0.2.1
  DESCRIPTION  "High performance sparse fixed-radius library"
  LANGUAGES    C)

# Include required CMake modules
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

#-------------------------------------------------------------------------------
# DEPENDENCY configuration
#-------------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

add_subdirectory(${PROJECT_SOURCE_DIR}/GKlib)

#-------------------------------------------------------------------------------
# LIBRARY configuration
#-------------------------------------------------------------------------------
add_library(${PROJECT_NAME})

target_compile_definitions(${PROJECT_NAME}
  PRIVATE L2PS RS4 L2CG PSCV DP5 L2CV TL1 TL2 MKPL)

target_sources(${PROJECT_NAME}
  PRIVATE ap.c cmdline.c idxjoin.c l2ap-t2.c memory.c select.c l2ap-m.c l2ap.c
          mkj.c shim.c sort.c csr.c l2ap-t.c mmjoin.c util.c)

target_include_directories(${PROJECT_NAME}
  PRIVATE ${GKlib_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} PUBLIC GKlib OpenMP::OpenMP_C)

add_library(L2AP::l2ap ALIAS ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME} PROPERTIES
  SOVERSION ${PROJECT_VERSION_MAJOR}
  VERSION   ${PROJECT_VERSION})

#-------------------------------------------------------------------------------
# PACKAGE configuration
#-------------------------------------------------------------------------------
configure_package_config_file(L2APConfig.cmake.in cmake/L2APConfig.cmake
  INSTALL_DESTINATION lib/cmake/L2AP)

write_basic_package_version_file(cmake/L2APConfigVersion.cmake
  VERSION       ${PROJECT_VERSION}
  COMPATIBILITY ExactVersion)

install(TARGETS ${PROJECT_NAME} GKlib EXPORT L2APTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
           COMPONENT   L2AP_Runtime
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           COMPONENT          L2AP_Runtime
           NAMELINK_COMPONENT L2AP_Development
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           COMPONENT   L2AP_Development
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT L2APTargets
  FILE        L2APTargets.cmake
  NAMESPACE   L2AP::
  DESTINATION lib/cmake/L2AP
  COMPONENT   L2AP_Development)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/"
  DESTINATION    ${CMAKE_INSTALL_INCLUDEDIR}/l2ap
  COMPONENT      L2AP_Development
  FILES_MATCHING
  PATTERN        "*[.]h"
  PATTERN        "build/GKlib"
  EXCLUDE
  PATTERN        "cmakebuild"
  EXCLUDE
  PATTERN        "GKlib"
  EXCLUDE)

install(DIRECTORY ${PROJECT_INCLUDE_DIRECTORIES}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/l2ap
  COMPONENT   L2AP_Development)

install(
  FILES       ${CMAKE_CURRENT_BINARY_DIR}/cmake/L2APConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/cmake/L2APConfigVersion.cmake
  DESTINATION lib/cmake/L2AP
  COMPONENT   L2AP_Development)
