cmake_minimum_required(VERSION 3.20)
project(sb2l VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

add_library(sb2l
    src/sb2l.cpp
)

# GMP (for symengine)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/3rd/symengine-0.14.0/cmake")
find_package(GMP REQUIRED)
target_include_directories(sb2l SYSTEM PRIVATE ${GMP_INCLUDE_DIRS})
target_link_libraries(sb2l PRIVATE gmp)

# Symengine
set(BUILD_TESTS no CACHE BOOL "Build SymEngine tests" FORCE)
set(BUILD_BENCHMARKS no CACHE BOOL "Build SymEngine benchmarks" FORCE)
set(WITH_SYMENGINE_TEUCHOS OFF CACHE BOOL "" FORCE)
set(WITH_SYMENGINE_RCP ON CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rd/symengine-0.14.0")
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${symengine_BINARY_DIR})
target_link_libraries(sb2l PUBLIC symengine)
target_compile_options(symengine PRIVATE -w)

# Dynibex
include(ExternalProject)
ExternalProject_Add(dynibex_ext
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/3rd/dynibex-2.3 ${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3
    CONFIGURE_COMMAND python2.7 ./waf configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3 --without-lp
    BUILD_COMMAND python2.7 ./waf build
    INSTALL_COMMAND python2.7 ./waf install
    BUILD_IN_SOURCE 1
)
add_library(filib STATIC IMPORTED GLOBAL)
set_target_properties(filib PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3/lib/libprim.a")
add_dependencies(filib dynibex_ext)
add_library(ibex STATIC IMPORTED GLOBAL)
set_target_properties(ibex PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3/lib/libibex.a")
add_dependencies(ibex dynibex_ext)
target_link_libraries(ibex INTERFACE filib)

target_include_directories(sb2l
    SYSTEM PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rd/symengine-0.14.0>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/3rd/symengine-0.14.0>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3/include>
)

target_include_directories(sb2l
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(sb2l PUBLIC ibex)

install(TARGETS sb2l EXPORT sb2lTargets)
install(DIRECTORY include/ DESTINATION include)
install(EXPORT sb2lTargets FILE sb2lTargets.cmake DESTINATION lib/cmake/sb2l)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/sb2lConfig.cmake" DESTINATION lib/cmake/sb2l)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3/include/ DESTINATION include)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3/lib/ DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/3rd/dynibex-2.3/share/ DESTINATION share)
# install cmakes to give access to GMP
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/3rd/symengine-0.14.0/cmake/LibFindMacros.cmake 
    FILES ${CMAKE_CURRENT_BINARY_DIR}/3rd/symengine-0.14.0/cmake/FindGMP.cmake 
    DESTINATION lib/cmake/sb2l)