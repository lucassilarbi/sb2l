include(CMakeFindDependencyMacro)

set(SymEngine_DIR "${CMAKE_CURRENT_LIST_DIR}/../symengine")
find_dependency(SymEngine REQUIRED)

# find gmp, necessary because symengine's cmake never finds it
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_dependency(GMP REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/sb2lTargets.cmake")

add_library(ibex STATIC IMPORTED)

set_target_properties(ibex PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_LIST_DIR}/../../../lib/libibex.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_LIST_DIR}/../../../include"
)

add_library(filib STATIC IMPORTED)

set_target_properties(filib PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_LIST_DIR}/../../../lib/libprim.a"
)

target_link_libraries(ibex INTERFACE filib)