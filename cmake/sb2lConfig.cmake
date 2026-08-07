include(CMakeFindDependencyMacro)

set(SymEngine_DIR "${CMAKE_CURRENT_LIST_DIR}/../symengine")
find_dependency(SymEngine REQUIRED)

# find gmp, necessary because symengine's cmake never finds it
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_dependency(GMP REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/sb2lTargets.cmake")

# Guarded, and GLOBAL: this file may be read again by a second find_package in
# another directory of the consuming project, and the targets have to be
# visible from every one of them.
if(NOT TARGET filib)
    add_library(filib STATIC IMPORTED GLOBAL)
    set_target_properties(filib PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_LIST_DIR}/../../../lib/libprim.a"
    )
endif()

if(NOT TARGET ibex)
    add_library(ibex STATIC IMPORTED GLOBAL)
    set_target_properties(ibex PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_LIST_DIR}/../../../lib/libibex.a"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_LIST_DIR}/../../../include"
    )
    target_link_libraries(ibex INTERFACE filib)
endif()
