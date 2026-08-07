include(CMakeFindDependencyMacro)

set(SymEngine_DIR "${CMAKE_CURRENT_LIST_DIR}/../symengine")
find_dependency(SymEngine REQUIRED)

# find gmp, necessary because symengine's cmake never finds it
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_dependency(GMP REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/sb2lTargets.cmake")

get_filename_component(SB2L_PREFIX "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)

# Guarded, and GLOBAL: this file may be read again by a second find_package in
# another directory of the consuming project, and the targets have to be
# visible from every one of them.
#
# DynIbex only builds and installs its own filib when it finds none on the
# system, so libprim.a is not always part of this package. When it is missing,
# the one DynIbex linked against is on the system and is asked for by name.
if(NOT TARGET filib)
    add_library(filib INTERFACE IMPORTED GLOBAL)
    if(EXISTS "${SB2L_PREFIX}/lib/libprim.a")
        set_target_properties(filib PROPERTIES
            INTERFACE_LINK_LIBRARIES "${SB2L_PREFIX}/lib/libprim.a"
        )
    else()
        set_target_properties(filib PROPERTIES INTERFACE_LINK_LIBRARIES prim)
    endif()
endif()

if(NOT TARGET ibex)
    add_library(ibex STATIC IMPORTED GLOBAL)
    set_target_properties(ibex PROPERTIES
        IMPORTED_LOCATION "${SB2L_PREFIX}/lib/libibex.a"
        INTERFACE_INCLUDE_DIRECTORIES "${SB2L_PREFIX}/include"
    )
    target_link_libraries(ibex INTERFACE filib)
endif()
