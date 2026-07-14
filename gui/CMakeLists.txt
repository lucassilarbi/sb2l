cmake_minimum_required(VERSION 3.20)

find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)

set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/3rd/imgui)

# Vendored Dear ImGui (core + GLFW/OpenGL3 backends). The opengl3 backend ships
# its own GL loader, so no glad/glew is needed.
add_library(imgui STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui SYSTEM PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends)
target_link_libraries(imgui PUBLIC glfw OpenGL::GL)
target_compile_options(imgui PRIVATE -w)

add_executable(sb2l_gui
    main.cpp
    app.cpp
    editor.cpp
)
target_compile_options(sb2l_gui PRIVATE -frounding-math)
# app.cpp builds/reads filib affine forms inline. GCC's optimizer (any -O >= 1,
# even with -frounding-math) miscompiles that interval/affine arithmetic and the
# zonotopes collapse, so build this one translation unit unoptimized. The heavy
# work stays in the optimized sb2l/ibex libraries.
set_source_files_properties(app.cpp PROPERTIES COMPILE_OPTIONS "-O0")
target_link_libraries(sb2l_gui PRIVATE sb2l imgui glfw OpenGL::GL)
