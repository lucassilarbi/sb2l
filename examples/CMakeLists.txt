cmake_minimum_required(VERSION 3.20)
project(example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(sb2l REQUIRED)

add_executable(basic_example basic_example.cpp vibes.cpp)
add_executable(advanced_example advanced_example.cpp vibes.cpp)

foreach(target basic_example advanced_example)
    target_link_libraries(${target}
        PRIVATE
            sb2l
            symengine
            gmp
    )
endforeach()