cmake_minimum_required(VERSION 3.20)
project(example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(sb2l REQUIRED)

add_executable(example example.cpp vibes.cpp)

target_link_libraries(example
    PRIVATE
        sb2l
        qhullcpp
        symengine
        gmp
)