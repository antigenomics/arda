cmake_minimum_required(VERSION 3.18)
project(arda LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Use pybind11's modern FindPython path (silences CMP0148 deprecation warnings).
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# Markup-transfer hot path: project reference region coordinates onto a query
# sequence by walking the alignment CIGAR. See src/_markup/markup.cpp.
pybind11_add_module(_markup src/_markup/markup.cpp)
target_compile_options(_markup PRIVATE -O3)

install(TARGETS _markup DESTINATION arda)
