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

# Exact k-mer prefilter: reject reads that cannot align before MMseqs2 sees them.
# See src/_prefilter/prefilter.cpp and OPTIMIZATION.md sections 3-4.
pybind11_add_module(_prefilter src/_prefilter/prefilter.cpp)
target_compile_options(_prefilter PRIVATE -O3)
find_package(Threads REQUIRED)
target_link_libraries(_prefilter PRIVATE Threads::Threads)

install(TARGETS _prefilter DESTINATION arda)

# Structure-aware segment mapper: best V and best J per read WITHOUT a homology search.
# Deliberately a SEPARATE module from _prefilter -- that one is shipped, measured and load-bearing,
# and this is opt-in until it is proven not to move the output.
# See src/_segmap/segmap.cpp and results/round10/aligner_options.md in the benchmark repo.
pybind11_add_module(_segmap src/_segmap/segmap.cpp)
target_compile_options(_segmap PRIVATE -O3)
target_link_libraries(_segmap PRIVATE Threads::Threads)

install(TARGETS _segmap DESTINATION arda)
