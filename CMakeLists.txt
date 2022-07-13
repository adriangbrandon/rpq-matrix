cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)
project(Ring-RPQ)
set(CMAKE_PROJECT_NAME "Ring RPQ")

#Change to "Debug" when you need it
set(CMAKE_BUILD_TYPE "Debug")

# C++11 compiler Check
if(NOT CMAKE_CXX_COMPILER_VERSION) # work around for cmake versions smaller than 2.8.10
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
endif()
if(CMAKE_CXX_COMPILER MATCHES ".*clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()
if( (CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.7) OR
(CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.2))
    message(FATAL_ERROR "Your C++ compiler does not support C++11. Please install g++ 4.7 (or greater) or clang 3.2 (or greater)")
else()
    message(STATUS "Compiler is recent enough to support C++11.")
endif()
if( CMAKE_COMPILER_IS_GNUCXX )
    append_cxx_compiler_flags("-fpermissive -std=c++11 -Wall -Wextra " "GCC" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3 -ffast-math -funroll-loops -fno-omit-frame-pointer -g" "GCC" CMAKE_CXX_FLAGS_RELEASE)
else()
    append_cxx_compiler_flags("-fpermissive -std=c++11" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3  -ffast-math -funroll-loops -DNDEBUG" "CLANG" CMAKE_CXX_FLAGS_RELEASE)
endif()
include(CheckSSE)
FindSSE ()
if( SSE4_2_FOUND )
    if( CMAKE_COMPILER_IS_GNUCXX )
        append_cxx_compiler_flags("-msse4.2" "GCC" CMAKE_CXX_FLAGS)
    else()
        append_cxx_compiler_flags("-msse4.2" "CLANG" CMAKE_CXX_FLAGS)
    endif()
    message(STATUS "CPU does support SSE4.2.")
else()
    message(STATUS "CPU does NOT support SSE4.2")
endif()

include_directories(~/include
        ${CMAKE_HOME_DIRECTORY}/nfa/include
        ${CMAKE_HOME_DIRECTORY}/external/malloc_count/include)

link_directories(~/lib)

add_subdirectory(nfa)
add_subdirectory(external/malloc_count)

add_executable(build-index src/build-index.cpp)
target_link_libraries(build-index sdsl divsufsort divsufsort64 nfa)

add_executable(query-index src/query-index.cpp)
target_link_libraries(query-index sdsl divsufsort divsufsort64 nfa)

add_executable(query-index-split src/query-index-split.cpp)
target_link_libraries(query-index-split sdsl divsufsort divsufsort64 nfa)

add_executable(split-benchmark src/split-benchmark.cpp)
target_link_libraries(split-benchmark sdsl divsufsort divsufsort64 nfa)

add_executable(split-benchmark_v3 src/split-benchmark_v3.cpp)
target_link_libraries(split-benchmark_v3 sdsl divsufsort divsufsort64 nfa)


add_executable(query-index-split_v3 src/query-index-split_v3.cpp)
target_link_libraries(query-index-split_v3 sdsl divsufsort divsufsort64 nfa)

add_executable(query-index-split_v3_v1 src/query-index-split_v3_1.cpp)
target_link_libraries(query-index-split_v3_v1 sdsl divsufsort divsufsort64 nfa)

add_executable(query-index-split_v4 src/query-index-split_v4.cpp)
target_link_libraries(query-index-split_v4 sdsl divsufsort divsufsort64 nfa)

add_executable(sort-dict src/sort-pos.cpp)
target_link_libraries(sort-dict sdsl divsufsort divsufsort64 nfa)

add_executable(sort-dict-bfs src/sort-pos-bfs.cpp)
target_link_libraries(sort-dict-bfs sdsl divsufsort divsufsort64 nfa)

add_executable(mandatory src/mandatory.cpp)
target_link_libraries(mandatory sdsl divsufsort divsufsort64 nfa)

add_executable(intersection_or src/intersection_or.cpp)
target_link_libraries(intersection_or sdsl divsufsort divsufsort64 nfa)

add_executable(pre-query src/pre-query.cpp)
target_link_libraries(pre-query sdsl divsufsort divsufsort64 nfa)

add_executable(prepare-yago2s src/prepare-yago2s.cpp)
target_link_libraries(prepare-yago2s sdsl divsufsort divsufsort64 nfa)

add_executable(prepare-queries-yago2s src/prepare-queries-yago2s.cpp)
target_link_libraries(prepare-queries-yago2s sdsl divsufsort divsufsort64 nfa)


add_executable(split-memory-bfs src/split-memory-bfs.cpp)
target_link_libraries(split-memory-bfs sdsl divsufsort divsufsort64 nfa malloc_count)

add_executable(split-memory-dfs src/split-memory-dfs.cpp)
target_link_libraries(split-memory-dfs sdsl divsufsort divsufsort64 nfa malloc_count)