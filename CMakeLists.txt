cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)
set(CMAKE_OSX_ARCHITECTURES x86_64)
project(RPQ-Matrix)
set(CMAKE_PROJECT_NAME "RPQ-Matrix")

#Change to "Debug" when you need it
set(CMAKE_BUILD_TYPE "Release")

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
        ${CMAKE_HOME_DIRECTORY}/include
        ${CMAKE_HOME_DIRECTORY}/nfa/include
        ${CMAKE_HOME_DIRECTORY}/lib/matrix_gn/include
        ${CMAKE_HOME_DIRECTORY}/lib/matrix_p_gn/include
        ${CMAKE_HOME_DIRECTORY}/lib/baseline_gn/include)

#link_directories(~/lib)

add_subdirectory(nfa)
add_subdirectory(lib/matrix_gn)
add_subdirectory(lib/matrix_p_gn)
add_subdirectory(lib/baseline_gn)

add_executable(k2_tree_query src/k2_tree_query.cpp)
target_link_libraries(k2_tree_query matrix_gn nfa)

add_executable(k2_tree_p_query src/k2_tree_p_query.cpp)
target_link_libraries(k2_tree_p_query matrix_p_gn nfa)

add_executable(baseline_query src/baseline_query.cpp)
target_link_libraries(baseline_query baseline_gn nfa)

add_executable(baseline_32_query src/baseline_32_query.cpp)
target_link_libraries(baseline_32_query baseline_gn nfa)

add_executable(create_pairs src/create_pairs.cpp)
target_link_libraries(create_pairs matrix_gn)

add_executable(k2_tree_build src/k2_tree_build.cpp)
target_link_libraries(k2_tree_build matrix_gn)

add_executable(baseline_32_build src/baseline_32_build.cpp)
target_link_libraries(baseline_32_build baseline_gn)

add_executable(baseline_build src/baseline_build.cpp)
target_link_libraries(baseline_build baseline_gn)