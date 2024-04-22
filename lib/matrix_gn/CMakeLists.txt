cmake_minimum_required(VERSION 2.8)
set(CMAKE_OSX_ARCHITECTURES x86_64)
project(MatrixGN)
set(CMAKE_PROJECT_NAME "MatrixGN")
#set(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_OSX_ARCHITECTURES x86_64)

include_directories(include)

set(SOURCES
        include/k2-tree/basics.h
        include/k2-tree/utilstime.h
        include/k2-tree/bitvector.h
        include/k2-tree/k2tree.h
        include/k2-tree/matrix.h
        src/k2-tree/basics.c
        src/k2-tree/bitvector.c
        src/k2-tree/k2tree.c
        src/k2-tree/matrix.c
        src/k2-tree/utilstime.c
        )

# Flags
set(CMAKE_CXX_FLAGS " -Wall -m64")
set(CMAKE_CXX_FLAGS_RELEASE " -std=c++0x -Wall -O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -O9 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -g3 -DDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -O3 -pg -DDEBUG")

# set default compilation and show details for debug build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif(NOT CMAKE_BUILD_TYPE)

# Create the library
add_library(matrix_gn STATIC ${SOURCES})

#add_executable(proba main2.cpp ${SOURCES})
