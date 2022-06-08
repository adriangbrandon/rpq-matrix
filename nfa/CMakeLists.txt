cmake_minimum_required(VERSION 2.8)
project(NFA)
set(PROJECT_VENDOR "Adrian Gomez Brandon")
set(CMAKE_PROJECT_NAME "NFA")
set(PROJECT_CONTACT "adrian.gbrandon@udc.es")
set(CMAKE_INSTALL_PREFIX /home/adrian/software)
set(CMAKE_COMPILER_IS_GNUCXX)

include_directories(./include)

set (SOURCES
        include/basics.hpp
        src/bitmasks.cpp
        include/bitmasks.hpp
        src/options.cpp
        include/options.hpp
        src/parser.cpp
        include/parser.hpp
        src/regular.cpp
        include/regular.hpp
        src/RpqAutomata.cpp
        include/RpqAutomata.hpp
        include/RpqTree.hpp
        src/RpqTree.cpp
        )


# Flags
set(CMAKE_CXX_FLAGS " -Wall -m64")
set(CMAKE_CXX_FLAGS_RELEASE " -std=c++0x -Wall -O9 -DNDEBUG")
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
add_library(nfa STATIC ${SOURCES})
install(TARGETS nfa DESTINATION lib)
install (DIRECTORY "include" DESTINATION src)