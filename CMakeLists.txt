cmake_minimum_required (VERSION 2.6)

project(ga_example)

if(!CMAKE_COMPILER_IS_GNUCC)
	message(SEND_ERROR "Please compile using GCC toolchain.")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_CXX_FLAGS}")
else (OPENMP_FOUND)
	message(SEND_ERROR "Your GCC does not support parallel API OpenMP.")
endif()

set(COMPILE_FLAGS "-Wall -std=c11")

configure_file (
	"src/config.h.in"
	"src/config.h"
)

add_executable(ga_example src/ga_example.c)

