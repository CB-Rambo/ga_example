cmake_minimum_required (VERSION 2.6)

project(ga_example)

if(NOT CMAKE_COMPILER_IS_GNUCC)
	message(SEND_ERROR "Please compile using GCC toolchain.")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Wall")
#	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_CXX_FLAGS}")
#	set (COMPILE_FLAGS "${COMPILE_FLAGS} ${OPENMP_C_FLAGS}")
else (OPENMP_FOUND)
	message(SEND_ERROR "Your GCC does not support parallel API OpenMP.")
endif(OPENMP_FOUND)

configure_file (
	"src/config.h.in"
	"src/config.h"
)

add_executable(ga_example src/ga_example.c)

find_package(Doxygen)
if (DOXYGEN_FOUND)
	add_custom_target(doc 
		COMMAND ${DOXYGEN_EXECUTABLE} 
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
	)
else (DOXYGEN_FOUND)
	message(WARNING "Did not find Doxygen. Not building documentation")
endif(DOXYGEN_FOUND)