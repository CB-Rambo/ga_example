cmake_minimum_required (VERSION 3.1)

project(gaalgo  C)

if(NOT(${CMAKE_C_COMPILER_ID} MATCHES "Clang|GNU"))
	message(SEND_ERROR "Please compile using GCC toolchain.")
endif()

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Wall -O3 -g")
#	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_CXX_FLAGS}")
#	set (COMPILE_FLAGS "${COMPILE_FLAGS} ${OPENMP_C_FLAGS}")
else (OPENMP_FOUND)
	message(SEND_ERROR "Your GCC does not support parallel API OpenMP.")
endif(OPENMP_FOUND)

#configure_file (
#	"src/config.h.in"
#	"src/config.h"
#)

include_directories("./include")
add_library(gaalgo STATIC ./src/ga_algo.c)
target_link_libraries(gaalgo m)

add_executable(ga_example "test/ga_example.c")
target_link_libraries(ga_example LINK_PUBLIC gaalgo)

#find_package(Doxygen)
#if (DOXYGEN_FOUND)
#	doxygen_add_docs(doc)
#else (DOXYGEN_FOUND)
#	message(WARNING "Did not find Doxygen. Not building documentation")
#endif(DOXYGEN_FOUND)
