cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(T1 Fortran C)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffloat-store -Wall") # -fopenmp"
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fcheck=all -fbacktrace")
endif()

#set(MOD_DIR ${CMAKE_BINARY_DIR}/mod)

set(src
	t1.f90
	mmio.f
	../../src/nitstb.f
	../../src/nittfq.f
	../../src/nitbt.f
	../../src/nitdflts.f
	../../src/nitdrv.f
	#../../src/nitcgm.f
	../../src/nitgm.f
	../../src/nitjv.f
	../../src/nitsol.f
)

add_executable(T1 ${src})
target_link_libraries(T1 PRIVATE BLAS::BLAS LAPACK::LAPACK)
