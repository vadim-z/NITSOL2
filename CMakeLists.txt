cmake_minimum_required(VERSION 2.6)
project(NITSOL2 Fortran)

if(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	set (LFLAGS "")
	set (BLAS_LIBRARIES "")
	set (LAPACK_LIBRARIES "")
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		set (FCFLAGS "/Qmkl /MT /fp:precise /QxHOST /O2 /warn:all") 
	else()
		set (FCFLAGS "/Qmkl /MT /fp:precise /Zi /Od /check /traceback /warn:all")
	endif()
	set (FCFLAGS_NOWARN "/warn:none")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	set (LFLAGS "-mkl -static-intel")
	set (BLAS_LIBRARIES "")
	set (LAPACK_LIBRARIES "")
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		set (FCFLAGS "-fp-model precise -xHOST -O2 -warn") 
	else()
		set (FCFLAGS "-fp-model precise -g -O0 -C -traceback -warn")
	endif()
	set (FCFLAGS_NOWARN "-nowarn")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
	if(MINGW)
		set (LFLAGS "-static-libgcc -static-libgfortran")
	else()		
		set (LFLAGS "")
	endif()
	find_package(BLAS REQUIRED)
	find_package(LAPACK REQUIRED)
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		set (FCFLAGS "-O2 -march=native -mfpmath=sse -Wall")
	else()
		set (FCFLAGS "-g -O0 -fcheck=all -fbacktrace -Wall")
	endif()
	set (FCFLAGS_NOWARN "-w")
endif()

set (libsrc
	src/nitdflts.f
	src/nitjv.f
	src/nitgm.f
	src/nitstb.f
	src/nittfq.f
	src/nitbt.f
	src/nitdrv.f
	src/nitsol.f
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(nitsol2 STATIC ${libsrc})
#target_link_libraries(somheat ${LUA_STATIC_LIB} ${BLAS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})
set_target_properties(nitsol2 PROPERTIES COMPILE_FLAGS "${FCFLAGS}")
#set_target_properties(somheat PROPERTIES LINK_FLAGS "${LFLAGS}")

# Ad-hoc
set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_RELEASE "")

# Add subprojects
add_subdirectory(examples/Bratu)
add_subdirectory(examples/Pormed)
add_subdirectory(examples/Cavity)
