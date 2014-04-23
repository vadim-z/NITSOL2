set (src
	bratu.f
	fbratu.f
	hwscrtd.f
	jacvbratu.f
	psolb.f
	)

add_executable(bratu ${src})
target_link_libraries(bratu nitsol2)
target_link_libraries(bratu ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
### ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})

set_source_files_properties(${src} PROPERTIES COMPILE_FLAGS ${FCFLAGS_NOWARN})
set_target_properties(bratu PROPERTIES COMPILE_FLAGS "${FCFLAGS}")
set_target_properties(bratu PROPERTIES LINK_FLAGS "${LFLAGS}")
