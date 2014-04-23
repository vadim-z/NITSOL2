set (src
	fpm.f
	getrow.f
	ilutf_aux.f
	ilutf.f
	jacvpm.f
	pormed.f
	)

add_executable(pormed ${src})
target_link_libraries(pormed nitsol2)
target_link_libraries(pormed ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
### ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})

set_source_files_properties(${src} PROPERTIES COMPILE_FLAGS ${FCFLAGS_NOWARN})
set_target_properties(pormed PROPERTIES COMPILE_FLAGS "${FCFLAGS}")
set_target_properties(pormed PROPERTIES LINK_FLAGS "${LFLAGS}")
