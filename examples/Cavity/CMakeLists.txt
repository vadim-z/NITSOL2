set (src
	bihsl.f
	cavity.f
	dbihar.f
	dppfa.f
	dppsl.f
	dspfa.f
	dspsl.f
	fcav.f
	init2.f
	init.f
	jacvcav.f
	psol.f
	resid.f
	wload.f
	)

add_executable(cavity ${src})
target_link_libraries(cavity nitsol2)
target_link_libraries(cavity ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
### ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})

set_source_files_properties(${src} PROPERTIES COMPILE_FLAGS ${FCFLAGS_NOWARN})
set_target_properties(cavity PROPERTIES COMPILE_FLAGS "${FCFLAGS}")
set_target_properties(cavity PROPERTIES LINK_FLAGS "${LFLAGS}")
