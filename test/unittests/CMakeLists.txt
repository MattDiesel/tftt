




set(test_sources
	main.cc
	treeid_ut.cc
	treecellref_ut.cc
	tftt_ut.cc
    treegroup_ut.cc
    gray_ut.cc
)
add_executable(alltests ${test_sources})
target_link_libraries(alltests libgtest libtftt)
