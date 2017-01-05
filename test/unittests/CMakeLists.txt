




set(test_sources
	main.cc
)
add_executable(alltests ${test_sources})
target_link_libraries(alltests libgtest libtftt)
