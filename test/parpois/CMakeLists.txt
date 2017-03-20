

set(parpois_sources
    main.cc
)
add_executable(parpois ${parpois_sources})
target_link_libraries(parpois libtftt util ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
