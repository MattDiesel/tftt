

set(parcircle_sources
    main.cc
)
add_executable(parcircle ${parcircle_sources})
target_link_libraries(parcircle libtftt util ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
