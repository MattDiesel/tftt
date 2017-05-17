

set(circle2_sources
    main.cc
)
add_executable(circle2 ${circle2_sources})
target_link_libraries(circle2 libtftt util ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
