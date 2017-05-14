

set(circle_sources
    main.cc
)
add_executable(circle ${circle_sources})
target_link_libraries(circle libtftt util ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
