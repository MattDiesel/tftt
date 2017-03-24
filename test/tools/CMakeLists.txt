
set(tool_sources
    plottr.cc
)
add_executable(plottr ${tool_sources})
target_link_libraries(plottr libtftt util ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})

set(pos2id_sources
    pos2id.cc
)
add_executable(pos2id ${pos2id_sources})
target_link_libraries(pos2id libtftt util ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
