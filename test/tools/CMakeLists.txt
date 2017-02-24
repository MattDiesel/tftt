
set(tool_sources
    plottr.cc
)
add_executable(plottr ${tool_sources})
target_link_libraries(plottr libtftt)
