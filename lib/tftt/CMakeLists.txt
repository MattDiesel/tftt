
set(libtftt_SRCS
    tftt.cc
    treeid.cc
    faceref.cc 
    vertexref.cc 
    tfttops.cc 
    parallel.cc 
    adapt.cc 
    fttcore.cc 
    config.cc 
    cellref_static.cc
)

add_subdirectory(iter)
add_subdirectory(structure)
add_subdirectory(io)

# Declare the library
add_library(libtftt STATIC
    ${libtftt_SRCS}
)
target_link_libraries(libtftt tfttio tfttstruct tfttiter ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})

target_include_directories(libtftt PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
