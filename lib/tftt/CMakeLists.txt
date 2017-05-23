
set(libtftt_SRCS
    tftt.cc
    treeid.cc
    cellref.cc
    gray.cc
    hilbert.cc 
    faceref.cc 
    vertexref.cc 
    tfttops.cc 
    parallel.cc 
    adapt.cc 
    fttcore.cc 
    config.cc 
)

add_subdirectory(iter)
add_subdirectory(structure)
add_subdirectory(io)

# Declare the library
add_library(libtftt STATIC
    ${libtftt_SRCS}
)
target_link_libraries(libtftt tfttio tfttstruct tfttiter)

target_include_directories(libtftt PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
