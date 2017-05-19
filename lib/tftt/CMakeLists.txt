
set(libtftt_SRCS
    tftt.h
    tftt.cc
    treeid.h 
    treeid.cc
    cellref.h
    cellref.cc
    gray.h
    gray.cc
    hilbert.h 
    hilbert.cc 
    faceref.h 
    faceref.cc 
    vertexref.h 
    vertexref.cc 
    tfttops.h 
    tfttops.cc 
    parallel.h 
    parallel.cc 
    adapt.h 
    adapt.cc 
    fttcore.h 
    fttcore.cc 
    config.h 
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
