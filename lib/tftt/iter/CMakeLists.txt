
set(tfttiter_SRCS
    leaves.h 
    leaves.cc 
    curve.h 
    curve.cc 
    activecurve.h 
    activecurve.cc 
    boundarycells.h 
    boundarycells.cc 
    leaforthos.h 
    leaforthos.cc 
    neighbours.h 
    neighbours.cc 
    cellrefiterator.h
    cellrefiterator.cc 
)

add_library(tfttiter STATIC
    ${tfttiter_SRCS}
)
