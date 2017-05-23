
set(tfttiter_SRCS
    leaves.cc 
    curve.cc 
    activecurve.cc 
    boundarycells.cc 
    leaforthos.cc 
    neighbours.cc 
    cellrefiterator.cc 
)

add_library(tfttiter STATIC
    ${tfttiter_SRCS}
)
