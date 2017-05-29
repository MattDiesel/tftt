
set(tfttio_SRCS
    plot.cc
    plot3d.cc
    pnm.cc
    serialize.cc
    tikz.cc
)

add_library(tfttio STATIC
    ${tfttio_SRCS}
)

target_link_libraries(tfttio libtftt)
