
set(tfttio_SRCS
    tfttio.cc
    serialize.cc
)

add_library(tfttio STATIC
    ${tfttio_SRCS}
)

target_link_libraries(tfttio libtftt)
