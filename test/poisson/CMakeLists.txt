

set(poisson_sources
    main.cc
)
add_executable(poisson ${poisson_sources})
target_link_libraries(poisson libtftt util)
