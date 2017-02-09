

set(circle_sources
    main.cc
    formatstring.h 
)
add_executable(circle ${circle_sources})
target_link_libraries(circle libtftt)
