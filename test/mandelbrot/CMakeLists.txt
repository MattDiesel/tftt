

set(mbrot_sources
	main.cc 
    pars.cc 
    pars.h 
    formatstring.h 
)
add_executable(mandelbrot ${mbrot_sources})
target_link_libraries(mandelbrot libtftt)
