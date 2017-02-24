

set(mbrot_sources
	main.cc 
)
add_executable(mandelbrot ${mbrot_sources})
target_link_libraries(mandelbrot libtftt util)
