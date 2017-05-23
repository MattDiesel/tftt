
set(util_SRCS
    pars.cc 
)

# Declare the library
add_library(util STATIC
    ${util_SRCS}
)

# Specify here the include directories exported
# by this library
target_include_directories(util PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
