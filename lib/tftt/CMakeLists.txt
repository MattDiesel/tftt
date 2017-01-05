



set(TFTT_DATA_T                     "int")
set(TFTT_MODEL_T                    "usermodel")
set(TFTT_ID_DATA_T                  "uint32_t")
set(TFTT_DIM                        2)
set(TFTT_MAXDEPTH                   12)
set(TFTT_MINDEPTH                   2)
set(TFTT_BALANCE_MINCHANGE          20)
set(TFTT_BALANCE_CHECKINTERVAL      10)


configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
        "${PROJECT_SOURCE_DIR}/include/tftt/config.h"
    )
include_directories(${PROJECT_SOURCE_DIR}/include)


set(libtftt_SRCS
    treeid.h 
    treeid.cc
)

# Declare the library
add_library(libtftt STATIC
    ${libtftt_SRCS}
)

# Specify here the include directories exported
# by this library
target_include_directories(libtftt PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
