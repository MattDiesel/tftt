cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(LinkedTree)

if($ENV{BUILD_NUMBER})
  set(BUILD_NUMBER $ENV{BUILD_NUMBER})
else()
  set(BUILD_NUMBER 0)
endif()

set(TFTT_NAME "LinkedTree")
set(TFTT_VERSION "0.0.${BUILD_NUMBER}")
set(TFTT_BRIEF "${TFTT_NAME} implements a parallel Threaded Fully-Threaded-Tree structure.")

set(COMMON_FLAGS "-std=c++11")
set(DEBUG_FLAGS "${COMMON_FLAGS} -g -O0 -Wall")
set(RELEASE_FLAGS "${COMMON_FLAGS} -O3")

set(CMAKE_CXX_FLAGS ${RELEASE_FLAGS})
set(CMAKE_C_FLAGS ${RELEASE_FLAGS})

set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${GCOV_FLAGS})
set(CMAKE_C_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_C_FLAGS_DEBUG ${GCOV_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS_DEBUG ${GCOV_FLAGS})

set(CMAKE_CONFIGURATION_TYPES Debug Release)



include_directories(./include)

add_subdirectory(lib)

include_directories(./include)
add_subdirectory(test)
