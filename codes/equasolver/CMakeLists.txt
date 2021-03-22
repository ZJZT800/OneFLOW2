cmake_minimum_required (VERSION 2.8)

SET (CURRENT_DIR ${PROJECT_SOURCE_DIR})
SET (CMAKE_CXX_COMPILER g++)
SET (CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS}  -std=c++11 -Wall -g -O0")
SET (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
SET (CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS}  -std=c++11 -O3 ")
SET(CMAKE_LIBRARY_PATH /opt/intel/mkl/lib/intel64 /usr/lcoal/lib)

INCLUDE_DIRECTORIES(/opt/intel/mkl/include /usr/local/include)

AUX_SOURCE_DIRECTORY(${CURRENT_DIR} SRC  )
ADD_EXECUTABLE(GMRes ${SRC})


FIND_LIBRARY(lib_mkl_rt			mkl_rt		)

TARGET_LINK_LIBRARIES(GMRes 	 ${lib_mkl_rt} )
