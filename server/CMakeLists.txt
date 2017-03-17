project(persistence)

list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")	# enables C++ 11 features
set(CMAKE_BUILD_TYPE Release)

find_package(Boost 1.36.0 REQUIRED serialization)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(${persistence_SOURCE_DIR})

add_executable(toposerver toposerver.cc mongoose.c)
target_link_libraries( toposerver topology geometry ANN ${Boost_LIBRARIES} )

add_executable(test test.cc)
target_link_libraries( test topology geometry ANN ${Boost_LIBRARIES} )

subdirs(geometry topology)
