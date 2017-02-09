project(persistence)

list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")	# enables C++ 11 features
set(CMAKE_BUILD_TYPE Release)

include_directories(${persistence_SOURCE_DIR})

add_executable(toposerver toposerver.cc mongoose.c)
target_link_libraries(toposerver topology geometry ANN)

subdirs(geometry topology)
