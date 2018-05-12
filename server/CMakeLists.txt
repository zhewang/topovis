project(persistence)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -pthread)
set(CMAKE_BUILD_TYPE Release)

include_directories(${persistence_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/ann_1.1.2/include)
link_directories(${CMAKE_SOURCE_DIR}/ann_1.1.2/lib)

subdirs(geometry topology)

add_executable(toposerver toposerver.cc mongoose.c topocubes.cc)
target_link_libraries( toposerver topology geometry ANN )

add_executable(example_gen example_gen.cc)
target_link_libraries( example_gen topology geometry ANN )

add_executable(preprocess preprocessing.cc)
target_link_libraries( preprocess topology geometry ANN )
