project(DelDec)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

find_package(boost)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/local/include)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -O3 -F/Library/Frameworks -fsanitize-undefined-trap-on-error -fsanitize=return")

set(CMAKE_EXE_LINKER_FLAGS "-framework OpenCL")

set(LINK_FLAGS ${LINK_FLAGS} "-L/opt/local/lib/ -lboost_unit_test_framework-mt -F/Library/Frameworks -framework OpenCL")

add_definitions(-DBOOST_TEST_DYN_LINK)

set(del_src
    main.cpp)
add_executable(delaunay ${del_src})

set(test_src
    ball_test.h
    point_test.h
    subspan_test.h
    test.cpp)
add_executable(test ${test_src})

set_target_properties(test PROPERTIES
    LINK_FLAGS "-L/opt/local/lib/ -lboost_unit_test_framework-mt -F/Library/Frameworks -framework OpenCL")

# let Qt Creator show all the sources
file(GLOB dummy_src "*.cpp")
file(GLOB dummy_hdr "*.h")
file(GLOB dummy_krn "*.cl")

add_custom_target(dummy
    WORKING_DIRECTORY .
    SOURCES ${dummy_hdr} ${dummy_src} ${dummy_krn})


