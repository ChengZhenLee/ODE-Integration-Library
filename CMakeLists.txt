cmake_minimum_required(VERSION 3.16)
project(odelib CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(odelib INTERFACE)
target_include_directories(odelib INTERFACE include)

add_executable(decay_euler examples/decay_euler.cpp)
target_link_libraries(decay_euler PRIVATE odelib)

add_executable(decay_rk4 examples/decay_rk4.cpp)
target_link_libraries(decay_rk4 PRIVATE odelib)

add_executable(scaled_decay_rk4 examples/scaled_decay_rk4.cpp)
target_link_libraries(scaled_decay_rk4 PRIVATE odelib)

add_executable(decay_rk4_vector examples/decay_rk4_vector.cpp)
target_link_libraries(decay_rk4_vector PRIVATE odelib)

enable_testing()
add_subdirectory(tests)