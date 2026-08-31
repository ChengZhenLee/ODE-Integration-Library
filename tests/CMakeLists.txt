add_executable(test_euler test_euler.cpp)
target_link_libraries(test_euler PRIVATE odelib)

add_executable(test_rk4 test_rk4.cpp)
target_link_libraries(test_rk4 PRIVATE odelib)

add_test(NAME test_euler COMMAND test_euler)
add_test(NAME test_rk4 COMMAND test_rk4)