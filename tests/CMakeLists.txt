add_executable(test_euler test_euler.cpp)
target_link_libraries(test_euler PRIVATE odelib)

add_test(NAME test_euler COMMAND test_euler)