function(add_benchmark name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE odelib)
    target_compile_definitions(${name} PRIVATE PROJECT_ROOT_DIR="${CMAKE_SOURCE_DIR}")
endfunction()


add_benchmark(fixed_step_convergence)
add_benchmark(adaptive_step_trace)
add_benchmark(adaptive_step_tolerance)