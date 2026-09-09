function(add_benchmark name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE odelib)
    target_compile_definitions(${name} PRIVATE PROJECT_ROOT_DIR="${CMAKE_SOURCE_DIR}")
endfunction()


add_benchmark(fixed_step_convergence)
add_benchmark(adaptive_step_trace)
add_benchmark(adaptive_step_tolerance)
add_benchmark(strong_convergence)
add_benchmark(monte_carlo_n)
add_benchmark(thread_speedup)
add_benchmark(walk_comparison)
add_benchmark(heston_model_comparison)

find_program(PYTHON_EXECUTABLE
    NAMES python.exe python3 python
    PATHS "${CMAKE_SOURCE_DIR}/.venv/Scripts" "${CMAKE_SOURCE_DIR}/.venv/bin"
    NO_DEFAULT_PATH
)
if (NOT PYTHON_EXECUTABLE)
    find_program(PYTHON_EXECUTABLE NAMES python3 python)
endif()

add_custom_target(run_benchmarks
    COMMAND fixed_step_convergence
    COMMAND adaptive_step_trace
    COMMAND adaptive_step_tolerance
    COMMAND strong_convergence
    COMMAND monte_carlo_n
    COMMAND thread_speedup
    COMMAND walk_comparison
    COMMAND heston_model_comparison
    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/plots/plot_all.py"
    DEPENDS fixed_step_convergence adaptive_step_trace adaptive_step_tolerance
            strong_convergence monte_carlo_n thread_speedup walk_comparison
    COMMENT "Running all benchmarks and regenerating plots"
    VERBATIM
)