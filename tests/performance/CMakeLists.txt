# Performance Tests CTest Configuration
# Tests for parallel efficiency and scaling behavior

set(PERFORMANCE_TEST_LABEL "performance")

# ============================================================================
# Parallel Efficiency Tests
# ============================================================================

# Test parallel efficiency for Poisson solver
add_test(
    NAME parallel_efficiency_poisson
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/parallel_efficiency.sh 
            $<TARGET_FILE:unified_pde_solver_2D>
            poisson
            128
            128
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(parallel_efficiency_poisson PROPERTIES 
    LABELS "${PERFORMANCE_TEST_LABEL};parallel;poisson"
    TIMEOUT 600
    PASS_REGULAR_EXPRESSION "PASS: Parallel efficiency is acceptable"
)

# Test parallel efficiency for Heat solver (smaller grid, shorter time)
add_test(
    NAME parallel_efficiency_heat
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/quick_parallel_test.sh 
            $<TARGET_FILE:unified_pde_solver_2D>
            heat
            64
            64
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(parallel_efficiency_heat PROPERTIES 
    LABELS "${PERFORMANCE_TEST_LABEL};parallel;heat"
    TIMEOUT 300
    PASS_REGULAR_EXPRESSION "PASS: Parallel efficiency is acceptable"
)

# Test with larger grid for more realistic parallel behavior
add_test(
    NAME parallel_efficiency_poisson_large
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/parallel_efficiency.sh 
            $<TARGET_FILE:unified_pde_solver_2D>
            poisson
            256
            256
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(parallel_efficiency_poisson_large PROPERTIES 
    LABELS "${PERFORMANCE_TEST_LABEL};parallel;poisson;large"
    TIMEOUT 1200
    PASS_REGULAR_EXPRESSION "PASS: Parallel efficiency is acceptable"
)
