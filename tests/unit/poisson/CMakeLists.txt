# Poisson Unit Tests CTest Configuration
# Register CLI-based tests for Poisson equation validation

# Label for filtering Poisson tests
set(POISSON_TEST_LABEL "poisson")

# T012: Smoke test - basic functionality check
add_test(
    NAME poisson_smoke_poly2
    COMMAND poisson_tests_cli --grid 8 8 --pair poly2
)
set_tests_properties(poisson_smoke_poly2 PROPERTIES 
    LABELS "${POISSON_TEST_LABEL}"
    TIMEOUT 30
)

# T013-T014: User Story 1 tests
add_test(
    NAME poisson_smoke_sinpi
    COMMAND poisson_tests_cli --grid 16 16 --pair sinpi
)
set_tests_properties(poisson_smoke_sinpi PROPERTIES 
    LABELS "${POISSON_TEST_LABEL}"
    TIMEOUT 30
)

add_test(
    NAME poisson_multiple_pairs
    COMMAND poisson_tests_cli --grid 8 8 --pair poly2 --pair sinpi --pair cospi
)
set_tests_properties(poisson_multiple_pairs PROPERTIES 
    LABELS "${POISSON_TEST_LABEL}"
    TIMEOUT 60
)

# T017: Shell script smoke test
add_test(
    NAME poisson_cli_smoke_script
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson_cli.sh $<TARGET_FILE:poisson_tests_cli>
)
set_tests_properties(poisson_cli_smoke_script PROPERTIES 
    LABELS "${POISSON_TEST_LABEL}"
    TIMEOUT 90
)

# ============================================================================
# DMStag Tests - Cell-centered discretization
# ============================================================================

add_test(
    NAME poisson_dmstag_smoke_poly2
    COMMAND poisson_dmstag_tests_cli --grid 8 8 --pair poly2
)
set_tests_properties(poisson_dmstag_smoke_poly2 PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};dmstag"
    TIMEOUT 30
)

add_test(
    NAME poisson_dmstag_smoke_sinpi
    COMMAND poisson_dmstag_tests_cli --grid 16 16 --pair sinpi
)
set_tests_properties(poisson_dmstag_smoke_sinpi PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};dmstag"
    TIMEOUT 30
)

add_test(
    NAME poisson_dmstag_multiple_pairs
    COMMAND poisson_dmstag_tests_cli --grid 8 8 --pair poly2 --pair sinpi --pair cospi
)
set_tests_properties(poisson_dmstag_multiple_pairs PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};dmstag"
    TIMEOUT 60
)

# DMStag shell script smoke test
add_test(
    NAME poisson_dmstag_cli_smoke_script
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_poisson_dmstag_cli.sh $<TARGET_FILE:poisson_dmstag_tests_cli>
)
set_tests_properties(poisson_dmstag_cli_smoke_script PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};dmstag"
    TIMEOUT 90
)

# ============================================================================
# Convergence Tests - Verify spatial convergence order
# ============================================================================

# T020: Poisson DMStag convergence test (verifies 2nd-order accuracy)
add_test(
    NAME poisson_dmstag_convergence
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/poisson_convergence.sh $<TARGET_FILE:poisson_DMStag_2D>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(poisson_dmstag_convergence PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};dmstag;convergence"
    TIMEOUT 300
    PASS_REGULAR_EXPRESSION "Expected spatial convergence rate: ~2.0"
)

# Vertex-centered Poisson convergence test (verifies 2nd-order accuracy)
add_test(
    NAME poisson_vertex_convergence
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/poisson_vertex_convergence.sh $<TARGET_FILE:poisson_vertex_2D>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(poisson_vertex_convergence PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};vertex;convergence"
    TIMEOUT 300
    PASS_REGULAR_EXPRESSION "Expected spatial convergence rate: ~2.0"
)

# Staggered (mixed) Poisson convergence test
add_test(
    NAME poisson_staggered_convergence
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/poisson_staggered_convergence.sh $<TARGET_FILE:poisson_staggered_2D>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(poisson_staggered_convergence PROPERTIES 
    LABELS "${POISSON_TEST_LABEL};staggered;convergence"
    TIMEOUT 300
    PASS_REGULAR_EXPRESSION "Expected spatial convergence rate: ~2.0"
)

# Additional tests will be added in subsequent phases (US2-US3)
