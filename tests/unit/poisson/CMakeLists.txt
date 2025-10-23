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

# Additional tests will be added in subsequent phases (US2-US3)
