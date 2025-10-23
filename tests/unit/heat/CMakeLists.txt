# Heat equation unit tests using CTest

# Convergence test for heat equation DMStag solver (spatial convergence)
add_test(
  NAME heat_dmstag_convergence
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/heat_convergence.sh $<TARGET_FILE:heat_DMStag_2D>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(heat_dmstag_convergence PROPERTIES
  LABELS "heat;dmstag;convergence;spatial"
  TIMEOUT 300
  PASS_REGULAR_EXPRESSION "Expected spatial convergence rate: ~2.0"
)

# Temporal convergence test for heat equation (time step convergence, Implicit Euler)
add_test(
  NAME heat_dmstag_temporal_convergence
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/heat_temporal_convergence.sh $<TARGET_FILE:heat_DMStag_2D>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(heat_dmstag_temporal_convergence PROPERTIES
  LABELS "heat;dmstag;convergence;temporal;euler"
  TIMEOUT 300
  PASS_REGULAR_EXPRESSION "Expected temporal convergence rate: ~1.0"
)

# Temporal convergence test for heat equation (Crank-Nicolson)
add_test(
  NAME heat_dmstag_CN_temporal_convergence
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/heat_CN_temporal_convergence.sh $<TARGET_FILE:heat_DMStag_CN_2D>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(heat_dmstag_CN_temporal_convergence PROPERTIES
  LABELS "heat;dmstag;convergence;temporal;crank-nicolson"
  TIMEOUT 3000
  PASS_REGULAR_EXPRESSION "Expected temporal convergence rate: ~2.0"
)
