# Stokes equation convergence tests

# Make convergence test script executable
file(CHMOD ${CMAKE_CURRENT_SOURCE_DIR}/stokes_convergence.sh
     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                 GROUP_READ GROUP_EXECUTE
                 WORLD_READ WORLD_EXECUTE)

# Stokes projection method convergence test
add_test(
    NAME stokes_projection_convergence
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/stokes_convergence.sh 
            ${CMAKE_BINARY_DIR}/stokes_DMStag_projection_2D
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(stokes_projection_convergence PROPERTIES
    LABELS "stokes;dmstag;convergence;projection;spatial"
    TIMEOUT 600
)
