setup_test(c_int_defect_srh test_c_defect_srh.c MPI LABELS unit)

# Copy HR.dat/HI.dat from c_int_elph_deph for the end-to-end solve test
if(NOT "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  add_custom_target(
    c_int_defect_srh-elph-data
    COMMENT "Copying HR.dat/HI.dat for c_int_defect_srh end-to-end test."
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_SOURCE_DIR}/tests/c_int_elph_deph/HR.dat
      ${CMAKE_CURRENT_BINARY_DIR}/HR.dat
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_SOURCE_DIR}/tests/c_int_elph_deph/HI.dat
      ${CMAKE_CURRENT_BINARY_DIR}/HI.dat
  )
  add_dependencies(c_int_defect_srh c_int_defect_srh-elph-data)
endif()
