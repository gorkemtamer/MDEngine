# Runs PROGRAM and compares its std output with the file REFERENCE
# cmake -DPROGRAM=<exe> -DREFERENCE=<file> -P compare_output.cmake

execute_process(COMMAND ${PROGRAM}
                OUTPUT_VARIABLE actual
                ERROR_VARIABLE errors
                RESULT_VARIABLE exit_code)

if(NOT exit_code EQUAL 0)
    message(FATAL_ERROR "${PROGRAM} exited with code ${exit_code}\n${errors}")
endif()

file(READ ${REFERENCE} expected)

if(NOT actual STREQUAL expected)
    get_filename_component(ref_name ${REFERENCE} NAME)
    set(actual_file ${CMAKE_CURRENT_BINARY_DIR}/actual_${ref_name})
    file(WRITE ${actual_file} "${actual}")
    message(FATAL_ERROR "output differs from the reference\n"
                        " see the difference with:\n"
                        " diff ${REFERENCE} ${actual_file}")
endif()