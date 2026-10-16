cmake_minimum_required(VERSION 3.20)
project(MDEngine LANGUAGES CXX)

# carried from the old makefile
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  #-std=c++20
add_compile_options(-Wall -Wextra -O2)

#executables go to build/bin/ same structure as before;
#output/ must exist
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)

#library, everything in lib/, compiled once and linked
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_SOURCE_DIR}/lib/*.cpp
     ${CMAKE_SOURCE_DIR}/lib/potentials/*.cpp
     ${CMAKE_SOURCE_DIR}/lib/integrators/*.cpp
     ${CMAKE_SOURCE_DIR}/lib/thermostats/*.cpp)
add_library(mdengine STATIC ${LIB_SOURCES})
target_include_directories(mdengine PUBLIC ${CMAKE_SOURCE_DIR}/lib)
target_compile_definitions(mdengine PUBLIC 
    INPUT_DIR="${CMAKE_SOURCE_DIR}/input/"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/output/")

#programs
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE mdengine)
add_executable(nacl src/nacl.cpp)
target_link_libraries(nacl PRIVATE mdengine)

find_package(raylib QUIET)
if(raylib_FOUND)
    add_executable(viewer src/viewer.cpp)
    target_link_libraries(viewer PRIVATE mdengine raylib)
else()
    message(STATUS "raylib not found: the viewer is not built")
endif()

add_executable(xyzio_test tests/xyzio_test.cpp)
target_link_libraries(xyzio_test PRIVATE mdengine)

#tests
enable_testing()

function(add_reference_test name target reference)
    add_test(NAME ${name}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${target}>
                -DREFERENCE=${CMAKE_SOURCE_DIR}/reference/${reference}
                -P ${CMAKE_SOURCE_DIR}/tests/compare_output.cmake)
endfunction()

function(add_physics_test name)
    add_executable(${name} tests/${name}.cpp)
    target_link_libraries(${name} PRIVATE mdengine)
    add_test(NAME ${name} COMMAND ${name})
endfunction()


add_reference_test(main_reference main main_seed42.txt)
add_reference_test(xyzio_roundtrip xyzio_test xyzio_test.txt)

add_physics_test(rng_test)
add_physics_test(force_test)
add_physics_test(newton_test)
add_physics_test(energy_test)
add_physics_test(reversibility_test)
add_physics_test(trap_test)
add_physics_test(coulomb_test)
add_physics_test(thermostat_test)