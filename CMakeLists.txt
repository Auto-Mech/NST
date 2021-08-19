project(NST Fortran)
cmake_minimum_required(VERSION 2.8)

find_package(LAPACK REQUIRED)
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

add_library(nstlibs
    ${PROJECT_SOURCE_DIR}/src/msxfreq.f
    ${PROJECT_SOURCE_DIR}/src/opt.f
    ${PROJECT_SOURCE_DIR}/src/proj.f
    ${PROJECT_SOURCE_DIR}/src/airy.f)

add_library(ddmolprolibs
    ${PROJECT_SOURCE_DIR}/src/dd-molpro.f)
add_library(ddgaussianlibs
    ${PROJECT_SOURCE_DIR}/src/dd-gaussian.f)

add_executable(nst-gaussian.x
    ${PROJECT_SOURCE_DIR}/src/nst.f)
target_link_libraries(nst-gaussian.x 
    ${LAPACK_LIBRARIES} nstlibs ddgaussianlibs)

add_executable(nst-molpro.x
    ${PROJECT_SOURCE_DIR}/src/nst.f)
target_link_libraries(nst-molpro.x 
    ${LAPACK_LIBRARIES} nstlibs ddmolprolibs)

install(TARGETS nst-gaussian.x DESTINATION bin)
install(TARGETS nst-molpro.x DESTINATION bin)
