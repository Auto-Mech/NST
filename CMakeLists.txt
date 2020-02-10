project(NST Fortran)
cmake_minimum_required(VERSION 2.8)

find_package(LAPACK REQUIRED)

add_library(nstlibs
    ${PROJECT_SOURCE_DIR}/src/msxfreq.f
    ${PROJECT_SOURCE_DIR}/src/opt.f
    ${PROJECT_SOURCE_DIR}/src/proj.f
    ${PROJECT_SOURCE_DIR}/src/dd.f
    ${PROJECT_SOURCE_DIR}/src/airy.f)

message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

add_executable(nst.x
    ${PROJECT_SOURCE_DIR}/src/nst.f)

target_link_libraries(nst.x 
        nstlibs ${LAPACK_LIBRARIES})

install(TARGETS nst.x DESTINATION bin)
