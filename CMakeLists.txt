cmake_minimum_required(VERSION 3.9)
set(CMAKE_C_STANDARD 11)

project(gmmem C)

add_definitions("-DDEBUG -DHAVE_INLINE")
file(GLOB SRCFILES *.c)
file(GLOB INCLUDES *.h)

add_executable(gmmem ${SRCFILES} ${INCLUDES})
target_link_libraries(gmmem gsl blas)
