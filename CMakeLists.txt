
# Se compila todo el proyecto usando gfortran
cmake_minimum_required(VERSION 3.13)

project(dbrown)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    set(dialect "-ffree-form -std=f2008 -fimplicit-none -fPIC -Ofast -march=native -funroll-loops")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(dialect "-fast -stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

# Se exige que LAPACK este instalado
find_package(LAPACK REQUIRED)

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compilar
#
# Esta linea junta a todos los archivos y al programa principal
file(GLOB_RECURSE sources include/*.f90 src/*.f90)
add_executable(dbrown ${sources})
# Aqui se añade la libreria de LAPACK y BLAS
target_link_libraries(dbrown ${LAPACK_LIBRARIES})