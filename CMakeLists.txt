cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
Project(Polinomial)


file(GLOB project_SRCS src/tmp.cpp src/util.cpp)

add_executable( Poly ${project_SRCS} )


