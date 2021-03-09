cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
Project(Polinomial)

include_directories("/usr/local/Cellar/boost/1.75.0_1/include")


add_executable(Poly src/main.cpp ) 
add_executable(Tylor src/tylortest.cpp)
