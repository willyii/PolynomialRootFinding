cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
Project(Polinomial)

# include path
include_directories(include)

# set library
file( GLOB LIB_SOURCES src/*.cpp )
file( GLOB LIB_HEADERS include/*.h )
add_library( POLYLIB ${LIB_SOURCES} ${LIB_HEADERS} )

add_executable( Poly main/main.cpp )
# Make sure YourLib is linked to each app
target_link_libraries( Poly  POLYLIB )

add_executable( Test test/test.cpp )
# Make sure YourLib is linked to each app
target_link_libraries( Test POLYLIB )

set(CMAKE_CXX_STANDARD_REQUIRED YES)  
find_package(GSL REQUIRED)   
target_link_libraries(Test GSL::gsl)


