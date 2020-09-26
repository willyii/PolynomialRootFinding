cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

Project(Polinomial)

# Add include path
include_directories("include")

# Main function source file
file(GLOB SOURCES "src/*.cpp")

# Build executable
add_executable(Test ${SOURCES})




