# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xinlong/Documents/Github/PolynomialRootFinding

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xinlong/Documents/Github/PolynomialRootFinding/build

# Include any dependencies generated for this target.
include CMakeFiles/Tylor.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Tylor.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Tylor.dir/flags.make

CMakeFiles/Tylor.dir/src/tylortest.cpp.o: CMakeFiles/Tylor.dir/flags.make
CMakeFiles/Tylor.dir/src/tylortest.cpp.o: ../src/tylortest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xinlong/Documents/Github/PolynomialRootFinding/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Tylor.dir/src/tylortest.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Tylor.dir/src/tylortest.cpp.o -c /home/xinlong/Documents/Github/PolynomialRootFinding/src/tylortest.cpp

CMakeFiles/Tylor.dir/src/tylortest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tylor.dir/src/tylortest.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xinlong/Documents/Github/PolynomialRootFinding/src/tylortest.cpp > CMakeFiles/Tylor.dir/src/tylortest.cpp.i

CMakeFiles/Tylor.dir/src/tylortest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tylor.dir/src/tylortest.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xinlong/Documents/Github/PolynomialRootFinding/src/tylortest.cpp -o CMakeFiles/Tylor.dir/src/tylortest.cpp.s

# Object files for target Tylor
Tylor_OBJECTS = \
"CMakeFiles/Tylor.dir/src/tylortest.cpp.o"

# External object files for target Tylor
Tylor_EXTERNAL_OBJECTS =

Tylor: CMakeFiles/Tylor.dir/src/tylortest.cpp.o
Tylor: CMakeFiles/Tylor.dir/build.make
Tylor: CMakeFiles/Tylor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xinlong/Documents/Github/PolynomialRootFinding/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Tylor"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Tylor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Tylor.dir/build: Tylor

.PHONY : CMakeFiles/Tylor.dir/build

CMakeFiles/Tylor.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Tylor.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Tylor.dir/clean

CMakeFiles/Tylor.dir/depend:
	cd /home/xinlong/Documents/Github/PolynomialRootFinding/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xinlong/Documents/Github/PolynomialRootFinding /home/xinlong/Documents/Github/PolynomialRootFinding /home/xinlong/Documents/Github/PolynomialRootFinding/build /home/xinlong/Documents/Github/PolynomialRootFinding/build /home/xinlong/Documents/Github/PolynomialRootFinding/build/CMakeFiles/Tylor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Tylor.dir/depend
