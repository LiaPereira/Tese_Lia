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


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp

# Include any dependencies generated for this target.
include CMakeFiles/cmTC_9cc7d.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cmTC_9cc7d.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cmTC_9cc7d.dir/flags.make

CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.o: CMakeFiles/cmTC_9cc7d.dir/flags.make
CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.o: testCXXCompiler.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.o -c /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp/testCXXCompiler.cxx

CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.i: cmake_force
	@echo "Preprocessing CXX source to CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp/testCXXCompiler.cxx > CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.i

CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.s: cmake_force
	@echo "Compiling CXX source to assembly CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp/testCXXCompiler.cxx -o CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.s

# Object files for target cmTC_9cc7d
cmTC_9cc7d_OBJECTS = \
"CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.o"

# External object files for target cmTC_9cc7d
cmTC_9cc7d_EXTERNAL_OBJECTS =

cmTC_9cc7d: CMakeFiles/cmTC_9cc7d.dir/testCXXCompiler.cxx.o
cmTC_9cc7d: CMakeFiles/cmTC_9cc7d.dir/build.make
cmTC_9cc7d: CMakeFiles/cmTC_9cc7d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cmTC_9cc7d"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cmTC_9cc7d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cmTC_9cc7d.dir/build: cmTC_9cc7d

.PHONY : CMakeFiles/cmTC_9cc7d.dir/build

CMakeFiles/cmTC_9cc7d.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cmTC_9cc7d.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cmTC_9cc7d.dir/clean

CMakeFiles/cmTC_9cc7d.dir/depend:
	cd /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp /home/ldcp/geant4-10.7.2/examples/pdb4proteinas_only/CMakeFiles/CMakeTmp/CMakeFiles/cmTC_9cc7d.dir/DependInfo.cmake
.PHONY : CMakeFiles/cmTC_9cc7d.dir/depend
