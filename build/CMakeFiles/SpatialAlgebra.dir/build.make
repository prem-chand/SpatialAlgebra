# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/premchand/Documents/GitHub/SpatialAlgebra

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/premchand/Documents/GitHub/SpatialAlgebra/build

# Include any dependencies generated for this target.
include CMakeFiles/SpatialAlgebra.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/SpatialAlgebra.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/SpatialAlgebra.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SpatialAlgebra.dir/flags.make

CMakeFiles/SpatialAlgebra.dir/codegen:
.PHONY : CMakeFiles/SpatialAlgebra.dir/codegen

CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ArticulatedBodyInertia.cpp
CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ArticulatedBodyInertia.cpp

CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ArticulatedBodyInertia.cpp > CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ArticulatedBodyInertia.cpp -o CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/AxialScrewTransform.cpp
CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/AxialScrewTransform.cpp

CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/AxialScrewTransform.cpp > CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/AxialScrewTransform.cpp -o CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ForceVector.cpp
CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ForceVector.cpp

CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ForceVector.cpp > CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/ForceVector.cpp -o CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/MotionVector.cpp
CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/MotionVector.cpp

CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/MotionVector.cpp > CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/MotionVector.cpp -o CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/PluckerTransform.cpp
CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/PluckerTransform.cpp

CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/PluckerTransform.cpp > CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/PluckerTransform.cpp -o CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/RigidBodyInertia.cpp
CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/RigidBodyInertia.cpp

CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/RigidBodyInertia.cpp > CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/RigidBodyInertia.cpp -o CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/Rotation.cpp
CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/Rotation.cpp

CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/Rotation.cpp > CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/Rotation.cpp -o CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialOperations.cpp
CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialOperations.cpp

CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialOperations.cpp > CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialOperations.cpp -o CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialVector.cpp
CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialVector.cpp

CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialVector.cpp > CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/SpatialVector.cpp -o CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.s

CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o: CMakeFiles/SpatialAlgebra.dir/flags.make
CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o: /Users/premchand/Documents/GitHub/SpatialAlgebra/src/main.cpp
CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o: CMakeFiles/SpatialAlgebra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o -MF CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o.d -o CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o -c /Users/premchand/Documents/GitHub/SpatialAlgebra/src/main.cpp

CMakeFiles/SpatialAlgebra.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SpatialAlgebra.dir/src/main.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/premchand/Documents/GitHub/SpatialAlgebra/src/main.cpp > CMakeFiles/SpatialAlgebra.dir/src/main.cpp.i

CMakeFiles/SpatialAlgebra.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SpatialAlgebra.dir/src/main.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/premchand/Documents/GitHub/SpatialAlgebra/src/main.cpp -o CMakeFiles/SpatialAlgebra.dir/src/main.cpp.s

# Object files for target SpatialAlgebra
SpatialAlgebra_OBJECTS = \
"CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o" \
"CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o"

# External object files for target SpatialAlgebra
SpatialAlgebra_EXTERNAL_OBJECTS =

libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/ArticulatedBodyInertia.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/AxialScrewTransform.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/ForceVector.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/MotionVector.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/PluckerTransform.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/RigidBodyInertia.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/Rotation.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/SpatialOperations.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/SpatialVector.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/src/main.cpp.o
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/build.make
libSpatialAlgebra.a: CMakeFiles/SpatialAlgebra.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX static library libSpatialAlgebra.a"
	$(CMAKE_COMMAND) -P CMakeFiles/SpatialAlgebra.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SpatialAlgebra.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SpatialAlgebra.dir/build: libSpatialAlgebra.a
.PHONY : CMakeFiles/SpatialAlgebra.dir/build

CMakeFiles/SpatialAlgebra.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SpatialAlgebra.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SpatialAlgebra.dir/clean

CMakeFiles/SpatialAlgebra.dir/depend:
	cd /Users/premchand/Documents/GitHub/SpatialAlgebra/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/premchand/Documents/GitHub/SpatialAlgebra /Users/premchand/Documents/GitHub/SpatialAlgebra /Users/premchand/Documents/GitHub/SpatialAlgebra/build /Users/premchand/Documents/GitHub/SpatialAlgebra/build /Users/premchand/Documents/GitHub/SpatialAlgebra/build/CMakeFiles/SpatialAlgebra.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/SpatialAlgebra.dir/depend

