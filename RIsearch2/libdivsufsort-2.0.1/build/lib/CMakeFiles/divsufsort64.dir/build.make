# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build

# Include any dependencies generated for this target.
include lib/CMakeFiles/divsufsort64.dir/depend.make

# Include the progress variables for this target.
include lib/CMakeFiles/divsufsort64.dir/progress.make

# Include the compile flags for this target's objects.
include lib/CMakeFiles/divsufsort64.dir/flags.make

lib/CMakeFiles/divsufsort64.dir/divsufsort.o: lib/CMakeFiles/divsufsort64.dir/flags.make
lib/CMakeFiles/divsufsort64.dir/divsufsort.o: ../lib/divsufsort.c
	$(CMAKE_COMMAND) -E cmake_progress_report /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/CMakeFiles/divsufsort64.dir/divsufsort.o"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/divsufsort.o   -c /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/divsufsort.c

lib/CMakeFiles/divsufsort64.dir/divsufsort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/divsufsort.i"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/divsufsort.c > CMakeFiles/divsufsort64.dir/divsufsort.i

lib/CMakeFiles/divsufsort64.dir/divsufsort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/divsufsort.s"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/divsufsort.c -o CMakeFiles/divsufsort64.dir/divsufsort.s

lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires:
.PHONY : lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires

lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides: lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires
	$(MAKE) -f lib/CMakeFiles/divsufsort64.dir/build.make lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides.build
.PHONY : lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides

lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides.build: lib/CMakeFiles/divsufsort64.dir/divsufsort.o

lib/CMakeFiles/divsufsort64.dir/sssort.o: lib/CMakeFiles/divsufsort64.dir/flags.make
lib/CMakeFiles/divsufsort64.dir/sssort.o: ../lib/sssort.c
	$(CMAKE_COMMAND) -E cmake_progress_report /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/CMakeFiles/divsufsort64.dir/sssort.o"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/sssort.o   -c /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/sssort.c

lib/CMakeFiles/divsufsort64.dir/sssort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/sssort.i"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/sssort.c > CMakeFiles/divsufsort64.dir/sssort.i

lib/CMakeFiles/divsufsort64.dir/sssort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/sssort.s"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/sssort.c -o CMakeFiles/divsufsort64.dir/sssort.s

lib/CMakeFiles/divsufsort64.dir/sssort.o.requires:
.PHONY : lib/CMakeFiles/divsufsort64.dir/sssort.o.requires

lib/CMakeFiles/divsufsort64.dir/sssort.o.provides: lib/CMakeFiles/divsufsort64.dir/sssort.o.requires
	$(MAKE) -f lib/CMakeFiles/divsufsort64.dir/build.make lib/CMakeFiles/divsufsort64.dir/sssort.o.provides.build
.PHONY : lib/CMakeFiles/divsufsort64.dir/sssort.o.provides

lib/CMakeFiles/divsufsort64.dir/sssort.o.provides.build: lib/CMakeFiles/divsufsort64.dir/sssort.o

lib/CMakeFiles/divsufsort64.dir/trsort.o: lib/CMakeFiles/divsufsort64.dir/flags.make
lib/CMakeFiles/divsufsort64.dir/trsort.o: ../lib/trsort.c
	$(CMAKE_COMMAND) -E cmake_progress_report /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/CMakeFiles/divsufsort64.dir/trsort.o"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/trsort.o   -c /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/trsort.c

lib/CMakeFiles/divsufsort64.dir/trsort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/trsort.i"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/trsort.c > CMakeFiles/divsufsort64.dir/trsort.i

lib/CMakeFiles/divsufsort64.dir/trsort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/trsort.s"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/trsort.c -o CMakeFiles/divsufsort64.dir/trsort.s

lib/CMakeFiles/divsufsort64.dir/trsort.o.requires:
.PHONY : lib/CMakeFiles/divsufsort64.dir/trsort.o.requires

lib/CMakeFiles/divsufsort64.dir/trsort.o.provides: lib/CMakeFiles/divsufsort64.dir/trsort.o.requires
	$(MAKE) -f lib/CMakeFiles/divsufsort64.dir/build.make lib/CMakeFiles/divsufsort64.dir/trsort.o.provides.build
.PHONY : lib/CMakeFiles/divsufsort64.dir/trsort.o.provides

lib/CMakeFiles/divsufsort64.dir/trsort.o.provides.build: lib/CMakeFiles/divsufsort64.dir/trsort.o

lib/CMakeFiles/divsufsort64.dir/utils.o: lib/CMakeFiles/divsufsort64.dir/flags.make
lib/CMakeFiles/divsufsort64.dir/utils.o: ../lib/utils.c
	$(CMAKE_COMMAND) -E cmake_progress_report /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object lib/CMakeFiles/divsufsort64.dir/utils.o"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/utils.o   -c /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/utils.c

lib/CMakeFiles/divsufsort64.dir/utils.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/utils.i"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/utils.c > CMakeFiles/divsufsort64.dir/utils.i

lib/CMakeFiles/divsufsort64.dir/utils.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/utils.s"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib/utils.c -o CMakeFiles/divsufsort64.dir/utils.s

lib/CMakeFiles/divsufsort64.dir/utils.o.requires:
.PHONY : lib/CMakeFiles/divsufsort64.dir/utils.o.requires

lib/CMakeFiles/divsufsort64.dir/utils.o.provides: lib/CMakeFiles/divsufsort64.dir/utils.o.requires
	$(MAKE) -f lib/CMakeFiles/divsufsort64.dir/build.make lib/CMakeFiles/divsufsort64.dir/utils.o.provides.build
.PHONY : lib/CMakeFiles/divsufsort64.dir/utils.o.provides

lib/CMakeFiles/divsufsort64.dir/utils.o.provides.build: lib/CMakeFiles/divsufsort64.dir/utils.o

# Object files for target divsufsort64
divsufsort64_OBJECTS = \
"CMakeFiles/divsufsort64.dir/divsufsort.o" \
"CMakeFiles/divsufsort64.dir/sssort.o" \
"CMakeFiles/divsufsort64.dir/trsort.o" \
"CMakeFiles/divsufsort64.dir/utils.o"

# External object files for target divsufsort64
divsufsort64_EXTERNAL_OBJECTS =

lib/libdivsufsort64.a: lib/CMakeFiles/divsufsort64.dir/divsufsort.o
lib/libdivsufsort64.a: lib/CMakeFiles/divsufsort64.dir/sssort.o
lib/libdivsufsort64.a: lib/CMakeFiles/divsufsort64.dir/trsort.o
lib/libdivsufsort64.a: lib/CMakeFiles/divsufsort64.dir/utils.o
lib/libdivsufsort64.a: lib/CMakeFiles/divsufsort64.dir/build.make
lib/libdivsufsort64.a: lib/CMakeFiles/divsufsort64.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library libdivsufsort64.a"
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && $(CMAKE_COMMAND) -P CMakeFiles/divsufsort64.dir/cmake_clean_target.cmake
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/divsufsort64.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/CMakeFiles/divsufsort64.dir/build: lib/libdivsufsort64.a
.PHONY : lib/CMakeFiles/divsufsort64.dir/build

lib/CMakeFiles/divsufsort64.dir/requires: lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires
lib/CMakeFiles/divsufsort64.dir/requires: lib/CMakeFiles/divsufsort64.dir/sssort.o.requires
lib/CMakeFiles/divsufsort64.dir/requires: lib/CMakeFiles/divsufsort64.dir/trsort.o.requires
lib/CMakeFiles/divsufsort64.dir/requires: lib/CMakeFiles/divsufsort64.dir/utils.o.requires
.PHONY : lib/CMakeFiles/divsufsort64.dir/requires

lib/CMakeFiles/divsufsort64.dir/clean:
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib && $(CMAKE_COMMAND) -P CMakeFiles/divsufsort64.dir/cmake_clean.cmake
.PHONY : lib/CMakeFiles/divsufsort64.dir/clean

lib/CMakeFiles/divsufsort64.dir/depend:
	cd /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1 /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/lib /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib /var/www/sites/rth.dk/htdocs/devresources/risearch/RIsearch/RIsearch2/libdivsufsort-2.0.1/build/lib/CMakeFiles/divsufsort64.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/CMakeFiles/divsufsort64.dir/depend

