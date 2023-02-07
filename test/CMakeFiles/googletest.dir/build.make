# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_SOURCE_DIR = /opt/apollo/contrib/hephaestus

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /opt/apollo/contrib/hephaestus

# Utility rule file for googletest.

# Include any custom commands dependencies for this target.
include test/CMakeFiles/googletest.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/googletest.dir/progress.make

test/CMakeFiles/googletest: test/CMakeFiles/googletest-complete

test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-install
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-mkdir
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-download
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-update
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-patch
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-configure
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-build
test/CMakeFiles/googletest-complete: test/googletest-prefix/src/googletest-stamp/googletest-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/CMakeFiles
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/CMakeFiles/googletest-complete
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-done

test/googletest-prefix/src/googletest-stamp/googletest-build: test/googletest-prefix/src/googletest-stamp/googletest-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing build step for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build && $(MAKE)
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-build

test/googletest-prefix/src/googletest-stamp/googletest-configure: test/googletest-prefix/tmp/googletest-cfgcmd.txt
test/googletest-prefix/src/googletest-stamp/googletest-configure: test/googletest-prefix/src/googletest-stamp/googletest-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Performing configure step for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build && /usr/local/bin/cmake "-GUnix Makefiles" /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-configure

test/googletest-prefix/src/googletest-stamp/googletest-download: test/googletest-prefix/src/googletest-stamp/googletest-gitinfo.txt
test/googletest-prefix/src/googletest-stamp/googletest-download: test/googletest-prefix/src/googletest-stamp/googletest-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (git clone) for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src && /usr/local/bin/cmake -P /opt/apollo/contrib/hephaestus/test/googletest-prefix/tmp/googletest-gitclone.cmake
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-download

test/googletest-prefix/src/googletest-stamp/googletest-install: test/googletest-prefix/src/googletest-stamp/googletest-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "No install step for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build && /usr/local/bin/cmake -E echo_append
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-install

test/googletest-prefix/src/googletest-stamp/googletest-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-build
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix/tmp
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix/src
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E make_directory /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-mkdir

test/googletest-prefix/src/googletest-stamp/googletest-patch: test/googletest-prefix/src/googletest-stamp/googletest-update
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E echo_append
	cd /opt/apollo/contrib/hephaestus/test && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-patch

test/googletest-prefix/src/googletest-stamp/googletest-update: test/googletest-prefix/src/googletest-stamp/googletest-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/opt/apollo/contrib/hephaestus/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "No update step for 'googletest'"
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest && /usr/local/bin/cmake -E echo_append
	cd /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest && /usr/local/bin/cmake -E touch /opt/apollo/contrib/hephaestus/test/googletest-prefix/src/googletest-stamp/googletest-update

googletest: test/CMakeFiles/googletest
googletest: test/CMakeFiles/googletest-complete
googletest: test/googletest-prefix/src/googletest-stamp/googletest-build
googletest: test/googletest-prefix/src/googletest-stamp/googletest-configure
googletest: test/googletest-prefix/src/googletest-stamp/googletest-download
googletest: test/googletest-prefix/src/googletest-stamp/googletest-install
googletest: test/googletest-prefix/src/googletest-stamp/googletest-mkdir
googletest: test/googletest-prefix/src/googletest-stamp/googletest-patch
googletest: test/googletest-prefix/src/googletest-stamp/googletest-update
googletest: test/CMakeFiles/googletest.dir/build.make
.PHONY : googletest

# Rule to build all files generated by this target.
test/CMakeFiles/googletest.dir/build: googletest
.PHONY : test/CMakeFiles/googletest.dir/build

test/CMakeFiles/googletest.dir/clean:
	cd /opt/apollo/contrib/hephaestus/test && $(CMAKE_COMMAND) -P CMakeFiles/googletest.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/googletest.dir/clean

test/CMakeFiles/googletest.dir/depend:
	cd /opt/apollo/contrib/hephaestus && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/apollo/contrib/hephaestus /opt/apollo/contrib/hephaestus/test /opt/apollo/contrib/hephaestus /opt/apollo/contrib/hephaestus/test /opt/apollo/contrib/hephaestus/test/CMakeFiles/googletest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/googletest.dir/depend
