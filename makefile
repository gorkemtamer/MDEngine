CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -Ilib

# Layout:
#   src/     programs (files with main)
#   lib/     everything else, compiled into every program
#   bin/     executables (created by make)
#   input/   input .xyz files
#   output/  generated files (created by make)
#
# Absolute folder paths are baked in (see lib/paths.h), so the programs
# work whether you start them from the project root or from bin/.
DIRS = -DINPUT_DIR='"$(CURDIR)/input/"' -DOUTPUT_DIR='"$(CURDIR)/output/"'

# every .c/.h in lib/ is picked up automatically
# (makefile is a prerequisite too: changing flags or paths here rebuilds everything)
LIB_SRC = $(wildcard lib/*.cpp)
LIB_H   = $(wildcard lib/*.h)

all: bin/main bin/viewer bin/xyzio_test

bin/main: src/main.cpp $(LIB_SRC) $(LIB_H) makefile | bin output
	$(CXX) $(CXXFLAGS) $(DIRS) src/main.cpp $(LIB_SRC) -o $@ -lm

bin/viewer: src/viewer.cpp $(LIB_SRC) $(LIB_H) makefile | bin output
	$(CXX) $(CXXFLAGS) $(DIRS) src/viewer.cpp $(LIB_SRC) -o $@ -lraylib -lm

bin/xyzio_test: $(LIB_SRC) $(LIB_H) makefile | bin output
	$(CXX) $(CXXFLAGS) $(DIRS) -DXYZIO_TEST $(LIB_SRC) -o $@ -lm

bin output:
	mkdir -p $@

# short names: "make viewer" etc.
main viewer xyzio_test: %: bin/%

# removes the executables only; output/ is kept (it may hold results)
clean:
	rm -rf bin

.PHONY: all clean main viewer xyzio_test
