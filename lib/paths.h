#ifndef PATHS_H
#define PATHS_H

// ---------------------------------------------------------------------
// Project folders. The makefile passes absolute paths
// (-DINPUT_DIR=".../input/" -DOUTPUT_DIR=".../output/"), so the programs
// find their files no matter which directory they are started from.
// The fallbacks below are only used when compiling without the makefile.
//
// Use with string-literal concatenation:  INPUT_DIR "test.xyz"
// ---------------------------------------------------------------------
#ifndef INPUT_DIR
#define INPUT_DIR "input/"
#endif

#ifndef OUTPUT_DIR
#define OUTPUT_DIR "output/"
#endif

#endif
