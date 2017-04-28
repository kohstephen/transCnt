Getting Started
How to compile this library
To use this library, you need to install a few packages:

	sudo apt update && sudo apt install libboost-all-dev gnuplot -y

Assume you have the all the code in transcnt directory and g++ 6.2+ installed:
	
	cd transcnt
	git clone git@github.com:dstahlke/gnuplot-iostream.git
	make all

How to run test program
To add the shared library to your PATH:

	export LD_LIBRARY_PATH=$(pwd):LD_LIBRARY_PATH

To run the test program:

	make test
	./test

How to use this library in user program
Include #include “calculation.h” in your files.

And compile with:

$(CXX) $(CPPFLAGS) -L. -ltranscnt app.cpp -o app \ -lboost_iostreams -lboost_system -lboost_filesystem

Or customize the Makefile.
