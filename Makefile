CXX = g++
CPPFLAGS = --std=c++17 -g -o2

all:
	$(CXX) $(CPPFLAGS) -shared -o libtranscnt.so -fPIC *.cpp
	export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH
test:
	$(CXX) $(CPPFLAGS) -L. -ltranscnt test.cpp -o test -lboost_iostreams -lboost_system -lboost_filesystem
val:
	valgrind --leak-check=full --track-origins=yes ./test
clean:
	rm libtranscnt.so test
