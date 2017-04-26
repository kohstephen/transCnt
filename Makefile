CXX = g++-6
CPPFLAGS = --std=c++17 -g -o2

all:
	$(CXX) $(CPPFLAGS) -shared -o libtranscnt.so -fPIC *.cpp
test:
	$(CXX) $(CPPFLAGS) -L. -ltranscnt test.cpp -o test 
val:
	valgrind --leak-check=full --track-origins=yes ./test
clean:
	rm libtranscnt.so test
