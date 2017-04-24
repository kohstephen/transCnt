CXX = g++
CPPFLAGS = --std=c++17 -g -o2

all:
	$(CXX) $(CPPFLAGS) -shared -o libtranscnt.so -fPIC *.cpp
spheretest:
	$(CXX) $(CPPFLAGS) -L. -ltranscnt spheretest.cpp -o spheretest
val:
	valgrind --leak-check=full --track-origins=yes ./spheretest
clean:
	rm spheretest
