 output: test.o Multigrid.o
	g++ -std=c++11 -Wall test.o Multigrid.o -o output

 test.o: test.cpp Multigrid.cpp
	g++ -c test.cpp -std=c++11

 Multigrid.o: Multigrid.cpp Multigrid.h
	g++ -c Multigrid.cpp -std=c++11

run :
	./output

 clean:
	rm *.o output
