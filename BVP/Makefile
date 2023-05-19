 output: test.o Poisson.o
	g++ -std=c++11 -Wall test.o Poisson.o -o output

 test.o: test.cpp Poisson.cpp
	g++ -c test.cpp -std=c++11

 Poisson.o: Poisson.cpp Poisson.h
	g++ -c Poisson.cpp -std=c++11

run :
	./output

 clean:
	rm *.o output
