all: main

main: main.o ACO.o EVRP.o stats.o
	g++ -std=c++11 -static-libgcc -static-libstdc++ main.o ACO.o EVRP.o stats.o -o main

main.o: main.cpp
	g++ -c main.cpp

ACO.o: ACO.cpp
	g++ -c ACO.cpp

DBGP.o: EVRP.cpp
	g++ -c DBGP.cpp

stats.o: stats.cpp
	g++ -c stats.cpp

clean: 
	rm -rt *o main
