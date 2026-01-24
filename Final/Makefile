CC = g++
CFlags = -std=c++17 -O3

all: main

main: main.o
	$(CC) $(CFlags) -o Staple_Insertion main.o

main.o: main.cpp
	$(CC) $(CFlags) -c main.cpp

clean:
	rm -rf main.o
	rm -rf Staple_Insertion