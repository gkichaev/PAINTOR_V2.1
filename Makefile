CC = g++
OPTS = -std=c++0x -Wall -O3

curr = "$(PWD)"

all: PAINTOR

PAINTOR: main.cpp
	$(CC) $(OPTS) main.cpp -fopenmp -lm -lnlopt -I/$(curr)/include -I/$(curr)/eigen/Eigen -L/$(curr)/lib -o PAINTOR 

