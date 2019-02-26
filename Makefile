CC=g++ 
CFLAGS = -std=c++11 -pthread -Wall

cpu : encoder.cpp main.cpp encoder.h defines.h
	g++ -pthread -std=c++11 main.cpp encoder.cpp -o enc

.PHONY: clean 
clean: 
	rm enc
