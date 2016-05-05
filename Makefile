GXX:= g++
CFLAGS:= -O3
all: main

main: algs.o bool_func.o bool_vec.o main.o
	$(GXX) $(CFLAGS) algs.o bool_func.o bool_vec.o main.o -o main

algs.o: algs.cpp
	$(GXX) $(CFLAGS) -c algs.cpp

bool_func.o: bool_func.cpp
	$(GXX) $(CFLAGS) -c bool_func.cpp

bool_vec.o: bool_vec.cpp
	$(GXX) $(CFLAGS) -c bool_vec.cpp

main.o: main.cpp
	$(GXX) $(CFLAGS) -c main.cpp

clear:
	rm -f *.o main
