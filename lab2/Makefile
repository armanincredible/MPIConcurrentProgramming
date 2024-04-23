all: collective simple

collective: collective.cpp 
	g++ collective.cpp -lm -pthread -O2 -Wall -std=c++14 -o collective.out 

simple: simple.cpp
	g++ simple.cpp -lm -pthread -O2 -Wall -std=c++14 -o simple.out 

global_stack: global_stack.cpp
	g++ global_stack.cpp -lm -pthread -O2 -Wall -std=c++14 -o stack.out 


clean:
	rm -f *.out