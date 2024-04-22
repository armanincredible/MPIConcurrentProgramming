all: parallel sequence

parallel: parallel.cpp 
	mpic++ parallel.cpp solver.cpp -lraylib -o par.out

sequence: sequence.cpp
	mpic++ sequence.cpp solver.cpp -lraylib -o seq.out

clean:
	rm -f *.out