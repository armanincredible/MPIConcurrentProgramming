#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

//Estimating the value of Pi using Monte Carlo algorithm
void estimate_pi_value(int argc, char* argv[])
{
    unsigned long long total_iterations = 100000000;
    int rank, size;

    size_t points_in_circle = 0;
    size_t total_points_in_circle = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Init generator of random numbers by different seed on each mpi process
    srand(time(NULL) + rank);
    
    int iterations_per_process = total_iterations / size;

    for(size_t i = 0; i < iterations_per_process; i++) {
        double x = (double)rand() / RAND_MAX;
        double y = (double)rand() / RAND_MAX;
        double distance_squared = x * x + y * y;
        // Here we suggest that RNG will not generate a lot of simple pair of numbers 
        if (distance_squared <= 1) 
        {
            points_in_circle++;
        }
    }

    MPI_Reduce(&points_in_circle, &total_points_in_circle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        double pi_estimate = 4.0 * total_points_in_circle / total_iterations;
        printf("Estimated pi value: %f\n", pi_estimate);
    }

    MPI_Finalize();
}

// Estimating the value of Pi using Monte Carlo
int main(int argc, char* argv[])
{    
    estimate_pi_value(argc, argv);
    return 0;
}