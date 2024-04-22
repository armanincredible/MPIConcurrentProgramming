#include <mpi.h>
#include <stdio.h>
#include <string.h>

const int kElem_count = 10000000;
const size_t kElem_size = sizeof(char);

int Start_job(int rank)
{
    printf("Node %d: Started...\n", rank);

    char one_elem_buffer[kElem_size] = {};

    if ( rank == 0 )
    {
        double start_time = MPI_Wtime();

        for ( int i = 0; i != kElem_count; ++i )
        {
            MPI_Ssend(one_elem_buffer, 1, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
        }

        double end_time = MPI_Wtime();

        printf("Elapsed time: %lf on %d transimitons\n", end_time - start_time, kElem_count);

    } else if ( rank == 1 )
    {
        MPI_Status status;
        for ( int i = 0; i != kElem_count; ++i )
        {
            MPI_Recv(one_elem_buffer, 1, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
    }

    return 0;
}

int main(int argc, char *argv[])
{
    MPI_Init( &argc , &argv);
    MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("size %u\n", size);
    if (size < 2)
    {
        printf("Necessary to have 2 processors\n");
        return -1;
    }
    Start_job(rank);

    MPI_Finalize();
    return 0;
}
