#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#include "types.h"
#include "solver.h"

static double **dataAlloc(uint32_t pts_x, uint32_t pts_y) {
    double **data = (double **) calloc(pts_x, sizeof(double*));
    if (data == NULL) return NULL;
    
    double* all_data = (double *)calloc(pts_x * pts_y, sizeof(double));
    if (all_data == NULL) {
        printf("Error");
        free(data); 
        return NULL;
    }

    data[0] = all_data;
    for (uint32_t i = 1; i < pts_x; i++) {
        data[i] = data[0] + i * pts_y; 
    }

    return data;
}

static void calculateRow(data_xy *state, uint32_t row, uint32_t from, uint32_t to, uint32_t rank, uint32_t size) {
    for (uint32_t i = from; i < to; i++) {
        state->data[i][row] = calculateNewPoint(state, i, row - 1);
    }

    double send_left = 0.0, recv_left = 0.0;
    double send_right = 0.0, recv_right = 0.0;

    if (rank > 0) {
        send_left = state->data[from][row];
        MPI_Sendrecv(&send_left, 1, MPI_DOUBLE, rank - 1, 0,
                     &recv_left, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        state->data[from - 1][row] = recv_left;
    }

    if (rank < size - 1) {
        send_right = state->data[to - 1][row];
        MPI_Sendrecv(&send_right, 1, MPI_DOUBLE, rank + 1, 0,
                     &recv_right, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        state->data[to][row] = recv_right;
    }
}


static void dataFree(data_xy *state) {
    free(state->data[0]);
    free(state->data);
}

int main(int argc, char **argv) {
    FILE* output = fopen("output.txt", "w");
    if (output == NULL) {
        fprintf(stderr, "error in opening output file\n");
    }

    int process_rank, size_of_cluster;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    dim x = {100.0, 100000, 0.0};
    x.step = x.len / x.pts_num;

    dim t = {0.01, 2000, 0.0};
    t.step = t.len / t.pts_num;

    const dimensions dims = {x, t};

    double **data = dataAlloc(x.pts_num, t.pts_num);

    initStartCoordDim(data, x.pts_num, x.step);

    initStartTimeDim(data, t.pts_num, t.step);

    initConditions(data, t.pts_num, x.pts_num, t.step);

    data_xy state = {data, dims, 1e-2};

    MPI_Barrier(MPI_COMM_WORLD);

    uint32_t chunk_size = (x.pts_num / (size_of_cluster));
    
    uint32_t init_i = 1 + chunk_size * process_rank;
    uint32_t end_i = init_i + chunk_size;

    if ((process_rank == size_of_cluster - 1) || (end_i > x.pts_num - 1)) {
        end_i = x.pts_num - 1;
    }

    printf("Process %d started [%d; %d]\n", process_rank, init_i, end_i - 1);

    double start = MPI_Wtime();

    for (uint32_t row = 1; row < t.pts_num; row++) {
        calculateRow(&state, row, init_i, end_i, process_rank, size_of_cluster);
    }
    double end = MPI_Wtime();
    printf("Process %d finished, time: %lf\n", process_rank, end-start);

    MPI_Barrier(MPI_COMM_WORLD);

    for (int rank = 0; rank < size_of_cluster; rank++) {
        if (rank == process_rank) {
            for (uint32_t i = init_i; i < end_i; i++) {
                fprintf(output,"data[%d][%d] = %lf\n", i, t.pts_num - 1, data[i][t.pts_num - 1]);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    dataFree(&state);
    fclose(output);
    return 0;
}
