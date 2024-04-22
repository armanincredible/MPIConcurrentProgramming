#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "types.h"
#include "solver.h"
#include <chrono>
#include <iostream>


static double **dataAlloc(uint32_t pts_x, uint32_t pts_y) {
    double **data = (double **) calloc(pts_x, sizeof(data[0]));
    
    for (uint32_t i = 0; i < pts_x; i++) {
        data[i] = (double *) calloc(pts_y, sizeof(data[i][0]));
    }

    return data;
}

static void calculateRow(data_xy *state, uint32_t row) {
    if ((row < 1) || (row >= state->d.t.pts_num)) return;
    for (uint32_t i = 1; i < state->d.x.pts_num - 1; i++) {
        state->data[i][row] = calculateNewPoint(state, i, row - 1);
    }
}

static void dataFree(data_xy *state) {
    for (uint32_t i = 0; i < state->d.x.pts_num; i++) {
        free(state->data[i]);
    }

    free(state->data);
}


int main(int argc, const char **argv) {
    FILE* output = fopen("output_seq.txt", "w");
    if (output == NULL) {
        fprintf(stderr, "error in opening output file\n");
    }

    dim x = {1.0, 1000, 0.0};
    x.step = x.len / x.pts_num;

    dim t = {20.0, 1000000, 0.0};
    t.step = t.len / t.pts_num;

    const dimensions dims = {x, t};

    double **data = dataAlloc(x.pts_num, t.pts_num);

    initStartCoordDim(data, x.pts_num, x.step);

    initStartTimeDim(data, t.pts_num, t.step);

    // initial conditions at the end of the grid for the time space
    initConditions(data, t.pts_num, x.pts_num, t.step);

    data_xy state = {data, dims, 1e-3};

    auto start_time = std::chrono::steady_clock::now();

    for (uint32_t row = 1; row < t.pts_num; row++) {
        calculateRow(&state, row);
    }

    auto stop_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (stop_time - start_time).count();
    std::cout << "Elapsed: " << elapsed / 1000. << " sec \n";

    for (uint32_t i = 1; i < x.pts_num; i++) {
        fprintf(output, "data[%d][%d] = %lf\n", i, t.pts_num - 1, data[i][t.pts_num - 1]);
    }
    draw_results(&state);
    dataFree(&state);
    fclose(output);
    return 0;
}
