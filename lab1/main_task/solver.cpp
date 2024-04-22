#include <stdint.h>
#include "types.h"
#include <math.h>
#include <raylib.h>
#include <stdio.h>

void initStartCoordDim(double** data, uint32_t pts_num, double step) {
    for (uint32_t i = 0; i < pts_num; i++) {
        data[i][0] = sin(step * i * 5);
    }
}

void initStartTimeDim(double** data, uint32_t pts_num, double step) {
    for (uint32_t i = 0; i < pts_num; i++) {
        data[0][i] = data[0][0] + (step * i);
    }
}

void initConditions(double** data, uint32_t t_pts_num, uint32_t x_pts_num, double step){
    for (uint32_t i = 1; i < t_pts_num; i++) {
        data[x_pts_num - 1][i] = data[x_pts_num - 1][0] - (step * i);
    }
}

double calculateNewPoint(data_xy *state, uint32_t x_pos, uint32_t y_pos) {
    double **data = state->data;

    double r_side = (data[x_pos + 1][y_pos] - 2 * data[x_pos][y_pos] + data[x_pos - 1][y_pos]) / state->d.x.step / state->d.x.step;

    return data[x_pos][y_pos] + state->d.t.step * state->d_coeff * r_side;
}

static double max_abs(data_xy *state) {
    double res = 0.0;

    for (uint32_t i = 0; i < state->d.x.pts_num; i++) {
        for (uint32_t j = 0; j < state->d.t.pts_num; j++) {
            double abs_data = abs(state->data[i][j]);
            if (abs_data > res) {
                res = abs_data;
            }
        }
    }

    return res ? res : 1.0;
}

void draw_results(data_xy *state) {
    const double scale_x = 1000.0 / state->d.t.pts_num;
    const double scale_y = 600.0 / state->d.x.pts_num;
    const double padding = 10;

    InitWindow((state->d.t.pts_num + padding) * scale_x, (state->d.x.pts_num + padding) * scale_y, "Temperature");
    SetTargetFPS(1);

    double max_val = max_abs(state);

    uint32_t step_x = (uint32_t)(1.0/scale_x);
    uint32_t step_y = (uint32_t)(1.0/scale_y);

    printf("step_x = %d, step_y = %d\n", step_x, step_y);

    step_x = step_x == 0 ? 1 : step_x;
    step_y = step_y == 0 ? 1 : step_y;

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(WHITE);

        for (uint32_t j = padding / 2; j < state->d.x.pts_num / step_y; j++) {
            for (uint32_t i = padding / 2; i < state->d.t.pts_num / step_x; i++) {
                double normalized = state->data[j * step_y][i * step_x] / max_val / 2;
                double hue = normalized * 180 + 180;

                DrawPixel(i, j, ColorFromHSV(hue, 1.0, 1.0));
            }
        }

        EndDrawing();
    }

    CloseWindow();
}
