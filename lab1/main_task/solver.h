#ifndef SOVLER_H
#define SOVLER_H

#include "types.h"
#include <math.h>


void initStartCoordDim(double** data, uint32_t pts_num, double step);
void initStartTimeDim(double** data, uint32_t pts_num, double step);
void initConditions(double** data, uint32_t t_pts_num, uint32_t x_pts_num, double step);
double calculateNewPoint(data_xy *state, uint32_t x_pos, uint32_t y_pos);
void draw_results(data_xy *state);

#endif