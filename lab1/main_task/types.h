#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>

typedef struct {
    double len;
    uint32_t pts_num;
    double step;
} dim;

typedef struct {
    dim x;
    dim t;
} dimensions;

typedef struct {
    double **data;
    dimensions d;
    double d_coeff;
} data_xy;

#endif