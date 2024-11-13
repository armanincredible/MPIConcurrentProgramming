#ifndef MATRIX_H_
#define MATRIX_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

enum TYPE_EXECUTION
{
    NO_PARALLEL,
    PARALLEL,
};

typedef struct _Matrix
{
    uint32_t rows;
    uint32_t cols;
    uint32_t *values;
} Matrix;

#define MXY_NOT_POINTER_(M_, R_, C_) ((int(*)[M_.rows])(M_.values))[R_][C_]

#define GET_VAL_FROM_MATRIX(mat, i, j) (mat->values[(i) * mat->cols + (j)])


// mx as a pointer
#define MXY(M_, R_, C_) MXY_NOT_POINTER_(((Matrix)*M_), R_, C_)

#define FOR(C_, N_) for (uint32_t C_ = 0; C_ < N_; ++C_)

/*
    Matrix mx = { 3,2, calloc(3 * 2, sizeof(int))  } ;
    MXY(mx,1,1) = 42;
    assert( MXY(mx,1,1) == 42) ;
*/

static inline void matrix_free(Matrix **ptr)
{
    assert(*ptr && (*ptr)->values);
    free((*ptr)->values);
    free(*ptr);
    *ptr = 0;
}

#define matrix_autofree __attribute__((cleanup(matrix_free)))


Matrix* matrix_new(uint32_t rows, uint32_t cols);
void matrix_print(const char *prompt, Matrix *mx);

Matrix* matrix_ijk_matmul(Matrix* dest, const Matrix* srcA, const Matrix* srcB);

void matrix_foreach(Matrix *mx, Matrix *(*callback)(Matrix *, unsigned, unsigned));

int matrix_equal(const Matrix *a, const Matrix *b);
Matrix* matrix_addition(Matrix *res, const Matrix *a, const Matrix *b);
Matrix* matrix_subtraction(Matrix *res, const Matrix *a, const Matrix *b);
Matrix* submatrix(const Matrix* src, uint32_t row_start, uint32_t col_start, uint32_t size);

void merge_submatrices(Matrix* dest, const Matrix* c11, const Matrix* c12, const Matrix* c21, const Matrix* c22);


Matrix* strassen_iterative(Matrix* dest,  Matrix* srcA,  Matrix* srcB, uint32_t rows);
Matrix* strassen_omp(Matrix* dest, const Matrix* srcA, const Matrix* srcB, uint32_t rows);

#endif /* MATRIX_H_ */