#include "Matrix.h"

Matrix* matrix_addition(Matrix *res, const Matrix *a, const Matrix *b)
{
    assert(res && a && b);
    FOR(i, res->rows)
    FOR(j, res->cols)
    // res->values[i][j] = a->values[i][j] + b->values[i][j];
    MXY(res, i, j) = MXY(a, i, j) + MXY(b, i, j);

    return res;
}

Matrix* matrix_subtraction(Matrix *res, const Matrix *a, const Matrix *b)
{
    assert(res && a && b);
    FOR(i, res->rows)
    FOR(j, res->cols)
    MXY(res, i, j) = MXY(a, i, j) - MXY(b, i, j);

    return res;
}

int matrix_equal(const Matrix *a, const Matrix *b)
{
    assert(a && b);

    if (a->rows != b->rows)
        return 0;
    if (a->cols != b->cols)
        return 0;

    FOR(i, a->rows)
    FOR(j, a->cols)
    if (MXY(a, i, j) != MXY(b, i, j))
        return 0;

    return 1;
}

Matrix* matrix_ijk_matmul(Matrix* dest, const Matrix* srcA, const Matrix* srcB) {
    if (srcA->cols != srcB->rows) {
        #ifdef DEBUG
            fprintf(stderr, "[ERROR]: srcA->cols != srcB->rows");
        #endif

        return NULL; 
    }

    if (dest->rows != srcA->rows || dest->cols != srcB->cols) {
        #ifdef DEBUG
            fprintf(stderr, "[ERROR]: dest->rows != srcA->rows || dest->cols != srcB->cols\n");
            fprintf(stderr, "[ERROR]: dest->rows = %d; dest->cols = %d;  srcA->rows = %d; srcB->cols = %d; \n", dest->rows, dest->cols, srcA->rows, srcB->cols);
        #endif

        return NULL; 
    }

    for (uint32_t i = 0; i < dest->rows; i++) {
        for (uint32_t j = 0; j < dest->cols; j++) {
            GET_VAL_FROM_MATRIX(dest, i, j) = 0;
        }
    }

    for (uint32_t i = 0; i < srcA->rows; i++) {
        for (uint32_t j = 0; j < srcB->cols; j++) {
            for (uint32_t k = 0; k < srcA->cols; k++) {
                GET_VAL_FROM_MATRIX(dest, i, j) += GET_VAL_FROM_MATRIX(srcA, i, k) * GET_VAL_FROM_MATRIX(srcB, k, j);
            }
        }
    }

    return dest;
}

inline Matrix* matrix_new(uint32_t rows, uint32_t cols)
{
    Matrix *matrix = calloc(1, sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->values = calloc(rows * cols, sizeof(int));
    return matrix;
}


inline Matrix* submatrix(const Matrix* src, uint32_t row_start, uint32_t col_start, uint32_t size) {
    Matrix* result = matrix_new(size, size);

    FOR(i, size) {
        FOR(j, size) {
            MXY(result, i, j) = MXY(src, i + row_start, j + col_start);
        }
    }

    return result;
}


void merge_submatrices(Matrix* dest, const Matrix* c11, const Matrix* c12, const Matrix* c21, const Matrix* c22) {
    uint32_t len = c11->rows;  // Предполагается, что размеры подматриц одинаковы

    // Вставляем c11 в верхний левый угол
    FOR(i, len) {
        FOR(j, len) {
            MXY(dest, i, j) = MXY(c11, i, j);
        }
    }

    // Вставляем c12 в верхний правый угол
    FOR(i, len) {
        FOR(j, len) {
            MXY(dest, i, j + len) = MXY(c12, i, j);
        }
    }

    // Вставляем c21 в нижний левый угол
    FOR(i, len) {
        FOR(j, len) {
            MXY(dest, i + len, j) = MXY(c21, i, j);
        }
    }

    // Вставляем c22 в нижний правый угол
    FOR(i, len) {
        FOR(j, len) {
            MXY(dest, i + len, j + len) = MXY(c22, i, j);
        }
    }
}

void matrix_print(const char *prompt, Matrix *mx)
{
    assert(prompt && mx);
    printf("%s\n", prompt);

    #ifdef DEBUG
        printf("rows num: %d; column num: %d; value: %p;\n", mx->rows, mx->cols, mx->values);
    #endif

    FOR(i, mx->rows)
    {
        FOR(j, mx->cols)
        {
            printf(" %5d ", GET_VAL_FROM_MATRIX(mx, i, j));
        }
        printf("\n");
    }
}

void matrix_foreach(Matrix *mx, Matrix *(*callback)(Matrix *, unsigned, unsigned))
{
    FOR(i, mx->rows)
    {
        FOR(j, mx->cols)
        {
            callback(mx, i, j);
        }
    }
}