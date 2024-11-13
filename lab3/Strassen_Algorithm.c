#include "Matrix.h"
#include <assert.h>
#include <omp.h>


static Matrix* strassen(Matrix* dest, const Matrix* srcA, const Matrix* srcB, uint32_t rows);


double start_time_parallel;

static void omp_start_time(const char* str)
{
    printf("%s", str);
    start_time_parallel = omp_get_wtime();  // Запуск таймера OpenMP
}

static void omp_end_time()
{
    double end_time_parallel = omp_get_wtime();  // Остановка таймера OpenMP

    double time_spent = end_time_parallel - start_time_parallel;  // Расчёт времени
    printf("Running time:   %f seconds\n", time_spent);
}


Matrix* do_strassen(Matrix* A, Matrix* B, enum TYPE_EXECUTION type) {
    assert(A != NULL);
    assert(B != NULL);

    int isPowerOfTwo = A->rows && !(A->rows & (A->rows - 1));

    if ((A->rows != A->cols) || (B->cols != B->rows) || (B->cols != A->cols) ||  !(isPowerOfTwo))
    {
        fprintf(stderr, "[ERROR]: Failed to calculate: wrong matrix size\nMatrixA.rows=%d, MatrixA.columns=%d\n"
                        "MatrixB.rows=%d. MatrixB.columns=%d\n", A->rows, A->cols, B->rows, B->cols);

        return NULL;
    }

    Matrix* dstMatrix = matrix_new(A->rows, B->cols);

    if (type == NO_PARALLEL)
    {
        omp_start_time("Start strassen no parallel\n");
        dstMatrix = strassen(dstMatrix, A, B, A->rows);
    }
    else
    {
        omp_start_time("Start strassen parallel\n");
        dstMatrix = strassen_omp(dstMatrix, A, B, A->rows);
    }

    omp_end_time();

    return dstMatrix;
}

static Matrix* strassen(Matrix* dest, const Matrix* srcA, const Matrix* srcB, uint32_t rows){

    if (rows == 2) {
        #ifdef DEBUG
            fprintf(stderr, "[DEBUG]: ROWS == 2, calculating matrix_ijk_matmul\n");
        #endif

        return matrix_ijk_matmul(dest, srcA, srcB);
    }

    uint32_t len = rows / 2;

    matrix_autofree Matrix *a11 = matrix_new(len, len), *a12 = matrix_new(len, len), *a21 = matrix_new(len, len), *a22 = matrix_new(len, len),
                           *b11 = matrix_new(len, len), *b12 = matrix_new(len, len), *b21 = matrix_new(len, len), *b22 = matrix_new(len, len), 
                           *c11 = matrix_new(len, len), *c12 = matrix_new(len, len), *c21 = matrix_new(len, len), *c22 = matrix_new(len, len), 
                           *m1 = matrix_new(len, len), *m2 = matrix_new(len, len), *m3 = matrix_new(len, len), *m4 = matrix_new(len, len), 
                           *m5 = matrix_new(len, len), *m6 = matrix_new(len, len), *m7 = matrix_new(len, len), 
                           *temp1 = matrix_new(len, len), *temp2 = matrix_new(len, len);

    /* Divide matrix into four parts */
    FOR(i, len)
    {
        FOR(j, len)
        {
            MXY(a11, i, j) = MXY(srcA, i, j);
            MXY(a12, i, j) = MXY(srcA, i, j + len);
            MXY(a21, i, j) = MXY(srcA, i + len, j);
            MXY(a22, i, j) = MXY(srcA, i + len, j + len);

            MXY(b11, i, j) = MXY(srcB, i, j);
            MXY(b12, i, j) = MXY(srcB, i, j + len);
            MXY(b21, i, j) = MXY(srcB, i + len, j);
            MXY(b22, i, j) = MXY(srcB, i + len, j + len);
        }
    }

    /* Calculate seven formulas of strassen Algorithm */
    strassen(m1, matrix_addition(temp1, a11, a22), matrix_addition(temp2, b11, b22), len);
    strassen(m2, matrix_addition(temp1, a21, a22), b11, len);
    strassen(m3, a11, matrix_subtraction(temp1, b12, b22), len);
    strassen(m4, a22, matrix_subtraction(temp1, b21, b11), len);
    strassen(m5, matrix_addition(temp1, a11, a12), b22, len);
    strassen(m6, matrix_subtraction(temp1, a21, a11), matrix_addition(temp2, b11, b12), len);
    strassen(m7, matrix_subtraction(temp1, a12, a22), matrix_addition(temp2, b21, b22), len);

     /* Merge the answer of matrix dest */
    /* c11 = m1 + m4 - m5 + m7 = m1 + m4 - (m5 - m7) */
    matrix_subtraction(c11, matrix_addition(temp1, m1, m4), matrix_subtraction(temp2, m5, m7));
    matrix_addition(c12, m3, m5);
    matrix_addition(c21, m2, m4);
    matrix_addition(c22, matrix_subtraction(temp1, m1, m2), matrix_addition(temp2, m3, m6));

    /* Store the answer of matrix multiplication */
    FOR(i, len)
    {
        FOR(j, len)
        {
            MXY(dest, i, j) = MXY(c11, i, j);
            MXY(dest, i, j + len) = MXY(c12, i, j);
            MXY(dest, i + len, j) = MXY(c21, i, j);
            MXY(dest, i + len, j + len) = MXY(c22, i, j);
        }
    }

    return dest;
}


Matrix* strassen_omp(Matrix* dest, const Matrix* srcA, const Matrix* srcB, uint32_t rows) {
    omp_set_num_threads(8);
    if (rows <= 64) {  // Ограничиваем переход на обычное умножение для малых матриц
        return matrix_ijk_matmul(dest, srcA, srcB);
    }

    uint32_t len = rows / 2;

    // Создаём вспомогательные матрицы
    matrix_autofree Matrix *m1 = matrix_new(len, len), *m2 = matrix_new(len, len), *m3 = matrix_new(len, len),
                           *m4 = matrix_new(len, len), *m5 = matrix_new(len, len), *m6 = matrix_new(len, len), *m7 = matrix_new(len, len),
                           *c11 = matrix_new(len, len), *c12 = matrix_new(len, len), *c21 = matrix_new(len, len), *c22 = matrix_new(len, len),
                           *temp1 = matrix_new(len, len);

    // Разбиваем исходные матрицы на подматрицы
    Matrix *a11 = submatrix(srcA, 0, 0, len);
    Matrix *a12 = submatrix(srcA, 0, len, len);
    Matrix *a21 = submatrix(srcA, len, 0, len);
    Matrix *a22 = submatrix(srcA, len, len, len);

    Matrix *b11 = submatrix(srcB, 0, 0, len);
    Matrix *b12 = submatrix(srcB, 0, len, len);
    Matrix *b21 = submatrix(srcB, len, 0, len);
    Matrix *b22 = submatrix(srcB, len, len, len);

    // Используем OpenMP для параллельного вычисления только самых трудоёмких операций — умножений матриц
    #pragma omp parallel
    {
        #pragma omp single // Гарантируем, что каждая из задач будет выполнена только одним потоком
        {
            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len), *temp2 = matrix_new(len, len);
                matrix_addition(temp1, a11, a22);
                matrix_addition(temp2, b11, b22);
                strassen_omp(m1, temp1, temp2, len);  // M1 = (A11 + A22) * (B11 + B22)
            }

            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len);
                matrix_addition(temp1, a21, a22);
                strassen_omp(m2, temp1, b11, len);    // M2 = (A21 + A22) * B11
            }

            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len);
                matrix_subtraction(temp1, b12, b22);
                strassen_omp(m3, a11, temp1, len);    // M3 = A11 * (B12 - B22)
            }

            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len);
                matrix_subtraction(temp1, b21, b11);
                strassen_omp(m4, a22, temp1, len);    // M4 = A22 * (B21 - B11)
            }

            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len);
                matrix_addition(temp1, a11, a12);
                strassen_omp(m5, temp1, b22, len);    // M5 = (A11 + A12) * B22
            }

            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len), *temp2 = matrix_new(len, len);
                matrix_subtraction(temp1, a21, a11);
                matrix_addition(temp2, b11, b12);
                strassen_omp(m6, temp1, temp2, len);  // M6 = (A21 - A11) * (B11 + B12)
            }

            #pragma omp task
            {
                matrix_autofree Matrix *temp1 = matrix_new(len, len), *temp2 = matrix_new(len, len);
                matrix_subtraction(temp1, a12, a22);
                matrix_addition(temp2, b21, b22);
                strassen_omp(m7, temp1, temp2, len);  // M7 = (A12 - A22) * (B21 + B22)
            }

            #pragma omp taskwait
        }
        #pragma omp taskwait  // Ждём завершения всех задач
    }

    // Вычисляем C11, C12, C21, C22
    matrix_subtraction(temp1, m5, m7);
    matrix_addition(c11, m1, m4);
    matrix_subtraction(c11, c11, temp1);  // C11 = M1 + M4 - (M5 - M7)

    matrix_addition(c12, m3, m5);         // C12 = M3 + M5
    matrix_addition(c21, m2, m4);         // C21 = M2 + M4

    matrix_subtraction(temp1, m1, m2);
    matrix_addition(c22, temp1, m3);
    matrix_addition(c22, c22, m6);        // C22 = M1 - M2 + M3 + M6

    // Объединяем подматрицы C11, C12, C21, C22 в результирующую матрицу dest
    merge_submatrices(dest, c11, c12, c21, c22);

    return dest;
}

