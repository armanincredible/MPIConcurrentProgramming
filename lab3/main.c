#include "Matrix.h"
#include "time.h"
#include <omp.h>


Matrix* do_strassen(Matrix* A, Matrix* B, enum TYPE_EXECUTION type);

void initTestData(Matrix* mtx1, Matrix* mtx2) {
    assert(mtx1 != NULL);
    assert(mtx2 != NULL);


    for (int i = 0; i < mtx1->rows; i++) {
        for (int j = 0; j < mtx1->cols; j++) {
            mtx1->values[j + i*mtx1->cols] = j + i*mtx1->cols;
        }
    }
    for (int i = 0; i < mtx2->rows; i++) {
        for (int j = 0; j < mtx2->cols; j++) {
            mtx2->values[j + i*mtx2->cols] = mtx2->cols*mtx2->rows - j - i*mtx2->cols;
        }
    }
}

Matrix* do_classic(Matrix* A, Matrix* B)
{
    Matrix* mtxClassic = matrix_new(A->cols, A->cols);

    printf("Start classic matrix solving\n");

    clock_t start_time = clock();

    mtxClassic = matrix_ijk_matmul(mtxClassic, A, B);

    clock_t end_time = clock();
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Running time: %f seconds\n", time_spent);

    return mtxClassic;
}

void test(unsigned size){
    printf ("-----------Results for size %u-----------\n\n", size);

    Matrix* mtx1 = matrix_new(size, size);
    Matrix* mtx2 = matrix_new(size, size);

    initTestData(mtx1, mtx2);

    // matrix_print("Matrix1:", mtx1);
    // matrix_print("Matrix2:", mtx2);

    Matrix* mtxClassic = do_classic(mtx1, mtx2);

    //matrix_print("Answer Classic:", mtxClassic);

    //printf("\n\n");

    Matrix* mtx3 = do_strassen(mtx1, mtx2, NO_PARALLEL);

    Matrix* mtx4 = do_strassen(mtx1, mtx2, PARALLEL);
  
    //matrix_print("Answer Strassen:", mtx3);

    //printf("\n\n");

    //matrix_print("Answer Strassen:", mtx4);

    matrix_free(&mtx1);
    matrix_free(&mtx2);
    matrix_free(&mtx3);
    matrix_free(&mtxClassic);

    printf("\n");
}

int main() {

    test(256);
    test(512);
    test(1024);
    test(2048);
}