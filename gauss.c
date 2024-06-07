#include <stdio.h>
#include <stdlib.h>
#include "gauss.h"


void permute_rows(struct matrix* m, size_t row1, size_t row2) {
    if (row1 == row2) return;

    for (size_t i = 0; i < m->cols; ++i) {
        double temp = *get_element(m, row1, i);
        *get_element(m, row1, i) = *get_element(m, row2, i);
        *get_element(m, row2, i) = temp;
    }
}

void multiply_row2(struct matrix* m, double scalar, size_t row) {
    for (size_t i = 0; i < m->cols; ++i) {
        *get_element(m, row, i) *= scalar;
    }
}

struct matrix* add_multiply_row_by_scalar(struct matrix* m, size_t target_row, size_t source_row, double scalar) {
    for (size_t i = 0; i < m->cols; ++i) {
        *get_element(m, target_row, i) += scalar * *get_element(m, source_row, i);
    }
    return m;
}

struct matrix* inverse_matrix(struct matrix* m) {
    if (m->rows != m->cols) {
        return NULL;
    }

    struct matrix* result = matrix_alloc_square_unit(m->rows);
    if (result == NULL) {
        return NULL;
    }

    struct matrix* temp = matrix_clone(m);
    if (temp == NULL) {
        matrix_free(result);
        return NULL;
    }

    for (size_t i = 0; i < temp->rows; ++i) {
        size_t max_row = i;
        for (size_t row = i + 1; row < temp->rows; ++row) {
            if (fabs(*get_element(temp, row, i)) > fabs(*get_element(temp, max_row, i))) {
                max_row = row;
            }
        }

        if (*get_element(temp, max_row, i) == 0) {
            matrix_free(result);
            matrix_free(temp);
            return NULL;
        }

        permute_rows(temp, i, max_row);
        permute_rows(result, i, max_row);

        double element = *get_element(temp, i, i);
        multiply_row2(temp, 1.0 / element, i);
        multiply_row2(result, 1.0 / element, i);

        for (size_t j = 0; j < temp->rows; ++j) {
            if (j != i) {
                double scalar = -*get_element(temp, j, i);
                add_multiply_row_by_scalar(temp, j, i, scalar);
                add_multiply_row_by_scalar(result, j, i, scalar);
            }
        }
    }

    matrix_free(temp);
    return result;
}

struct matrix* solve_system(struct matrix* A, struct matrix* B) {
    if (A->rows != A->cols || A->rows != B->rows) {
        return NULL;
    }

    struct matrix *A_inv = inverse_matrix(A);
    if (A_inv == NULL) {
        return NULL;
    }

    struct matrix *x = matrix_multiplication(A_inv, B);
    matrix_free(A_inv);
    return x;
}

