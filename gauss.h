#ifndef GAUSS_H_INCLUDED
#define GAUSS_H_INCLUDED
#include <stddef.h>

#include "matrix.h"
#include "operations_raws.h"
#include "operations.h"

void permute_rows(struct matrix* m, size_t row1, size_t row2);

void multiply_row2(struct matrix* m, double scalar, size_t row);

struct matrix* add_multiply_row_by_scalar(struct matrix* m, size_t k, size_t t, double scalar);

struct matrix* inverse_matrix(struct matrix *m);

struct matrix* solve_system(struct matrix* A, struct matrix* B);

#endif // GAUSS_H_INCLUDED
