#ifndef AUX_H
#define AUX_H

#include <mkl.h>
#include "mkl_dfti.h"
#include <math.h>
#include "stdio.h"
#include "stdlib.h"
#include "complex.h"
#include "time.h"
#include "mkl.h"
#include "omp.h"

#define RAND_MAX_F 2147483647.0f
#define PI 3.14159265358979323846

void init_complex_matrix(MKL_Complex8 **matrix, size_t rows, size_t columns);
void init_float_vector(float **v, size_t size);
void free_matrix(MKL_Complex8 *matrix);
void show_matrix(MKL_Complex8 *matrix, size_t rows, size_t columns);
void show_ram_allocation(size_t rows, size_t columns);
void fill(MKL_Complex8 *matrix, size_t rows, size_t columns, unsigned int seed);
void compute_fft2D_column_row(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_fft2D_column_row_2(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_periodic_border_B(MKL_Complex8 *I_t, MKL_Complex8 *B_t, size_t rows, size_t columns);
void compute_fft2D_of_B(MKL_Complex8 *B_t_B_w, size_t rows, size_t columns);
void compute_smooth_component_S(MKL_Complex8 *B_S, size_t rows, size_t columns);
void compute_smooth_component_S_2(MKL_Complex8 *B_S, size_t rows, size_t columns);
void compute_periodic_component_P(MKL_Complex8 *I_w, MKL_Complex8 *S, size_t rows, size_t columns);
void read_binary(float *v, size_t size);
void read_matrix(MKL_Complex8 *matrix, size_t rows, size_t columns);

#endif