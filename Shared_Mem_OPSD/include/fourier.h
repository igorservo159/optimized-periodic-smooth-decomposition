#ifndef FOURIER_H
#define FOURIER_H

#include "common.h"

void compute_fft2d(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_ifft2d(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_fft2d_column_row(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_fft2d_column_row_2(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_periodic_border_B(MKL_Complex8 *I_t, MKL_Complex8 *B_t, size_t rows, size_t columns);
void compute_fft2d_of_border_B(MKL_Complex8 *B_t_B_w, size_t rows, size_t columns);
void compute_smooth_component_S(MKL_Complex8 *B_S, size_t rows, size_t columns);
void compute_smooth_component_S_2(MKL_Complex8 *B_S, size_t rows, size_t columns);
void compute_periodic_component_P(MKL_Complex8 *I_w, MKL_Complex8 *S, size_t rows, size_t columns);
void compute_fftshift(MKL_Complex8 *vector, size_t rows, size_t columns);

#endif