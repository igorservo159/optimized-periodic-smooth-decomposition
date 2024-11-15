#ifndef FOURIER_H
#define FOURIER_H

#include "common.h"
#include "utils.h"

void compute_cfft2d(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_cifft2d(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_cfft2d_column_row(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns);
void compute_cperiodic_border_B(MKL_Complex8 *I_t, MKL_Complex8 *B_t, size_t rows, size_t columns);
void compute_cfft2d_of_border_B(MKL_Complex8 *B_t_B_w, size_t rows, size_t columns);
void compute_csmooth_component_S(MKL_Complex8 *B_S, size_t rows, size_t columns);
void compute_csmooth_component_S_2(MKL_Complex8 *B_S, size_t rows, size_t columns);
void compute_cperiodic_component_P(MKL_Complex8 *I_w, MKL_Complex8 *S, size_t rows, size_t columns);
void compute_cfftshift(MKL_Complex8 *vector, size_t rows, size_t columns);
void normalize_fvector(float *fvector, size_t size);

void compute_zfft2d(MKL_Complex16 *I_t_I_w, size_t rows, size_t columns);
void compute_zifft2d(MKL_Complex16 *I_t_I_w, size_t rows, size_t columns);
void compute_zfft2d_column_row(MKL_Complex16 *I_t_I_w, size_t rows, size_t columns);
void compute_zperiodic_border_B(MKL_Complex16 *I_t, MKL_Complex16 *B_t, size_t rows, size_t columns);
void compute_zfft2d_of_border_B(MKL_Complex16 *B_t_B_w, size_t rows, size_t columns);
void compute_zsmooth_component_S(MKL_Complex16 *B_S, size_t rows, size_t columns);
void compute_zsmooth_component_S_2(MKL_Complex16 *B_S, size_t rows, size_t columns);
void compute_zperiodic_component_P(MKL_Complex16 *I_w, MKL_Complex16 *S, size_t rows, size_t columns);
void compute_zfftshift(MKL_Complex16 *vector, size_t rows, size_t columns);
void normalize_dvector(double *dvector, size_t size);

#endif