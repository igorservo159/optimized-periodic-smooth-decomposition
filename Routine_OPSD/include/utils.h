#ifndef UTIL_H
#define UTIL_H

#include "common.h"

int check_args(const char *BIN, const char *ROUTINE, const char *PRECISION, const char *SAVE_VECTORS, const char *INPUT);

void init_cvector(MKL_Complex8 **vector, size_t size);
void free_cvector(MKL_Complex8 *vector);
void show_cmatrix(MKL_Complex8 *matrix, size_t rows, size_t columns);
void fill_cmatrix(MKL_Complex8 *matrix, size_t rows, size_t columns, unsigned int seed);
void read_cvector_bin(const char *filename, MKL_Complex8 *vector, size_t size);
void save_cvector_on_bin(const char *filename, MKL_Complex8 *vector, size_t size);
void copy_cvector_to_real_fvector(MKL_Complex8 *cvector, float *fvector, size_t size);
void copy_cvector_to_real_dvector(MKL_Complex8 *cvector, double *dvector, size_t size);

void init_zvector(MKL_Complex16 **vector, size_t size);
void free_zvector(MKL_Complex16 *vector);
void show_zmatrix(MKL_Complex16 *matrix, size_t rows, size_t columns);
void fill_zmatrix(MKL_Complex16 *matrix, size_t rows, size_t columns, unsigned int seed);
void read_zvector_bin(const char *filename, MKL_Complex16 *vector, size_t size);
void save_zvector_on_bin(const char *filename, MKL_Complex16 *vector, size_t size);
void copy_zvector_to_real_fvector(MKL_Complex16 *zvector, float *fvector, size_t size);
void copy_zvector_to_real_dvector(MKL_Complex16 *zvector, double *dvector, size_t size);

void init_fvector(float **vector, size_t size);
void free_fvector(float *vector);
void read_fvector_bin(const char *filename, float *vector, size_t size);
void save_fvector_on_bin(const char *filename, float *vector, size_t size);
void copy_fvector_to_cvector(MKL_Complex8 *cvector, float *fvector, size_t size);
void copy_fvector_to_zvector(MKL_Complex16 *zvector, float *fvector, size_t size);

void init_dvector(double **vector, size_t size);
void free_dvector(double *vector);
void read_dvector_bin(const char *filename, double *vector, size_t size);
void save_dvector_on_bin(const char *filename, double *vector, size_t size);
void copy_dvector_to_cvector(MKL_Complex8 *cvector, double *dvector, size_t size);
void copy_dvector_to_zvector(MKL_Complex16 *zvector, double *dvector, size_t size);

#endif