#ifndef UTIL_H
#define UTIL_H

#include "common.h"

void init_cvector(MKL_Complex8 **vector, size_t size);
void init_fvector(float **vector, size_t size);
void free_cvector(MKL_Complex8 *vector);
void free_fvector(float *vector);
void show_cmatrix(MKL_Complex8 *matrix, size_t rows, size_t columns);
void fill_cmatrix(MKL_Complex8 *matrix, size_t rows, size_t columns, unsigned int seed);
void read_fvector_bin(const char *filename, float *vector, size_t size);
void read_cvector_bin(const char *filename, MKL_Complex8 *vector, size_t size);
void save_cvector_on_bin(const char *filename, MKL_Complex8 *vector, size_t size);
void save_fvector_on_bin(const char *filename, float *vector, size_t size);
void copy_cvector_to_real_fvector(MKL_Complex8 *cvector, float *fvector, size_t size);
void copy_fvector_to_cvector(MKL_Complex8 *cvector, float *fvector, size_t size);

#endif