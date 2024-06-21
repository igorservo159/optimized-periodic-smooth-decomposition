#include "../include/fourier.h"

//CVECTOR FFT FUNCTIONS

void compute_cfft2d(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns){
    if (I_t_I_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[2] = {rows, columns};

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                                  DFTI_COMPLEX, 2, dim_sizes);
    
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeForward(desc_handle_dim1, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);
}

void compute_cifft2d(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns){
    if (I_t_I_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    double start = omp_get_wtime();

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[2] = {rows, columns};

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                                  DFTI_COMPLEX, 2, dim_sizes);
    
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeBackward(desc_handle_dim1, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);     
}

void compute_cfft2d_column_row(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns)
{
    if (I_t_I_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim2 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, rows);
    status = DftiCreateDescriptor(&desc_handle_dim2, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, columns);

    MKL_LONG stride[2] = {0, columns};

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, columns);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_DISTANCE, 1);
    status = DftiCommitDescriptor(desc_handle_dim1);

    stride[1] = 1;

    status = DftiSetValue(desc_handle_dim2, DFTI_NUMBER_OF_TRANSFORMS, rows);
    status = DftiSetValue(desc_handle_dim2, DFTI_INPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim2, DFTI_OUTPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim2, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim2, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim2);

    status = DftiComputeForward(desc_handle_dim1, I_t_I_w);
    status = DftiComputeForward(desc_handle_dim2, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);
    status = DftiFreeDescriptor(&desc_handle_dim2);
}

void compute_cperiodic_border_B(MKL_Complex8 *I_t, MKL_Complex8 *B_t, size_t rows, size_t columns)
{
    if (I_t == NULL || B_t == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    B_t[0].real = I_t[columns - 1].real - 2 * I_t[0].real + I_t[(rows - 1) * columns].real;
    B_t[0].imag = I_t[columns - 1].imag - 2 * I_t[0].imag + I_t[(rows - 1) * columns].imag;

    B_t[columns - 1].real = I_t[0].real - 2 * I_t[columns - 1].real + I_t[rows * columns - 1].real;
    B_t[columns - 1].imag = I_t[0].imag - 2 * I_t[columns - 1].imag + I_t[rows * columns - 1].imag;

    B_t[(rows - 1) * columns].real = I_t[0].real - 2 * I_t[(rows - 1) * columns].real + I_t[rows * columns - 1].real;
    B_t[(rows - 1) * columns].imag = I_t[0].imag - 2 * I_t[(rows - 1) * columns].imag + I_t[rows * columns - 1].imag;

    B_t[rows * columns - 1].real = I_t[columns - 1].real - 2 * I_t[rows * columns - 1].real + I_t[(rows - 1) * columns].real;
    B_t[rows * columns - 1].imag = I_t[columns - 1].imag - 2 * I_t[rows * columns - 1].imag + I_t[(rows - 1) * columns].imag;

    // linha 0 e rows-1
    for (int j = 1; j < columns - 1; j++)
    {
        B_t[j].real = I_t[j + (rows - 1) * columns].real - I_t[j].real;
        B_t[j].imag = I_t[j + (rows - 1) * columns].imag - I_t[j].imag;

        B_t[j + (rows - 1) * columns].real = B_t[j].real * (-1);
        B_t[j + (rows - 1) * columns].imag = B_t[j].imag * (-1);
    }

    // coluna 0 e columns-1
    for (int i = 1; i < rows - 1; i++)
    {
        B_t[i * columns].real = I_t[columns - 1 + i * columns].real - I_t[i * columns].real;
        B_t[i * columns].imag = I_t[columns - 1 + i * columns].imag - I_t[i * columns].imag;

        B_t[columns - 1 + i * columns].real = B_t[i * columns].real * (-1);
        B_t[columns - 1 + i * columns].imag = B_t[i * columns].imag * (-1);
    }
}

void compute_cfft2d_of_border_B(MKL_Complex8 *B_t_B_w, size_t rows, size_t columns)
{
    if (B_t_B_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    // Column-one FFT
    MKL_Complex8 *a = (MKL_Complex8 *)malloc(sizeof(MKL_Complex8));
    a->real = B_t_B_w[0].real + B_t_B_w[columns - 1].real;
    a->imag = B_t_B_w[0].imag + B_t_B_w[columns - 1].imag;

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, rows);

    MKL_LONG stride[2] = {0, columns};

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, 1);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeForward(desc_handle_dim1, B_t_B_w);
    status = DftiFreeDescriptor(&desc_handle_dim1);

    // Calcular cada elemento de v usando a fórmula W^k = exp(-i * 2 * PI * k / M)
    MKL_Complex8 *v = NULL;
    init_cvector(&v, rows);
    if (v == NULL)
    {
        printf("Erro ao alocar memória\n");
        return;
    }

    v[0].real = 0.00f;
    v[0].imag = 0.00f;

    for (int k = 1; k < rows; k++)
    {
        float theta = -2.0f * PI * (rows - k) / rows;
        v[k].real = 1.0f - cosf(theta);
        v[k].imag = sinf(theta) * (-1);
    }

    cblas_ccopy(rows, B_t_B_w, columns, &B_t_B_w[columns - 1], columns);
    cblas_csscal(rows, -1.0f, &B_t_B_w[columns - 1], columns);
    cblas_caxpy(rows, a, v, 1, &B_t_B_w[columns - 1], columns);

    for (int j = 1; j < columns - 1; j++)
    {
        a->real = B_t_B_w[j].real;
        a->imag = B_t_B_w[j].imag;
        cblas_ccopy(rows, v, 1, &B_t_B_w[j], columns);
        cblas_cscal(rows, a, &B_t_B_w[j], columns);
    }

    free(v);
    free(a);

    // Row-by-Row FFT
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim3 = NULL;
    status = DftiCreateDescriptor(&desc_handle_dim3, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, columns);

    stride[1] = 1;

    status = DftiSetValue(desc_handle_dim3, DFTI_NUMBER_OF_TRANSFORMS, rows);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim3);
    status = DftiComputeForward(desc_handle_dim3, B_t_B_w);
    status = DftiFreeDescriptor(&desc_handle_dim3);
}

void compute_csmooth_component_S(MKL_Complex8 *B_S, size_t rows, size_t columns)
{
    if (B_S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    MKL_Complex8 aux = {B_S[0].real, B_S[0].imag};

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            float denom = (2.0f * cosf(2.0f * PI * i / rows) + 2.0f * cosf(2.0f * PI * j / columns) - 4.0f);
            B_S[j + i * columns].real /= denom;
            B_S[j + i * columns].imag /= denom;
        }
    }

    B_S[0].real = aux.real;
    B_S[0].imag = aux.imag;
}

void compute_csmooth_component_S_2(MKL_Complex8 *B_S, size_t rows, size_t columns)
{
    if (B_S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    MKL_Complex8 aux = {B_S[0].real, B_S[0].imag};

    float *v1 = NULL;
    float *v2 = NULL;

    init_fvector(&v1, rows);
    init_fvector(&v2, columns);

    for (int i = 0; i < rows; i++)
    {
        v1[i] = (2.0f * PI * i) / rows;
    }
    
    for (int i = 0; i < columns; i++)
    {
        v2[i] = (2.0f * PI * (i)) / columns;
    }

    float *cos_1 = NULL;
    init_fvector(&cos_1, rows);
    float *cos_2 = NULL;
    init_fvector(&cos_2, columns);
    vsCos(rows, v1, cos_1);
    vsCos(columns, v2, cos_2);
    cblas_sscal(rows, 2.0f, cos_1, 1);
    cblas_sscal(columns, 2.0f, cos_2, 1);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            float denom = cos_1[i] + cos_2[j] - 4.0f;
            B_S[i + j * rows].real /= denom;
            B_S[i + j * rows].imag /= denom;
        }
    }

    B_S[0].real = aux.real;
    B_S[0].imag = aux.imag;
}


void compute_cperiodic_component_P(MKL_Complex8 *I_w, MKL_Complex8 *S, size_t rows, size_t columns)
{
    if (I_w == NULL || S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }
    
    vcSub(rows * columns, I_w, S, I_w);
}

void compute_cfftshift(MKL_Complex8 *vector, size_t rows, size_t columns) {
    size_t half_rows = rows / 2;
    size_t half_columns = columns / 2;
    
    // Troca os quadrantes (1,4) e (2,3)
    for (size_t i = 0; i < half_rows; i++) {
        for (size_t j = 0; j < half_columns; j++) {
            // Índices dos elementos a serem trocados
            size_t index1 = i * columns + j;
            size_t index2 = (i + half_rows) * columns + (j + half_columns);
            size_t index3 = i * columns + (j + half_columns);
            size_t index4 = (i + half_rows) * columns + j;

            // Swap (1,4)
            MKL_Complex8 temp = vector[index1];
            vector[index1] = vector[index2];
            vector[index2] = temp;
            
            // Swap (2,3)
            temp = vector[index3];
            vector[index3] = vector[index4];
            vector[index4] = temp;
        }
    }
}

void normalize_cvector(MKL_Complex8 *cvector, size_t size){
    float *magnitudes = NULL;
    init_fvector(&magnitudes, size);
    
    vcAbs(size, cvector, magnitudes);

    for (int i = 0; i < size; i++) {
        cvector[i].real = cvector[i].real / magnitudes[i];
        cvector[i].imag = cvector[i].imag / magnitudes[i];
    }
}

//ZVECTOR FFT FUNCTIONS

void compute_zfft2d(MKL_Complex16 *I_t_I_w, size_t rows, size_t columns){
    if (I_t_I_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[2] = {rows, columns};

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_DOUBLE,
                                  DFTI_COMPLEX, 2, dim_sizes);
    
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeForward(desc_handle_dim1, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);
}

void compute_zifft2d(MKL_Complex16 *I_t_I_w, size_t rows, size_t columns){
    if (I_t_I_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    double start = omp_get_wtime();

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[2] = {rows, columns};

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_DOUBLE,
                                  DFTI_COMPLEX, 2, dim_sizes);
    
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeBackward(desc_handle_dim1, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);     
}

void compute_zfft2d_column_row(MKL_Complex16 *I_t_I_w, size_t rows, size_t columns)
{
    if (I_t_I_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim2 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_DOUBLE,
                                  DFTI_COMPLEX, 1, rows);
    status = DftiCreateDescriptor(&desc_handle_dim2, DFTI_DOUBLE,
                                  DFTI_COMPLEX, 1, columns);

    MKL_LONG stride[2] = {0, columns};

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, columns);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_DISTANCE, 1);
    status = DftiCommitDescriptor(desc_handle_dim1);

    stride[1] = 1;

    status = DftiSetValue(desc_handle_dim2, DFTI_NUMBER_OF_TRANSFORMS, rows);
    status = DftiSetValue(desc_handle_dim2, DFTI_INPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim2, DFTI_OUTPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim2, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim2, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim2);

    status = DftiComputeForward(desc_handle_dim1, I_t_I_w);
    status = DftiComputeForward(desc_handle_dim2, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);
    status = DftiFreeDescriptor(&desc_handle_dim2);
}

void compute_zperiodic_border_B(MKL_Complex16 *I_t, MKL_Complex16 *B_t, size_t rows, size_t columns)
{
    if (I_t == NULL || B_t == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    B_t[0].real = I_t[columns - 1].real - 2 * I_t[0].real + I_t[(rows - 1) * columns].real;
    B_t[0].imag = I_t[columns - 1].imag - 2 * I_t[0].imag + I_t[(rows - 1) * columns].imag;

    B_t[columns - 1].real = I_t[0].real - 2 * I_t[columns - 1].real + I_t[rows * columns - 1].real;
    B_t[columns - 1].imag = I_t[0].imag - 2 * I_t[columns - 1].imag + I_t[rows * columns - 1].imag;

    B_t[(rows - 1) * columns].real = I_t[0].real - 2 * I_t[(rows - 1) * columns].real + I_t[rows * columns - 1].real;
    B_t[(rows - 1) * columns].imag = I_t[0].imag - 2 * I_t[(rows - 1) * columns].imag + I_t[rows * columns - 1].imag;

    B_t[rows * columns - 1].real = I_t[columns - 1].real - 2 * I_t[rows * columns - 1].real + I_t[(rows - 1) * columns].real;
    B_t[rows * columns - 1].imag = I_t[columns - 1].imag - 2 * I_t[rows * columns - 1].imag + I_t[(rows - 1) * columns].imag;

    // linha 0 e rows-1
    for (int j = 1; j < columns - 1; j++)
    {
        B_t[j].real = I_t[j + (rows - 1) * columns].real - I_t[j].real;
        B_t[j].imag = I_t[j + (rows - 1) * columns].imag - I_t[j].imag;

        B_t[j + (rows - 1) * columns].real = B_t[j].real * (-1);
        B_t[j + (rows - 1) * columns].imag = B_t[j].imag * (-1);
    }

    // coluna 0 e columns-1
    for (int i = 1; i < rows - 1; i++)
    {
        B_t[i * columns].real = I_t[columns - 1 + i * columns].real - I_t[i * columns].real;
        B_t[i * columns].imag = I_t[columns - 1 + i * columns].imag - I_t[i * columns].imag;

        B_t[columns - 1 + i * columns].real = B_t[i * columns].real * (-1);
        B_t[columns - 1 + i * columns].imag = B_t[i * columns].imag * (-1);
    }
}

void compute_zfft2d_of_border_B(MKL_Complex16 *B_t_B_w, size_t rows, size_t columns)
{
    if (B_t_B_w == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    // Column-one FFT
    MKL_Complex16 *a = (MKL_Complex16 *)malloc(sizeof(MKL_Complex16));
    a->real = B_t_B_w[0].real + B_t_B_w[columns - 1].real;
    a->imag = B_t_B_w[0].imag + B_t_B_w[columns - 1].imag;

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_DOUBLE,
                                  DFTI_COMPLEX, 1, rows);

    MKL_LONG stride[2] = {0, columns};

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, 1);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeForward(desc_handle_dim1, B_t_B_w);
    status = DftiFreeDescriptor(&desc_handle_dim1);

    // Calcular cada elemento de v usando a fórmula W^k = exp(-i * 2 * PI * k / M)
    MKL_Complex16 *v = NULL;
    init_zvector(&v, rows); 
    if (v == NULL)
    {
        printf("Erro ao alocar memória\n");
        return;
    }

    v[0].real = 0.00;
    v[0].imag = 0.00;

    for (int k = 1; k < rows; k++)
    {
        float theta = -2.0 * PI * (rows - k) / rows;
        v[k].real = 1.0 - cos(theta);
        v[k].imag = sin(theta) * (-1);
    }

    cblas_ccopy(rows, B_t_B_w, columns, &B_t_B_w[columns - 1], columns);
    cblas_csscal(rows, -1.0, &B_t_B_w[columns - 1], columns);
    cblas_caxpy(rows, a, v, 1, &B_t_B_w[columns - 1], columns);

    for (int j = 1; j < columns - 1; j++)
    {
        a->real = B_t_B_w[j].real;
        a->imag = B_t_B_w[j].imag;
        cblas_ccopy(rows, v, 1, &B_t_B_w[j], columns);
        cblas_cscal(rows, a, &B_t_B_w[j], columns);
    }

    free(v);
    free(a);

    // Row-by-Row FFT
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim3 = NULL;
    status = DftiCreateDescriptor(&desc_handle_dim3, DFTI_DOUBLE,
                                  DFTI_COMPLEX, 1, columns);

    stride[1] = 1;

    status = DftiSetValue(desc_handle_dim3, DFTI_NUMBER_OF_TRANSFORMS, rows);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_DISTANCE, columns);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim3);
    status = DftiComputeForward(desc_handle_dim3, B_t_B_w);
    status = DftiFreeDescriptor(&desc_handle_dim3);
}

void compute_zsmooth_component_S(MKL_Complex16 *B_S, size_t rows, size_t columns)
{
    if (B_S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    MKL_Complex16 aux = {B_S[0].real, B_S[0].imag};

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            float denom = (2.0 * cos(2.0 * PI * i / rows) + 2.0 * cos(2.0 * PI * j / columns) - 4.0);
            B_S[j + i * columns].real /= denom;
            B_S[j + i * columns].imag /= denom;
        }
    }

    B_S[0].real = aux.real;
    B_S[0].imag = aux.imag;
}

void compute_zsmooth_component_S_2(MKL_Complex16 *B_S, size_t rows, size_t columns)
{
    if (B_S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    MKL_Complex16 aux = {B_S[0].real, B_S[0].imag};

    double *v1 = NULL;
    double *v2 = NULL;
    
    init_dvector(&v1, rows); 
    init_dvector(&v2, columns); 

    for (int i = 0; i < rows; i++)
    {
        v1[i] = (2.0 * PI * i) / rows;
    }
    
    for (int i = 0; i < columns; i++)
    {
        v2[i] = (2.0 * PI * (i)) / columns;
    }

    double *cos_1 = NULL;
    double *cos_2 = NULL;
    init_dvector(&cos_1, rows); 
    init_dvector(&cos_2, columns);
    vdCos(rows, v1, cos_1);
    vdCos(columns, v2, cos_2);
    cblas_dscal(rows, 2.0, cos_1, 1);
    cblas_dscal(columns, 2.0, cos_2, 1);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            float denom = cos_1[i] + cos_2[j] - 4.0;
            B_S[i + j * rows].real /= denom;
            B_S[i + j * rows].imag /= denom;
        }
    }

    B_S[0].real = aux.real;
    B_S[0].imag = aux.imag;
}

void compute_zperiodic_component_P(MKL_Complex16 *I_w, MKL_Complex16 *S, size_t rows, size_t columns)
{
    if (I_w == NULL || S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }
    
    vzSub(rows * columns, I_w, S, I_w);
}

void compute_zfftshift(MKL_Complex16 *vector, size_t rows, size_t columns) {
    size_t half_rows = rows / 2;
    size_t half_columns = columns / 2;
    
    // Troca os quadrantes (1,4) e (2,3)
    for (size_t i = 0; i < half_rows; i++) {
        for (size_t j = 0; j < half_columns; j++) {
            // Índices dos elementos a serem trocados
            size_t index1 = i * columns + j;
            size_t index2 = (i + half_rows) * columns + (j + half_columns);
            size_t index3 = i * columns + (j + half_columns);
            size_t index4 = (i + half_rows) * columns + j;

            // Swap (1,4)
            MKL_Complex16 temp = vector[index1];
            vector[index1] = vector[index2];
            vector[index2] = temp;
            
            // Swap (2,3)
            temp = vector[index3];
            vector[index3] = vector[index4];
            vector[index4] = temp;
        }
    }
}

void normalize_zvector(MKL_Complex16 *zvector, size_t size){
    double *magnitudes = NULL;
    init_dvector(&magnitudes, size);
    
    vzAbs(size, zvector, magnitudes);

    for (int i = 0; i < size; i++) {
        zvector[i].real = zvector[i].real / magnitudes[i];
        zvector[i].imag = zvector[i].imag / magnitudes[i];
    }
}