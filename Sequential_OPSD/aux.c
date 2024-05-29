#include "aux.h"

int init(MKL_Complex8 **matrix, size_t rows, size_t columns)
{
    (*matrix) = (MKL_Complex8 *)malloc(rows * columns * sizeof(MKL_Complex8));
    if ((*matrix) == NULL)
    {
        printf("Allocation error!\n");
        free(*matrix);
        return -1;
    }

    return 0;
}

void free_matrix(MKL_Complex8 *matrix)
{
    if (matrix == NULL)
    {
        return;
    }

    free(matrix);
}

void show_matrix(MKL_Complex8 *matrix, size_t rows, size_t columns)
{
    printf("\nMatriz:\n");
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            printf("(%.2f, %.2f) ", matrix[i + j * rows].real, matrix[i + j * rows].imag);
        }
        printf("\n");
    }
    printf("\n");
}

void show_ram_allocation(size_t rows, size_t columns)
{
    printf("%.2f GB\n", rows * columns * sizeof(MKL_Complex8) / (1024.0f * 1024.0f * 1024.0f));
}

void fill(MKL_Complex8 *matrix, size_t rows, size_t columns, unsigned int seed)
{
    srand(seed);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            float random_real = (float)rand() / RAND_MAX_F;
            float scaled_real = random_real * 4 + 1;

            // float random_complex = (float)rand() / RAND_MAX_F;
            // float scaled_complex = random_real * 3 + 1;

            matrix[i + j * rows].real = scaled_real;
            matrix[i + j * rows].imag = 0.0f;
        }
    }
}

void compute_fft2D_column_row(MKL_Complex8 *I_t_I_w, size_t rows, size_t columns)
{
    double start = omp_get_wtime();

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim2 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, rows);
    status = DftiCreateDescriptor(&desc_handle_dim2, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, columns);

    MKL_LONG stride[2] = {0, rows};

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, columns);
    status = DftiSetValue(desc_handle_dim1, DFTI_INPUT_DISTANCE, rows);
    status = DftiSetValue(desc_handle_dim1, DFTI_OUTPUT_DISTANCE, rows);
    status = DftiCommitDescriptor(desc_handle_dim1);

    status = DftiSetValue(desc_handle_dim2, DFTI_NUMBER_OF_TRANSFORMS, rows);
    status = DftiSetValue(desc_handle_dim2, DFTI_INPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim2, DFTI_OUTPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim2, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim2, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim2);


    // Realizar as transformações FFT
    status = DftiComputeForward(desc_handle_dim1, I_t_I_w);
    status = DftiComputeForward(desc_handle_dim2, I_t_I_w);

    status = DftiFreeDescriptor(&desc_handle_dim1);
    status = DftiFreeDescriptor(&desc_handle_dim2);

    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("Tempo gasto no step A: %f s\n", time_spent);
}

void compute_periodic_border_B(MKL_Complex8 *I_t, MKL_Complex8 *B_t, size_t rows, size_t columns)
{
    double start = omp_get_wtime();

    B_t[0].real = I_t[rows - 1].real - 2 * I_t[0].real + I_t[(columns - 1) * rows].real;
    B_t[0].imag = I_t[rows - 1].imag - 2 * I_t[0].imag + I_t[(columns - 1) * rows].imag;

    B_t[rows - 1].real = I_t[0].real - 2 * I_t[rows - 1].real + I_t[rows - 1 + (columns - 1) * rows].real;
    B_t[rows - 1].imag = I_t[0].imag - 2 * I_t[rows - 1].imag + I_t[rows - 1 + (columns - 1) * rows].imag;

    B_t[(columns - 1) * rows].real = I_t[0].real - 2 * I_t[(columns - 1) * rows].real + I_t[rows - 1 + (columns - 1) * rows].real;
    B_t[(columns - 1) * rows].imag = I_t[0].imag - 2 * I_t[(columns - 1) * rows].imag + I_t[rows - 1 + (columns - 1) * rows].imag;

    B_t[rows - 1 + (columns - 1) * rows].real = I_t[rows - 1].real - 2 * I_t[rows - 1 + (columns - 1) * rows].real + I_t[(columns - 1) * rows].real;
    B_t[rows - 1 + (columns - 1) * rows].imag = I_t[rows - 1].imag - 2 * I_t[rows - 1 + (columns - 1) * rows].imag + I_t[(columns - 1) * rows].imag;

    // linha 0 e rows-1
    for (int j = 1; j < columns - 1; j++)
    {
        B_t[j * rows].real = I_t[j * rows + rows - 1].real - I_t[j * rows].real;
        B_t[j * rows].imag = I_t[j * rows + rows - 1].imag - I_t[j * rows].imag;

        B_t[j * rows + rows - 1].real = B_t[j * rows].real * (-1);
        B_t[j * rows + rows - 1].imag = B_t[j * rows].imag * (-1);
    }

    // coluna 0 e columns-1
    for (int i = 1; i < rows - 1; i++)
    {
        B_t[i].real = I_t[i + rows * (columns - 1)].real - I_t[i].real;
        B_t[i].imag = I_t[i + rows * (columns - 1)].imag - I_t[i].imag;

        B_t[i + rows * (columns - 1)].real = B_t[i].real * (-1);
        B_t[i + rows * (columns - 1)].imag = B_t[i].imag * (-1);
    }

    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("Tempo gasto no step B: %f s\n", time_spent);
}

void compute_fft2D_of_B(MKL_Complex8 *B_t_B_w, size_t rows, size_t columns)
{
    double start = omp_get_wtime();

    // Column-one FFT
    MKL_Complex8 *a = (MKL_Complex8 *)malloc(sizeof(MKL_Complex8));
    a->real = B_t_B_w[0].real + B_t_B_w[(columns - 1) * rows].real;
    a->imag = B_t_B_w[0].imag + B_t_B_w[(columns - 1) * rows].imag;

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, rows);

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, 1);
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeForward(desc_handle_dim1, B_t_B_w);
    status = DftiFreeDescriptor(&desc_handle_dim1);

    // Calcular cada elemento de v usando a fórmula W^k = exp(-i * 2 * PI * k / M)
    MKL_Complex8 *v = (MKL_Complex8 *)malloc(rows * sizeof(MKL_Complex8));
    if (v == NULL)
    {
        printf("Erro ao alocar memória\n");
        return;
    }

    for (int k = 1; k < rows; k++)
    {
        float theta = -2.0f * PI * (rows - k) / rows;
        v[k].real = 1.0f - cosf(theta);
        v[k].imag = sinf(theta) * (-1);
    }

    cblas_ccopy(rows, B_t_B_w, 1, &B_t_B_w[(columns - 1) * rows], 1);
    cblas_csscal(rows, -1.0f, &B_t_B_w[(columns - 1) * rows], 1);
    cblas_caxpy(rows, a, v, 1, &B_t_B_w[(columns - 1) * rows], 1);

    for (int j = 1; j < columns - 1; j++)
    {
        a->real = B_t_B_w[j * rows].real;
        a->imag = B_t_B_w[j * rows].imag;
        cblas_ccopy(rows, v, 1, &B_t_B_w[j * rows], 1);
        cblas_cscal(rows, a, &B_t_B_w[j * rows], 1);
    }

    free(v);
    free(a);

    // Row-by-Row FFT
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim3 = NULL;
    status = DftiCreateDescriptor(&desc_handle_dim3, DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, columns);

    MKL_LONG stride[2] = {0, rows};

    status = DftiSetValue(desc_handle_dim3, DFTI_NUMBER_OF_TRANSFORMS, rows);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim3);
    status = DftiComputeForward(desc_handle_dim3, B_t_B_w);
    status = DftiFreeDescriptor(&desc_handle_dim3);

    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("Tempo gasto no step C: %f s\n", time_spent);
}

void compute_smooth_component_S(MKL_Complex8 *B_S, size_t rows, size_t columns)
{
    double start = omp_get_wtime();

    MKL_Complex8 aux = {B_S[0].real, B_S[0].imag};

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            B_S[i + j * rows].real /= (2.0f * cosf(2.0f * PI * i / rows) + 2.0f * cosf(2.0f * PI * j / columns) - 4.0f);
            B_S[i + j * rows].imag /= (2.0f * cosf(2.0f * PI * i / rows) + 2.0f * cosf(2.0f * PI * j / columns) - 4.0f);
        }
    }

    B_S[0].real = aux.real;
    B_S[0].imag = aux.imag;

    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("Tempo gasto no step D: %f s\n", time_spent);
}

void compute_periodic_component_P(MKL_Complex8 *I_w, MKL_Complex8 *S, size_t rows, size_t columns)
{
    double start = omp_get_wtime();
    vcSub(rows * columns, I_w, S, I_w);
    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("Tempo gasto no step E: %f s\n", time_spent);
}
