#include <mkl.h>
#include <mkl_dfti.h>
#include <stdio.h>
#include <omp.h>

#define PI 3.14159265358979323846

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

void compute_csmooth_component_S_2(MKL_Complex8 *B_S, size_t rows, size_t columns)
{
    if (B_S == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    MKL_Complex8 aux = {B_S[0].real, B_S[0].imag};

    float *v1 = (float *)malloc(rows * sizeof(float));
    float *v2 = (float *)malloc(columns * sizeof(float));

    for (int i = 0; i < rows; i++)
    {
        v1[i] = (2.0f * PI * i) / rows;
    }
    
    for (int i = 0; i < columns; i++)
    {
        v2[i] = (2.0f * PI * (i)) / columns;
    }

    float *cos_1 = (float *)malloc(rows * sizeof(float));
    float *cos_2 = (float *)malloc(columns * sizeof(float));
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

int main(){
    MKL_Complex8 *vector = NULL;
    vector = (MKL_Complex8 *)malloc(10000 * 10000 * sizeof(MKL_Complex8));

    #pragma parallel for collapse(2)
    for(int i = 0; i < 10000; i++){
        for(int j = 0; j < 10000; j++){
            vector[j + i * 10000].real = 1.0f;
            vector[j + i * 10000].imag = 1.0f;
        }
    }

    compute_cfft2d(vector, 10000, 10000);

}