#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <math.h>
#include <mkl.h>
#include <omp.h>

#define RAND_MAX_F 2147483647.0f
#define PI 3.14159265358979323846

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        printf("Use: %s <M> <N>\n", argv[0]);
        return -1;
    }

    size_t M = atoi(argv[1]);
    size_t N = atoi(argv[2]);
    

    MKL_Complex8 *data = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));
    if (data == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    printf("%.2f GB\n", M * N * sizeof(MKL_Complex8)/(1024.0f * 1024.0f * 1024.0f));

    //srand(time(NULL));

    srand(1);

    //Exemplo de matriz B

    //Computando as quinas
    data[0].real = -2.01f;
    data[0].imag = 0.0f;

    data[M-1].real = 0.29f;
    data[M-1].imag = 0.0f;

    data[(N-1)*M].real = 3.40f;
    data[(N-1)*M].imag = 0.0f;

    data[M-1 + (N-1)*M].real = -1.68f;
    data[M-1 + (N-1)*M].imag = 0.0f;


    //linha 0 e M-1
    for(int j = 1; j < N-1; j++){
        float random_real = (float)rand() / RAND_MAX_F;
        float scaled_real = random_real * 4 + 1;

        data[j*M].real = scaled_real;
        data[j*M].imag = 0.0f;

        data[j*M + M-1].real = data[j*M].real*(-1);
        data[j*M + M-1].imag = data[j*M].imag*(-1);
    }

    //coluna 0 e N-1
    for(int i = 1; i < M-1; i++){
        float random_real = (float)rand() / RAND_MAX_F;
        float scaled_real = random_real * 4 + 1;

        data[i].real = scaled_real;
        data[i].imag = 0.0f;

        data[i + M*(N-1)].real = data[i].real*(-1);
        data[i + M*(N-1)].imag = data[i].imag*(-1);
    }

    for(int i = 1; i < M-1; i++){
        for(int j = 1; j < N-1; j++){
            data[i + j*M].real = 0.0f;
            data[i + j*M].imag = 0.0f;
        }
    }

    /*
    printf("\nValores de data antes da FFT:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            printf("(%.2f, %.2f) ", data[i + j*M].real, data[i + j*M].imag); 
        }
        printf("\n");
    }
    printf("\n");
    */   

    double start = omp_get_wtime();

    //Column-by-column FFT

    MKL_Complex8 *a = (MKL_Complex8 *)malloc(sizeof(MKL_Complex8));
    a->real = data[0].real + data[(N-1)*M].real;
    a->imag = data[0].imag + data[(N-1)*M].imag;

    DFTI_DESCRIPTOR_HANDLE desc_handle_dim1 = NULL;
    MKL_LONG status;

    status = DftiCreateDescriptor(&desc_handle_dim1, DFTI_SINGLE,
                              DFTI_COMPLEX, 1, M);

    status = DftiSetValue(desc_handle_dim1, DFTI_NUMBER_OF_TRANSFORMS, 1);
    status = DftiCommitDescriptor(desc_handle_dim1);
    status = DftiComputeForward(desc_handle_dim1, data);
    status = DftiFreeDescriptor(&desc_handle_dim1);



    // Calcular cada elemento de v usando a fórmula W^k = exp(-i * 2 * PI * k / M)
    MKL_Complex8 *v = (MKL_Complex8 *)malloc(M * sizeof(MKL_Complex8));
    if (v == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    /*
    float *t = (float *)malloc(M * sizeof(float));
    for (int k = 0; k < M; k++) {
        t[k] = -2.0f * PI * k / M;
    }
    vcCIS(M, t, v);
    */

    v[0].real = 0;
    v[0].imag = 0;
    for(int k = 1; k < M; k++) {
        float theta = -2.0f * PI * (M-k) / M;
        v[k].real = 1 - cosf(theta);
        v[k].imag = sinf(theta)*(-1);
    }    
    
    cblas_ccopy(M, data, 1, &data[(N-1)*M],1);
    cblas_csscal(M, -1, &data[(N-1)*M], 1);
    cblas_caxpy(M, a, v, 1, &data[(N-1)*M], 1);
    
    for(int j = 1; j < N-1; j++){
        a->real = data[j*M].real;
        a->imag = data[j*M].imag;
        cblas_ccopy(M, v, 1, &data[j*M],1);
        cblas_cscal(M, a, &data[j*M], 1);
    }


    free(v);
    free(a);
        
    /*
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim2 = NULL;
    status = DftiCreateDescriptor(&desc_handle_dim2, DFTI_SINGLE,
                              DFTI_COMPLEX, 1, M);
    status = DftiSetValue(desc_handle_dim2, DFTI_NUMBER_OF_TRANSFORMS, 1);
    status = DftiCommitDescriptor(desc_handle_dim2);
    status = DftiComputeForward(desc_handle_dim2, &data[(N-1)*M]);
    */

    //Row-by-Row FFT
    DFTI_DESCRIPTOR_HANDLE desc_handle_dim3 = NULL;
    status = DftiCreateDescriptor(&desc_handle_dim3, DFTI_SINGLE,
                              DFTI_COMPLEX, 1, N);

    MKL_LONG stride[2] = {0, M};

    status = DftiSetValue(desc_handle_dim3, DFTI_NUMBER_OF_TRANSFORMS, M);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_DISTANCE, 1);
    status = DftiSetValue(desc_handle_dim3, DFTI_INPUT_STRIDES, stride);
    status = DftiSetValue(desc_handle_dim3, DFTI_OUTPUT_STRIDES, stride);
    status = DftiCommitDescriptor(desc_handle_dim3);
    status = DftiComputeForward(desc_handle_dim3, data);
    status = DftiFreeDescriptor(&desc_handle_dim3);


    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("%f s\n", time_spent);


    printf("\nValores de data depois da FFT:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            printf("(%.2f, %.2f) ", data[i + j*M].real, data[i + j*M].imag); 
        }
        printf("\n");
    }
    printf("\n");

    free(data);
                            
}