#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <mkl.h>
#include <omp.h>

#define RAND_MAX_F 2147483647.0f

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        printf("Use: %s <M> <N>\n", argv[0]);
        return -1;
    }

    int M = atoi(argv[1]);
    int N = atoi(argv[2]);

    MKL_Complex8 *data = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));
    if (data == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    //srand(time(NULL));

    srand(1);

    //printf("\nValores da Imagem:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            float random_real = (float)rand() / RAND_MAX_F;
            float scaled_real = random_real * 4 + 1;

            data[i + j*M].real = scaled_real;
            data[i + j*M].imag = 0.0f;

            //printf("(%.2f, %.2f) ", data[i + j*M].real, data[i + j*M].imag);
        }
        //printf("\n");
    }
    printf("\n");


    //Fazendo o Vetor B

    MKL_Complex8 *VectorB = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));

    double start = omp_get_wtime();

    int i = 0;

    //Computando as quinas
    VectorB[0].real = data[M-1].real - 2*data[0].real + data[(N-1)*M].real;
    VectorB[0].imag = data[M-1].imag - 2*data[0].imag + data[(N-1)*M].imag;

    VectorB[M-1].real = data[0].real - 2*data[M-1].real + data[M-1 + (N-1)*M].real;
    VectorB[M-1].imag = data[0].imag - 2*data[M-1].imag + data[M-1 + (N-1)*M].imag;

    VectorB[(N-1)*M].real = data[0].real - 2*data[(N-1)*M].real + data[M-1 + (N-1)*M].real;
    VectorB[(N-1)*M].imag = data[0].imag - 2*data[(N-1)*M].imag + data[M-1 + (N-1)*M].imag;

    VectorB[M-1 + (N-1)*M].real = data[M-1].real - 2*data[M-1 + (N-1)*M].real + data[(N-1)*M].real;
    VectorB[M-1 + (N-1)*M].imag = data[M-1].imag - 2*data[M-1 + (N-1)*M].imag + data[(N-1)*M].imag;


    //linha 0 e M-1
    for(int j = 1; j < N-1; j++){
        VectorB[j*M].real = data[j*M + M-1].real - data[j*M].real;
        VectorB[j*M].imag = data[j*M + M-1].imag - data[j*M].imag;

        VectorB[j*M + M-1].real = VectorB[j*M].real*(-1);
        VectorB[j*M + M-1].imag = VectorB[j*M].imag*(-1);
    }

    //coluna 0 e N-1
    for(int i = 1; i < M-1; i++){
        VectorB[i].real = data[i + M*(N-1)].real - data[i].real;
        VectorB[i].imag = data[i + M*(N-1)].imag - data[i].imag;

        VectorB[i + M*(N-1)].real = VectorB[i].real*(-1);
        VectorB[i + M*(N-1)].imag = VectorB[i].imag*(-1);
    }

    for(int i = 1; i < M-1; i++){
        for(int j = 1; j < N-1; j++){
            VectorB[i + j*M].real = 0.0f;
            VectorB[i + j*M].imag = 0.0f;
        }
    }

    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("%f s\n", time_spent);
    
    /*
    printf("\nValores de B:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            printf("(%.2f, %.2f) ", VectorB[i + j*M].real, VectorB[i + j*M].imag);
        }
        printf("\n");
    }
    printf("\n");
    */

    free(data);
    free(VectorB);

    return 0;
}
