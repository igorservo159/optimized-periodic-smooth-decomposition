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

    size_t M = atoi(argv[1]);
    size_t N = atoi(argv[2]);

    MKL_Complex8 *data = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));
    if (data == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    printf("%.2f GB", M * N * sizeof(MKL_Complex8)/(1024.0f * 1024.0f * 1024.0f));


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


    MKL_Complex8 *VectorB = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));
    MKL_Complex8 *aux = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));

    double start = omp_get_wtime();

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            if(i == 0 || i == M-1){
                VectorB[i + j*M].real = data[M-1-i + j*M].real - data[i + j*M].real;
                VectorB[i + j*M].imag = data[M-1-i + j*M].imag - data[i + j*M].imag;
            } else{
                VectorB[i + j*M].real = 0.0f;
                VectorB[i + j*M].imag = 0.0f;
            }
            if(j == 0 || j == N-1){
                aux[i + j*M].real = data[i + (N-1-j)*M].real - data[i + j*M].real;
                aux[i + j*M].imag = data[i + (N-1-j)*M].imag - data[i + j*M].imag;
            } else{
                aux[i + j*M].real = 0.0f;
                aux[i + j*M].imag = 0.0f;
            }
        }
    }

    /*
    printf("\nValores de R:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            printf("(%.2f, %.2f) ", VectorB[i + j*M].real, VectorB[i + j*M].imag);
        }
        printf("\n");
    }
    printf("\n");

    printf("\nValores de C:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            printf("(%.2f, %.2f) ", aux[i + j*M].real, aux[i + j*M].imag);
        }
        printf("\n");
    }
    printf("\n");

    */

    vcAdd(N*M, (MKL_Complex8 *)aux, (MKL_Complex8 *)VectorB, (MKL_Complex8 *)VectorB);

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
    free(aux);
    free(VectorB);

    return 0;
}
