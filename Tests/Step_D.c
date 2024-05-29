#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <math.h>
#include <mkl.h>
#include <omp.h>

#define RAND_MAX_F 2147483647.0f
#define PI 3.14159265358979323846

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        printf("Use: %s <M> <N>\n", argv[0]);
        return -1;
    }

    size_t M = atoi(argv[1]);
    size_t N = atoi(argv[2]);

    MKL_Complex8 *data = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));
    if (data == NULL)
    {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    printf("%.2f GB\n", M * N * sizeof(MKL_Complex8) / (1024.0f * 1024.0f * 1024.0f));

    // srand(time(NULL));

    srand(1);

    
    data[0].real = 0.00f;
    data[1].real = 2.78f;
    data[2].real = -3.44f;
    data[3].real = -7.38f;
    

    // printf("\nValores da Imagem:\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            //float random_real = (float)rand() / RAND_MAX_F;
            //float scaled_real = random_real * 4 + 1;

            //data[i + j * M].real = scaled_real;

            printf("(%.2f, %.2f) ", data[i + j * M].real, data[i + j * M].imag);
        }
        printf("\n");
    }
    printf("\n");

    MKL_Complex8 aux = {data[0].real, data[0].imag};
    // MKL_Complex8 *D = (MKL_Complex8 *)malloc(M * N * sizeof(MKL_Complex8));

    double start = omp_get_wtime();

    /*
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            D[i + j * M].real = (2.0f * cosf(2.0f * PI * i / M) + 2.0f * cosf(2.0f * PI * j / N) - 4.0f);
            // printf("%.2f + i%.2f ", D[i + j * M].real, D[i + j * M].imag);
        }
        // printf("\n");
    }

    vcDiv(M * N, data, D, data);
     */

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            data[i + j * M].real /= (2.0f * cosf(2.0f * PI * i / M) + 2.0f * cosf(2.0f * PI * j / N) - 4.0f);
            data[i + j * M].imag /= (2.0f * cosf(2.0f * PI * i / M) + 2.0f * cosf(2.0f * PI * j / N) - 4.0f);
        }
    }

    data[0].real = aux.real;
    data[0].imag = aux.imag;

    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("%f s\n", time_spent);

    
    printf("\nValores de D:\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("(%.2f, %.2f) ", data[i + j * M].real, data[i + j * M].imag);
        }
        printf("\n");
    }
    printf("\n");
    

    free(data);

    return 0;
}
