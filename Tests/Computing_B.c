#include "stdio.h"
#include "stdlib.h"
#include "complex.h"
#include "time.h"
#include "mkl.h"

#define RAND_MAX_F 2147483647.0f

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        printf("Use: %s <M> <N>\n", argv[0]);
        return -1;
    }

    int M = atoi(argv[1]);
    int N = atoi(argv[2]);

    // Alocar memória dinamicamente para o array data
    _Complex float *data = (_Complex float *)malloc(M * N * sizeof(_Complex float));
    if (data == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    // Seed para gerar números aleatórios diferentes em cada execução
    //srand(time(NULL));

    srand(1);

    printf("\nValores de data antes da FFT:\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            float random_real = (float)rand() / RAND_MAX_F;
            // Normaliza o número para o intervalo [1, 5]
            float scaled_real = random_real * 4 + 1;

            //data[j + i*N] = scaled_real + 0.0f * I; // Row-Major Layout
            data[i + j*M] = scaled_real + 0.0f * I; // Column-Major Layout

            //printf("(%.2f, %.2f) ", crealf(data[j + i*N]), cimagf(data[j + i*N])); //Row-Major Layout
            printf("(%.2f, %.2f) ", crealf(data[i + j*M]), cimagf(data[i + j*M])); //Column-Major Layout
        }
        printf("\n");
    }
    printf("\n");

    // Liberar a memória alocada dinamicamente
    free(data);

    return 0;
}
