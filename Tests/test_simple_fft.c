#include "mkl_dfti.h"
#include "stdio.h"
#include "stdlib.h"
#include "complex.h"
#include "time.h"
#include "mkl.h"

#define RAND_MAX_F 2147483647.0f

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        printf("Use: %s <size_of_FFT>\n", argv[0]);
        return -1;
    }

    int size_of = atoi(argv[1]);

    // Alocar memória dinamicamente para o array r2c_data
    float *r2c_data = (float *)malloc(size_of * sizeof(float));
    if (r2c_data == NULL) {
        printf("Erro ao alocar memória\n");
        return 1;
    }

    DFTI_DESCRIPTOR_HANDLE my_desc2_handle = NULL;
    MKL_LONG status;
    clock_t start, end;
    double cpu_time_used;

    // Seed para gerar números aleatórios diferentes em cada execução
    srand(time(NULL));

    //printf("\nValores de r2c_data antes da FFT:\n");
    for(int i = 0; i < size_of; i++){
        float random_real = (float)rand() / RAND_MAX_F;
        //Normaliza o número para o intervalo [1, 5]
        float scaled_real = random_real * 4 + 1;

        r2c_data[i] = scaled_real;

        //printf("%.2f ", r2c_data[i]);
    }
    //printf("\n\n");

    status = DftiCreateDescriptor(&my_desc2_handle, DFTI_SINGLE,
                                  DFTI_REAL, 1, size_of);
    //status = DftiSetValue(my_desc2_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiCommitDescriptor(my_desc2_handle);

    start = clock();

    status = DftiComputeForward(my_desc2_handle, r2c_data);

    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    status = DftiFreeDescriptor(&my_desc2_handle);

    /*printf("\nValores de r2c_data depois da FFT:\n");
    for(int i = 0; i < size_of+2; i++){
        printf("%.2f ", r2c_data[i]);
    }

    printf("\n\n");
    */
    printf("Tempo de execução da FFT: %f segundos\n", cpu_time_used);

    // Liberar a memória alocada dinamicamente
    free(r2c_data);

    return 0;
}
