#include "../include/utils.h"
#include "../include/fourier.h"

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        printf("Use: %s <M> <N>\n", argv[0]);
        return -1;
    }

    size_t rows = atoi(argv[1]);
    size_t columns = atoi(argv[2]);
    size_t size = rows*columns;

    MKL_Complex8 *I_t = NULL;
    init_cvector(&I_t, size);

    MKL_Complex8 *B_t = NULL;
    init_cvector(&B_t, size);

    float *aux = NULL;
    init_fvector(&aux, size);
    read_fvector_bin("../bin/example/data.bin", aux, size);
    copy_fvector_to_cvector(I_t, aux, size);
    
    double start = omp_get_wtime();


    compute_periodic_border_B(I_t, B_t, rows, columns);
    compute_fft2d(I_t, rows, columns);

    //compute_fftshift(I_t, rows, columns);
    //save_cvector_on_bin("../bin/example/spectrum_shifted.bin", I_t, size);

    save_cvector_on_bin("../bin/example/spectrum.bin", I_t, size);

    compute_fft2d_of_border_B(B_t, rows, columns);
    compute_smooth_component_S(B_t, rows, columns);
  
    //compute_fftshift(B_t, rows, columns);
    //save_cvector_on_bin("../bin/example/smooth_shifted.bin", B_t, size);
   
    save_cvector_on_bin("../bin/example/smooth.bin", B_t, size);

    compute_periodic_component_P(I_t, B_t, rows, columns);

    save_cvector_on_bin("../bin/example/periodic.bin", I_t, size);
    
    //save_cvector_on_bin("../bin/example/periodic_shifted.bin", I_t, size);

    compute_ifft2d(I_t, rows, columns);

    copy_cvector_to_real_fvector(I_t, aux, size);
    save_fvector_on_bin("../bin/example/data_filtered.bin", aux, size);


    double end = omp_get_wtime();
    double time_spent = (end - start);
    printf("Tempo total gasto: %f s\n", time_spent);

    free_fvector(aux);
    free_cvector(B_t);
    free_cvector(I_t);

    return 0;
}
