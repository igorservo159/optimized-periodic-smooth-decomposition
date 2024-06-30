#include "../include/utils.h"
#include "../include/fourier.h"

int main(int argc, char const *argv[])
{
    if (argc != 10)
    {
        printf("Use: %s <rows> <columns> <routine> <precision> <save_vectors> <normalize> <input> <directory> <seed>\n", argv[0]);
        return -1;
    }

    size_t rows = atoi(argv[1]), columns = atoi(argv[2]), size = rows * columns;
    const char *ROUTINE = argv[3], *PRECISION = argv[4], *SAVE_VECTORS = argv[5], *NORMALIZE = argv[6], *INPUT = argv[7], *DIR = argv[8];
    int seed = atoi(argv[9]);

    int err_code = check_args(argv[0], ROUTINE, PRECISION, SAVE_VECTORS, INPUT, NORMALIZE);
    if (err_code)
        return err_code;

    char filepath[1024];
   
    snprintf(filepath, sizeof(filepath), "../bin/%s", DIR);

    ensure_directory_exists(filepath);

    if (!strcmp(ROUTINE, "ccr"))
    {
        if (!strcmp(PRECISION, "single"))
        {
            MKL_Complex8 *I_t = NULL;
            init_cvector(&I_t, size);

            MKL_Complex8 *B_t = NULL;
            init_cvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data.bin", DIR);
                read_fvector_bin(filepath, aux, size);

                if(!strcmp(NORMALIZE, "yes")){
                    normalize_fvector(aux, size);
                    save_fvector_on_bin(filepath, aux, size);
                }

                copy_fvector_to_cvector(I_t, aux, size);
                free_fvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_cmatrix(I_t, rows, columns, seed);
            }

            MKL_Complex8 *teste = NULL;
            init_cvector(&teste, size);            

            compute_cperiodic_border_B(I_t, B_t, rows, columns);
            compute_cfft2d(I_t, rows, columns);
            cblas_ccopy(size, I_t, 1, teste, 1);
            compute_cifft2d(teste, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/spectrum.bin", DIR);
                save_cvector_on_bin(filepath, I_t, size);
            }

            compute_cfft2d_of_border_B(B_t, rows, columns);
            compute_csmooth_component_S(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/smooth.bin", DIR);
                save_cvector_on_bin(filepath, B_t, size);
            }

            compute_cperiodic_component_P(I_t, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/periodic.bin", DIR);
                save_cvector_on_bin(filepath, I_t, size);
            }

            compute_cifft2d(I_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                copy_cvector_to_real_fvector(I_t, aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data_filtered.bin", DIR);
                save_fvector_on_bin(filepath, aux, size);
                free_fvector(aux);
            }


            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/data1.bin", DIR);
                save_cvector_on_bin(filepath, teste, size);
            }

            vcSub(size, teste, I_t, teste);
            
            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/compare.bin", DIR);
                save_cvector_on_bin(filepath, teste, size);
            }
            
            free_cvector(B_t);
            free_cvector(teste);
            free_cvector(I_t);
        }
        else if (!strcmp(PRECISION, "double"))
        {
            MKL_Complex16 *I_t = NULL;
            init_zvector(&I_t, size);

            MKL_Complex16 *B_t = NULL;
            init_zvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data.bin", DIR);
                read_dvector_bin(filepath, aux, size);

                if(!strcmp(NORMALIZE, "yes")){
                    normalize_dvector(aux, size);
                    save_dvector_on_bin(filepath, aux, size);
                }
                
                copy_dvector_to_zvector(I_t, aux, size);
                free_dvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_zmatrix(I_t, rows, columns, seed);
            }

            compute_zperiodic_border_B(I_t, B_t, rows, columns);
            compute_zfft2d(I_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/spectrum.bin", DIR);
                save_zvector_on_bin(filepath, I_t, size);
            }

            compute_zfft2d_of_border_B(B_t, rows, columns);
            compute_zsmooth_component_S(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/smooth.bin", DIR);
                save_zvector_on_bin(filepath, B_t, size);
            }

            compute_zperiodic_component_P(I_t, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/periodic.bin", DIR);
                save_zvector_on_bin(filepath, I_t, size);
            }

            compute_zifft2d(I_t, rows, columns);
            //normalize_zvector(I_t, size);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                copy_zvector_to_real_dvector(I_t, aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data_filtered.bin", DIR);
                save_dvector_on_bin(filepath, aux, size);
                free_dvector(aux);
            }

            free_zvector(B_t);
            free_zvector(I_t);
        }
    }
    else if (!strcmp(ROUTINE, "cts"))
    {
        if (!strcmp(PRECISION, "single"))
        {
            MKL_Complex8 *I_t = NULL;
            init_cvector(&I_t, size);

            MKL_Complex8 *B_t = NULL;
            init_cvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data.bin", DIR);
                read_fvector_bin(filepath, aux, size);

                if(!strcmp(NORMALIZE, "yes")){
                    normalize_fvector(aux, size);
                    save_fvector_on_bin(filepath, aux, size);
                }
                
                copy_fvector_to_cvector(I_t, aux, size);
                free_fvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_cmatrix(I_t, rows, columns, seed);
            }

            compute_cperiodic_border_B(I_t, B_t, rows, columns);
            compute_cfft2d(I_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/spectrum.bin", DIR);
                save_cvector_on_bin(filepath, I_t, size);
            }

            compute_cfft2d_of_border_B(B_t, rows, columns);
            compute_csmooth_component_S(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/smooth.bin", DIR);
                save_cvector_on_bin(filepath, B_t, size);
            }

            compute_cperiodic_component_P(I_t, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/periodic.bin", DIR);
                save_cvector_on_bin(filepath, I_t, size);
            }

            free_cvector(B_t);
            free_cvector(I_t);
        }
        else if (!strcmp(PRECISION, "double"))
        {
            MKL_Complex16 *I_t = NULL;
            init_zvector(&I_t, size);

            MKL_Complex16 *B_t = NULL;
            init_zvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data.bin", DIR);
                read_dvector_bin(filepath, aux, size);

                if(!strcmp(NORMALIZE, "yes")){
                    normalize_dvector(aux, size);
                    save_dvector_on_bin(filepath, aux, size);
                }
                
                copy_dvector_to_zvector(I_t, aux, size);
                free_dvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_zmatrix(I_t, rows, columns, seed);
            }

            compute_zperiodic_border_B(I_t, B_t, rows, columns);
            compute_zfft2d(I_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/spectrum.bin", DIR);
                save_zvector_on_bin(filepath, I_t, size);
            }

            compute_zfft2d_of_border_B(B_t, rows, columns);
            compute_zsmooth_component_S(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/smooth.bin", DIR);
                save_zvector_on_bin(filepath, B_t, size);
            }

            compute_zperiodic_component_P(I_t, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/periodic.bin", DIR);
                save_zvector_on_bin(filepath, I_t, size);
            }

            free_zvector(B_t);
            free_zvector(I_t);
        }
    }
    else if (!strcmp(ROUTINE, "css"))
    {
        if (!strcmp(PRECISION, "single"))
        {
            MKL_Complex8 *I_t = NULL;
            init_cvector(&I_t, size);

            MKL_Complex8 *B_t = NULL;
            init_cvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data.bin", DIR);
                read_fvector_bin(filepath, aux, size);

                if(!strcmp(NORMALIZE, "yes")){
                    normalize_fvector(aux, size);
                    save_fvector_on_bin(filepath, aux, size);
                }
                
                copy_fvector_to_cvector(I_t, aux, size);
                free_fvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_cmatrix(I_t, rows, columns, seed);
            }

            compute_cperiodic_border_B(I_t, B_t, rows, columns);
            compute_cfft2d(I_t, rows, columns);
            compute_cfftshift(I_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/spectrum_shifted.bin", DIR);
                save_cvector_on_bin(filepath, I_t, size);
            }

            compute_cfft2d_of_border_B(B_t, rows, columns);
            compute_csmooth_component_S(B_t, rows, columns);
            compute_cfftshift(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/smooth_shifted.bin", DIR);
                save_cvector_on_bin(filepath, B_t, size);
            }

            compute_cperiodic_component_P(I_t, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/periodic_shifted.bin", DIR);
                save_cvector_on_bin(filepath, I_t, size);
            }

            free_cvector(B_t);
            free_cvector(I_t);
        }
        else if (!strcmp(PRECISION, "double"))
        {
            MKL_Complex16 *I_t = NULL;
            init_zvector(&I_t, size);

            MKL_Complex16 *B_t = NULL;
            init_zvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/data.bin", DIR);
                read_dvector_bin(filepath, aux, size);

                if(!strcmp(NORMALIZE, "yes")){
                    normalize_dvector(aux, size);
                    save_dvector_on_bin(filepath, aux, size);
                }
                
                copy_dvector_to_zvector(I_t, aux, size);
                free_dvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_zmatrix(I_t, rows, columns, seed);
            }

            compute_zperiodic_border_B(I_t, B_t, rows, columns);
            compute_zfft2d(I_t, rows, columns);
            compute_zfftshift(I_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/spectrum_shifted.bin", DIR);
                save_zvector_on_bin(filepath, I_t, size);
            }

            compute_zfft2d_of_border_B(B_t, rows, columns);
            compute_zsmooth_component_S(B_t, rows, columns);
            compute_zfftshift(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/smooth_shifted.bin", DIR);
                save_zvector_on_bin(filepath, B_t, size);
            }

            compute_zperiodic_component_P(I_t, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/periodic_shifted.bin", DIR);
                save_zvector_on_bin(filepath, I_t, size);
            }

            free_zvector(B_t);
            free_zvector(I_t);
        }
    }
    return 0;
}
