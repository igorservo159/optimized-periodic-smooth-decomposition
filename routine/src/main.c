#include "../include/utils.h"
#include "../include/fourier.h"

int main(int argc, char const *argv[])
{
    if (argc != 9)
    {
        printf("Use: %s <rows> <columns> <routine> <precision> <save_vectors> <input> <directory> <seed>\n", argv[0]);
        return -1;
    }

    size_t rows = atoi(argv[1]), columns = atoi(argv[2]), size = rows * columns;
    const char *ROUTINE = argv[3], *PRECISION = argv[4], *SAVE_VECTORS = argv[5], *INPUT = argv[6], *DIR = argv[7];
    int seed = atoi(argv[8]);

    int err_code = check_args(argv[0], ROUTINE, PRECISION, SAVE_VECTORS, INPUT);
    if (err_code)
        return err_code;

    char filepath[1024];    
   
    snprintf(filepath, sizeof(filepath), "../bin/%s", DIR);

    ensure_directory_exists(filepath);

    //CONST PATHS
    const char *IMAGE = "image.bin";
    const char *SHIFTED_SPECTRUM = "shifted_spectrum.bin";
    const char *SHIFTED_SMOOTH_SPECTRUM = "shifted_smooth_spectrum.bin";
    const char *SHIFTED_PERIODIC_SPECTRUM = "shifted_periodic_spectrum.bin";
    const char *SPECTRUM = "spectrum.bin";
    const char *SMOOTH_SPECTRUM = "smooth_spectrum.bin";
    const char *PERIODIC_SPECTRUM = "periodic_spectrum.bin";
    const char *FILTERED_IMAGE = "filtered_image.bin";

    if (!strcmp(ROUTINE, "ccr"))
    {
        if (!strcmp(PRECISION, "single"))
        {
            MKL_Complex8 *image = NULL;
            init_cvector(&image, size);

            MKL_Complex8 *B_t = NULL;
            init_cvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, IMAGE);
                read_fvector_bin(filepath, aux, size);
                copy_fvector_to_cvector(image, aux, size);
                free_fvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_cmatrix(image, rows, columns, seed);
            }           

            compute_cperiodic_border_B(image, B_t, rows, columns);
            compute_cfft2d(image, rows, columns);

            compute_cfft2d_of_border_B(B_t, rows, columns);
            compute_csmooth_component_S(B_t, rows, columns);

            compute_cperiodic_component_P(image, B_t, rows, columns);

            compute_cifft2d(image, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                copy_cvector_to_real_fvector(image, aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, FILTERED_IMAGE);
                save_fvector_on_bin(filepath, aux, size);
                free_fvector(aux);
            }
            
            free_cvector(B_t);
            free_cvector(image);
        }
        else if (!strcmp(PRECISION, "double"))
        {
            MKL_Complex16 *image = NULL;
            init_zvector(&image, size);

            MKL_Complex16 *B_t = NULL;
            init_zvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, IMAGE);
                read_dvector_bin(filepath, aux, size);
                copy_dvector_to_zvector(image, aux, size);
                free_dvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_zmatrix(image, rows, columns, seed);
            }

            compute_zperiodic_border_B(image, B_t, rows, columns);
            compute_zfft2d(image, rows, columns);

            compute_zfft2d_of_border_B(B_t, rows, columns);
            compute_zsmooth_component_S(B_t, rows, columns);

            compute_zperiodic_component_P(image, B_t, rows, columns);

            compute_zifft2d(image, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                copy_zvector_to_real_dvector(image, aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, FILTERED_IMAGE);
                save_dvector_on_bin(filepath, aux, size);
                free_dvector(aux);
            }

            free_zvector(B_t);
            free_zvector(image);
        }
    }
    else if (!strcmp(ROUTINE, "cts"))
    {
        if (!strcmp(PRECISION, "single"))
        {
            MKL_Complex8 *image = NULL;
            init_cvector(&image, size);

            MKL_Complex8 *B_t = NULL;
            init_cvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, IMAGE);
                read_fvector_bin(filepath, aux, size);
                copy_fvector_to_cvector(image, aux, size);
                free_fvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_cmatrix(image, rows, columns, seed);
            }

            compute_cperiodic_border_B(image, B_t, rows, columns);
            compute_cfft2d(image, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SPECTRUM);
                save_cvector_on_bin(filepath, image, size);
            }

            compute_cfft2d_of_border_B(B_t, rows, columns);
            compute_csmooth_component_S(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SMOOTH_SPECTRUM);
                save_cvector_on_bin(filepath, B_t, size);
            }

            compute_cperiodic_component_P(image, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, PERIODIC_SPECTRUM);
                save_cvector_on_bin(filepath, image, size);
            }

            free_cvector(B_t);
            free_cvector(image);
        }
        else if (!strcmp(PRECISION, "double"))
        {
            MKL_Complex16 *image = NULL;
            init_zvector(&image, size);

            MKL_Complex16 *B_t = NULL;
            init_zvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, IMAGE);
                read_dvector_bin(filepath, aux, size);
                copy_dvector_to_zvector(image, aux, size);
                free_dvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_zmatrix(image, rows, columns, seed);
            }

            compute_zperiodic_border_B(image, B_t, rows, columns);
            compute_zfft2d(image, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SPECTRUM);
                save_zvector_on_bin(filepath, image, size);
            }

            compute_zfft2d_of_border_B(B_t, rows, columns);
            compute_zsmooth_component_S(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SMOOTH_SPECTRUM);
                save_zvector_on_bin(filepath, B_t, size);
            }

            compute_zperiodic_component_P(image, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, PERIODIC_SPECTRUM);
                save_zvector_on_bin(filepath, image, size);
            }

            free_zvector(B_t);
            free_zvector(image);
        }
    }
    else if (!strcmp(ROUTINE, "css"))
    {
        if (!strcmp(PRECISION, "single"))
        {
            MKL_Complex8 *image = NULL;
            init_cvector(&image, size);

            MKL_Complex8 *B_t = NULL;
            init_cvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                float *aux = NULL;
                init_fvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, IMAGE);
                read_fvector_bin(filepath, aux, size);
                copy_fvector_to_cvector(image, aux, size);
                free_fvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_cmatrix(image, rows, columns, seed);
            }

            compute_cperiodic_border_B(image, B_t, rows, columns);
            compute_cfft2d(image, rows, columns);
            compute_cfftshift(image, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SHIFTED_SPECTRUM);
                save_cvector_on_bin(filepath, image, size);
            }

            compute_cfft2d_of_border_B(B_t, rows, columns);
            compute_csmooth_component_S(B_t, rows, columns);
            compute_cfftshift(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SHIFTED_SMOOTH_SPECTRUM);
                save_cvector_on_bin(filepath, B_t, size);
            }

            compute_cperiodic_component_P(image, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SHIFTED_PERIODIC_SPECTRUM);
                save_cvector_on_bin(filepath, image, size);
            }

            free_cvector(B_t);
            free_cvector(image);
        }
        else if (!strcmp(PRECISION, "double"))
        {
            MKL_Complex16 *image = NULL;
            init_zvector(&image, size);

            MKL_Complex16 *B_t = NULL;
            init_zvector(&B_t, size);

            if (!strcmp(INPUT, "rb"))
            {
                double *aux = NULL;
                init_dvector(&aux, size);
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, IMAGE);
                read_dvector_bin(filepath, aux, size);
                copy_dvector_to_zvector(image, aux, size);
                free_dvector(aux);
            }
            else if (!strcmp(INPUT, "fm"))
            {
                fill_zmatrix(image, rows, columns, seed);
            }

            compute_zperiodic_border_B(image, B_t, rows, columns);
            compute_zfft2d(image, rows, columns);
            compute_zfftshift(image, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SHIFTED_SPECTRUM);
                save_zvector_on_bin(filepath, image, size);
            }

            compute_zfft2d_of_border_B(B_t, rows, columns);
            compute_zsmooth_component_S(B_t, rows, columns);
            compute_zfftshift(B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SHIFTED_SMOOTH_SPECTRUM);
                save_zvector_on_bin(filepath, B_t, size);
            }

            compute_zperiodic_component_P(image, B_t, rows, columns);

            if (!strcmp(SAVE_VECTORS, "yes"))
            {
                snprintf(filepath, sizeof(filepath), "../bin/%s/%s", DIR, SHIFTED_PERIODIC_SPECTRUM);
                save_zvector_on_bin(filepath, image, size);
            }

            free_zvector(B_t);
            free_zvector(image);
        }
    }
    return 0;
}
