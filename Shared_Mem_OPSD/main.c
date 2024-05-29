#include "aux.h"

#define RAND_MAX_F 2147483647.0f

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        printf("Use: %s <M> <N>\n", argv[0]);
        return -1;
    }

    size_t rows = atoi(argv[1]);
    size_t columns = atoi(argv[2]);

    MKL_Complex8 *I_t = NULL;
    init(&I_t, rows, columns);
    fill(I_t, rows, columns, 1);

    show_matrix(I_t, rows, columns);

    MKL_Complex8 *B_t = NULL;
    init(&B_t, rows, columns);

    double start_ = omp_get_wtime();

    compute_periodic_border_B(I_t, B_t, rows, columns);
    compute_fft2D_column_row(I_t, rows, columns);
    show_matrix(I_t, rows, columns);
    //compute_fft2D_column_row(B_t, rows, columns);
    compute_fft2D_of_B(B_t, rows, columns);
    compute_smooth_component_S(B_t, rows, columns);
    compute_periodic_component_P(I_t, B_t, rows, columns);

    double end_ = omp_get_wtime();
    double time_spent_ = (end_ - start_);
    printf("Tempo total gasto: %f s\n", time_spent_);


    free_matrix(B_t);
    free_matrix(I_t);
}