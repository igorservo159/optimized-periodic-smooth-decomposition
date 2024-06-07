#include "../include/utils.h"

void init_cvector(MKL_Complex8 **vector, size_t size)
{
    (*vector) = (MKL_Complex8 *)malloc(size * sizeof(MKL_Complex8));
    if ((*vector) == NULL)
    {
        printf("Allocation error!\n");
        free(*vector);
        return;
    }
}

void init_fvector(float **vector, size_t size)
{
    (*vector) = (float *)malloc(size * sizeof(float));
    if ((*vector) == NULL)
    {
        printf("Allocation error!\n");
        free(*vector);
        return;
    }
}

void free_cvector(MKL_Complex8 *vector)
{
    if (vector == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    free(vector);
}

void free_fvector(float *vector){
    if (vector == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    free(vector);
}

void show_cmatrix(MKL_Complex8 *matrix, size_t rows, size_t columns)
{
    if (matrix == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    printf("\nMatriz:\n");
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            printf("(%.2f, %.2f) ", matrix[j + i * columns].real, matrix[j + i * columns].imag);
        }
        printf("\n");
    }
    printf("\n");
}

void fill_cmatrix(MKL_Complex8 *matrix, size_t rows, size_t columns, unsigned int seed)
{
    if (matrix == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    srand(seed);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            float random_real = (float)rand() / RAND_MAX_F;
            float scaled_real = random_real * 4 + 1;

            // float random_complex = (float)rand() / RAND_MAX_F;
            // float scaled_complex = random_real * 3 + 1;

            // matrix[i + j * rows].real = 1.0f;
            matrix[i + j * rows].real = scaled_real;
            matrix[i + j * rows].imag = 0.0f;
        }
    }
}

void read_fvector_bin(const char *filename, float *vector, size_t size)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Error opening file");
        exit(1);
    }

    size_t elements_read = fread(vector, sizeof(float), size, fp);
    if (elements_read != size) {
        if (feof(fp)) {
            fprintf(stderr, "Error: unexpected end of file\n");
        } else if (ferror(fp)) {
            perror("Error reading file");
        }
        fclose(fp);
        exit(1);
    }

    fclose(fp);
}

void read_cvector_bin(const char *filename, MKL_Complex8 *vector, size_t size){
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Error opening file");
        exit(1);
    }

    size_t elements_read = fread(vector, sizeof(MKL_Complex8), size, fp);
    if (elements_read != size) {
        if (feof(fp)) {
            fprintf(stderr, "Error: unexpected end of file\n");
        } else if (ferror(fp)) {
            perror("Error reading file");
        }
        fclose(fp);
        exit(1);
    }

    fclose(fp);
}

void save_cvector_on_bin(const char *filename, MKL_Complex8 *vector, size_t size){
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Erro ao abrir o arquivo para escrita");
        exit(EXIT_FAILURE);
    }
    fwrite(vector, sizeof(MKL_Complex8), size, file);
    fclose(file);
}

void save_fvector_on_bin(const char *filename, float *vector, size_t size){
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Erro ao abrir o arquivo para escrita");
        exit(EXIT_FAILURE);
    }
    fwrite(vector, sizeof(float), size, file);
    fclose(file);
}

void copy_cvector_to_real_fvector(MKL_Complex8 *cvector, float *fvector, size_t size){
    for (int i = 0; i < size; i++)
    {
        fvector[i] = cvector[i].real;
    }
}

void copy_fvector_to_cvector(MKL_Complex8 *cvector, float *fvector, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        cvector[i].real = fvector[i];
    }
}
