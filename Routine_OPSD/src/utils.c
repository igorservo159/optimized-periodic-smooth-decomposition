#include "../include/utils.h"

int check_args(const char *BIN, const char *ROUTINE, const char *PRECISION, const char *SAVE_VECTORS, const char *INPUT)
{
    if(strcmp(ROUTINE, "ccr") && strcmp(ROUTINE, "css") && strcmp(ROUTINE, "cts")){
        printf("Use: %s <rows> <columns> <routine> <precision> <save_vectors> <input> <directory> <seed>\n", BIN);
        printf("Options to <routine>: 'ccr', 'cts', 'css'\n");
        return -2;

        //compute completed routine - ccr
        //compute shifted spectrums - css
        //compute tradicional spectrums - cts
    }

    if(strcmp(PRECISION, "single") && strcmp(PRECISION, "double")){
        printf("Use: %s <rows> <columns> <routine> <precision> <save_vectors> <input> <directory> <seed>\n", BIN);
        printf("Options to <precision>: 'single', 'double'\n");
        return -3;
    }

    if(strcmp(SAVE_VECTORS, "yes") && strcmp(SAVE_VECTORS, "no")){
        printf("Use: %s <rows> <columns> <routine> <precision> <save_vectors> <input> <directory> <seed>\n", BIN);
        printf("Options to <save_vectors>: 'yes', 'no'\n");
        return -4;
    }

    if(strcmp(INPUT, "rb") && strcmp(INPUT, "fm")){
        printf("Use: %s <rows> <columns> <routine> <precision> <save_vectors> <input> <directory> <seed>\n", BIN);
        printf("Options to <input>: 'rb', 'fm'\n");
        return -5;

        //read binary - rb
        //fill matrix - fm
    }

    return 0;
}

void ensure_directory_exists(const char *path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0700) != 0) {
            perror("mkdir");
            exit(EXIT_FAILURE);
        }
    }
}

//CVECTOR FUNCTIONS

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

void free_cvector(MKL_Complex8 *vector)
{
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

void copy_cvector_to_real_fvector(MKL_Complex8 *cvector, float *fvector, size_t size){
    for (int i = 0; i < size; i++)
    {
        fvector[i] = cvector[i].real;
    }
}

void copy_cvector_to_real_dvector(MKL_Complex8 *cvector, double *dvector, size_t size){
    for (int i = 0; i < size; i++)
    {
        dvector[i] = cvector[i].real;
    }
}

//ZVECTOR FUNCTIONS

void init_zvector(MKL_Complex16 **vector, size_t size)
{
    (*vector) = (MKL_Complex16 *)malloc(size * sizeof(MKL_Complex16));
    if ((*vector) == NULL)
    {
        printf("Allocation error!\n");
        free(*vector);
        return;
    }
}

void free_zvector(MKL_Complex16 *vector)
{
    if (vector == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    free(vector);
}

void show_zmatrix(MKL_Complex16 *matrix, size_t rows, size_t columns)
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

void fill_zmatrix(MKL_Complex16 *matrix, size_t rows, size_t columns, unsigned int seed)
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
            double random_real = (double)rand() / RAND_MAX_D;
            double scaled_real = random_real * 4 + 1;

            // double random_complex = (double)rand() / RAND_MAX_D;
            // double scaled_complex = random_real * 3 + 1;

            // matrix[i + j * rows].real = 1.0;
            matrix[i + j * rows].real = scaled_real;
            matrix[i + j * rows].imag = 0.0;
        }
    }
}

void read_zvector_bin(const char *filename, MKL_Complex16 *vector, size_t size){
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Error opening file");
        exit(1);
    }

    size_t elements_read = fread(vector, sizeof(MKL_Complex16), size, fp);
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

void save_zvector_on_bin(const char *filename, MKL_Complex16 *vector, size_t size){
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Erro ao abrir o arquivo para escrita");
        exit(EXIT_FAILURE);
    }
    fwrite(vector, sizeof(MKL_Complex16), size, file);
    fclose(file);
}

void copy_zvector_to_real_fvector(MKL_Complex16 *zvector, float *fvector, size_t size){
    for (int i = 0; i < size; i++)
    {
        fvector[i] = zvector[i].real;
    }
}

void copy_zvector_to_real_dvector(MKL_Complex16 *zvector, double *dvector, size_t size){
    for (int i = 0; i < size; i++)
    {
        dvector[i] = zvector[i].real;
    }
}

//FVECTOR FUNCTIONS

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

void free_fvector(float *vector){
    if (vector == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    free(vector);
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

void save_fvector_on_bin(const char *filename, float *vector, size_t size){
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Erro ao abrir o arquivo para escrita");
        exit(EXIT_FAILURE);
    }
    fwrite(vector, sizeof(float), size, file);
    fclose(file);
}

void copy_fvector_to_cvector(MKL_Complex8 *cvector, float *fvector, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        cvector[i].real = fvector[i];
    }
}

void copy_fvector_to_zvector(MKL_Complex16 *zvector, float *fvector, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        zvector[i].real = fvector[i];
    }
}

//DVECTOR FUNCTIONS

void init_dvector(double **vector, size_t size)
{
    (*vector) = (double *)malloc(size * sizeof(double));
    if ((*vector) == NULL)
    {
        printf("Allocation error!\n");
        free(*vector);
        return;
    }
}

void free_dvector(double *vector){
    if (vector == NULL)
    {
        printf("Matrix not found!\n");
        return;
    }

    free(vector);
}

void read_dvector_bin(const char *filename, double *vector, size_t size)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Error opening file");
        exit(1);
    }

    size_t elements_read = fread(vector, sizeof(double), size, fp);
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

void save_dvector_on_bin(const char *filename, double *vector, size_t size){
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Erro ao abrir o arquivo para escrita");
        exit(EXIT_FAILURE);
    }
    fwrite(vector, sizeof(double), size, file);
    fclose(file);
}

void copy_dvector_to_cvector(MKL_Complex8 *cvector, double *dvector, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        cvector[i].real = dvector[i];
    }
}

void copy_dvector_to_zvector(MKL_Complex16 *zvector, double *dvector, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        zvector[i].real = dvector[i];
    }
}