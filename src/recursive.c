#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "iniparser.h"

#define PARALOG_COUNTS_SECTION "paralog_counts"

/*
 * Given a configuration file, load the paralog copy number states for all
 * paralogs into a vector.
 */
gsl_vector_int* get_paralog_copy_numbers(char* configuration_filename) {
    gsl_vector_int* paralog_copy_numbers = NULL;
    dictionary* ini;
    int total_paralogs, i;
    char** paralog_keys;

    /*
     * Load the configuration file.
     */
    ini = iniparser_load(configuration_filename);
    if (ini == NULL) {
        fprintf(stderr, "Cannot open configuration file: %s\n", configuration_filename);
        return NULL;
    }

    /*
     * Get total number of paralogs defined in the configuration file and keys
     * for each paralog entry in the file.
     */
    total_paralogs = iniparser_getsecnkeys(ini, PARALOG_COUNTS_SECTION);
    paralog_keys = iniparser_getseckeys(ini, PARALOG_COUNTS_SECTION);
    printf("Found %i paralogs\n", total_paralogs);

    /*
     * Allocate a vector for the paralog copy numbers and populate the vector
     * from the configuration file using the configuration keys.
     */
    paralog_copy_numbers = gsl_vector_int_alloc(total_paralogs);
    for (i = 0; i < total_paralogs; i++) {
        gsl_vector_int_set(paralog_copy_numbers, i, iniparser_getint(ini, paralog_keys[i], -1));
        printf("Paralog %i: %i\n", i, gsl_vector_int_get(paralog_copy_numbers, i));
    }

    iniparser_freedict(ini);

    return paralog_copy_numbers;
}

/*
 * Given n values and an array x of integers for each n value, populate a matrix
 * with n columns and rows equal to the product of the length of the n-1 matrix
 * and x[n].
 */
gsl_matrix_int* populate_matrix(int n, gsl_vector_int* x) {
    int i, j, k;
    gsl_matrix_int* matrix;

    printf("n=%i (CN: %i)\n", n, gsl_vector_int_get(x, n - 1));

    if (n == 1) {
        printf("Allocate matrix (%i x %i)\n", gsl_vector_int_get(x, n - 1), n);
        matrix = gsl_matrix_int_alloc(gsl_vector_int_get(x, n - 1), n);

        for (i=0; i < gsl_vector_int_get(x, n - 1); i++) {
            gsl_matrix_int_set(matrix, i, n-1, i);
        }
    }
    else {
        gsl_matrix_int* previous_matrix = populate_matrix(n-1, x);
        int previous_rows = (int)previous_matrix->size1;
        printf("Previous matrix had %i rows\n", previous_rows);

        printf("Allocate matrix (%i x %i)\n", previous_rows*gsl_vector_int_get(x, n - 1), n);
        matrix = gsl_matrix_int_alloc(previous_rows*gsl_vector_int_get(x, n - 1), n);

        /*
         * Fill in the current matrix with i copies of the previous matrix's
         * values which fill up to n - 2. Populate the final column with values
         * for this paralog.
         */
        for (i = 0; i < gsl_vector_int_get(x, n - 1); i++) {
            for (j = 0; j < previous_rows; j++) {
                for(k = 0; k < n - 1; k++) {
                    gsl_matrix_int_set(
                        matrix,
                        previous_rows * i + j,
                        k,
                        gsl_matrix_int_get(previous_matrix, j, k)
                    );
                }

                // Populate the final column with values for this paralog.
                gsl_matrix_int_set(matrix, previous_rows * i + j, n-1, i);
            }
        }

        printf("Freeing previous matrix (%i)\n", n - 1);
        gsl_matrix_int_free(previous_matrix);
    }

    return matrix;
}

int main(int argc, char*argv[]) {
    gsl_matrix_int* matrix = NULL;
    gsl_vector_int* paralog_copy_numbers = NULL;
    int i, j;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <configuration.ini>\n", argv[0]);
        return -1;
    }

    // Load paralog copy numbers.
    paralog_copy_numbers = get_paralog_copy_numbers(argv[1]);

    // Populate copy states matrix.
    matrix = populate_matrix(paralog_copy_numbers->size, paralog_copy_numbers);

    for (i = 0; i < (int)matrix->size1; i++) {
        for (j = 0; j < (int)matrix->size2; j++) {
            printf("%i ", gsl_matrix_int_get(matrix, i, j));
        }
        printf("\n");
    }

    if (paralog_copy_numbers != NULL) {
        printf("Freeing copy numbers vector\n");
        gsl_vector_int_free(paralog_copy_numbers);
    }

    if (matrix != NULL) {
        printf("Freeing final matrix\n");
        gsl_matrix_int_free(matrix);
    }

    return 0;
}
