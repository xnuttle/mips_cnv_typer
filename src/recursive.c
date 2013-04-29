#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>

#define PARALOGS 2

/*
 * Given n values and an array x of integers for each n value, populate a matrix
 * with n columns and rows equal to the product of the length of the n-1 matrix
 * and x[n].
 */
gsl_matrix_int* populate_matrix(int n, int* x) {
    int i, j, k;
    gsl_matrix_int* matrix;

    printf("n=%i (CN: %i)\n", n, x[n-1]);

    if (n == 1) {
        printf("Allocate matrix (%i x %i)\n", x[n-1], n);
        matrix = gsl_matrix_int_alloc(x[n-1], n);

        for (i=0; i < x[n-1]; i++) {
            gsl_matrix_int_set(matrix, i, n-1, i);
        }
    }
    else {
        gsl_matrix_int* previous_matrix = populate_matrix(n-1, x);
        int previous_rows = (int)previous_matrix->size1;
        printf("Previous matrix had %i rows\n", previous_rows);

        printf("Allocate matrix (%i x %i)\n", previous_rows*x[n-1], n);
        matrix = gsl_matrix_int_alloc(previous_rows*x[n-1], n);

        /*
         * First iterator i is the number of copies for the current paralog.
         * Second iterator j is the length of the previous matrix.
         */
        for (i = 0; i < x[n-1]; i++) {
            for (j = 0; j < previous_rows; j++) {
                for(k = 0; k < n - 1; k++) {
                    printf("Set (%i, %i) from previous\n", j, k);
                    gsl_matrix_int_set(
                        matrix,
                        previous_rows * i + j,
                        k,
                        gsl_matrix_int_get(previous_matrix, j, k)
                    );
                }
                gsl_matrix_int_set(matrix, previous_rows * i + j, n-1, i);
            }
        }

        printf("Freeing previous matrix (%i)\n", n - 1);
        gsl_matrix_int_free(previous_matrix);
    }

    printf("\n");
    return matrix;
}

int main(int argc, char*argv[]) {
    gsl_matrix_int* matrix = NULL;
    int x[PARALOGS] = {2, 3};
    int i, j;
    matrix = populate_matrix(PARALOGS, x);

    for (i = 0; i < (int)matrix->size1; i++) {
        for (j = 0; j < (int)matrix->size2; j++) {
            printf("%i ", gsl_matrix_int_get(matrix, i, j));
        }
        printf("\n");
    }

    if (matrix != NULL) {
        printf("Freeing final matrix\n");
        gsl_matrix_int_free(matrix);
    }

    return 0;
}
