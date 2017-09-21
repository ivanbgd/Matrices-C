/* Matrices are represented as 1-D arrays in memory.
* That means they are contiguous in memory, flat arrays.
* Minimum dimension is 1, not 0, and internal dimensions must match. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Initializes vector or matrix, sequentially, with indices. */
void init_seq(double *a, const int n_rows_a, const int n_cols_a) {
    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = i*n_cols_a + j;
        }
    }
}

/* Initializes vector or matrix, randomly. */
void init_rand(double *a, const int n_rows_a, const int n_cols_a) {
    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = rand() / (double)RAND_MAX;
        }
    }
}

/* Mean value of an array */
double mean(const double *arr, const int size) {
    double sum = 0.;
    for (size_t i = 0; i < size; i++) {
        sum += arr[i];
    }
    return sum / size;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array. */
double *dot(const double *a, const int n_rows_a, const int n_cols_a, const double *b, const int n_rows_b, const int n_cols_b) {
    if (n_cols_a != n_rows_b) {
        printf("#columns A must be equal to #rows B!\n");
        system("pause");
        exit(-2);
    }

    double *c = malloc(n_rows_a * n_cols_b * sizeof(*c));
    if (c == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t k = 0; k < n_cols_b; k++) {
            double sum = 0.0;
            for (size_t j = 0; j < n_cols_a; j++) {
                sum += a[i*n_cols_a + j] * b[j*n_cols_b + k];
            }
            c[i*n_cols_b + k] = sum;
        }
    }

    return c;
}

/*  Adds two arrays, element-wise, and puts
    the result in an array that it allocates and returns.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *add_arrays(const double *a, const int n_a, const double *b, const int n_b) {
    /* Check lengths of both arrays, just in case. */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result = malloc(1 * sizeof(*result));
        if (result != NULL) {
            result[0] = *a + *b;
        }
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        result = malloc(n_b * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_b; i++) {
                result[i] = *a + b[i];
            }
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] + *b;
            }
        }
    }
    /* Neither a nor b are scalars. */
    else {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] + b[i];
            }
        }
    }

    if (result == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    return result;
}

/*  Subtracts the second array from the first one, element-wise, and puts
    the result in an array that it allocates and returns.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *subtract_arrays(const double *a, const int n_a, const double *b, const int n_b) {
    /* Check lengths of both arrays, just in case. */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result = malloc(1 * sizeof(*result));
        if (result != NULL) {
            result[0] = *a - *b;
        }
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        result = malloc(n_b * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_b; i++) {
                result[i] = *a - b[i];
            }
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] - *b;
            }
        }
    }
    /* Neither a nor b are scalars. */
    else {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] - b[i];
            }
        }
    }

    if (result == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    return result;
}

/*  Multiplies two arrays element-wise, and puts
    the result in an array that it allocates and returns.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *multiply_arrays(const double *a, const int n_a, const double *b, const int n_b) {
    /* Check lengths of both arrays, just in case. */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result = malloc(1 * sizeof(*result));
        if (result != NULL) {
            result[0] = *a * *b;
        }
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        result = malloc(n_b * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_b; i++) {
                result[i] = *a * b[i];
            }
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] * *b;
            }
        }
    }
    /* Neither a nor b are scalars. */
    else {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] * b[i];
            }
        }
    }

    if (result == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    return result;
}

/*  Divides two arrays element-wise, and puts
    the result in an array that it allocates and returns.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *divide_arrays(const double *a, const int n_a, const double *b, const int n_b) {
    /* Check lengths of both arrays, just in case. */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result = malloc(1 * sizeof(*result));
        if (result != NULL) {
            result[0] = *a / *b;
        }
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        result = malloc(n_b * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_b; i++) {
                result[i] = *a / b[i];
            }
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] / *b;
            }
        }
    }
    /* Neither a nor b are scalars. */
    else {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] / b[i];
            }
        }
    }

    if (result == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    return result;
}

/*  Allocates and returns a new matrix, which is a transpose of the input one.
    It's still flat in memory, i.e., 1-D. */
double *transpose(const double *m, const int n_rows_m, const int n_cols_m) {
    double *t = malloc(n_cols_m * n_rows_m * sizeof(*t));
    if (t == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    for (size_t i = 0; i < n_cols_m; i++) {
        for (size_t j = 0; j < n_rows_m; j++) {
            t[i*n_rows_m + j] = m[j*n_cols_m + i];
        }
    }

    /* Visual validation */
    const int validate = 0;
    if (validate) {
        for (size_t i = 0; i < n_cols_m; i++) {
            for (size_t j = 0; j < n_rows_m; j++) {
                printf("%8.3f ", t[i*n_rows_m + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    return t;
}

/* Prints vector, or matrix. */
void print(const double *a, const int n_rows_a, const int n_cols_a) {
    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            printf("%8.3f ", a[i*n_cols_a + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    /* Intializes random number generator */
    time_t t;
    srand((unsigned)time(&t));

    const double one =  1.0;
    const double two =  2.0;
    const double ten = 10.0;

    const int n_rows_a = 4;
    const int n_cols_a = 3;
    const int n_rows_b = 3;
    const int n_cols_b = 2;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_b * n_cols_b * sizeof(*b));
    double *c = NULL;
    double *d = NULL;
    double *e = NULL;
    double *f = NULL;

    if (!a || !b) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    init_rand(a, n_rows_a, n_cols_a);
    init_rand(b, n_rows_b, n_cols_b);

    c = dot(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b);                  // shape (4, 2)
    d = transpose(c, n_rows_a, n_cols_b);                                   // shape (2, 4)
    e = divide_arrays(&ten, 0, a, n_rows_a * n_cols_a);                     // shape (4, 3)
    free(e);
    e = divide_arrays(a, n_rows_a * n_cols_a, &ten, 0);                     // shape (4, 3)
    free(e);
    e = divide_arrays(&one, 0, &ten, 0);                                    // shape (4, 3)
    free(e);
    e = add_arrays(a, n_rows_a * n_cols_a, a, n_rows_a * n_cols_a);         // shape (4, 3)
    f = multiply_arrays(&two, 0, a, n_rows_a * n_cols_a);                   // shape (4, 3)

    printf("Matrix A:\n");
    print(a, n_rows_a, n_cols_a);
    printf("Matrix B:\n");
    print(b, n_rows_b, n_cols_b);
    printf("Matrix C:\n");
    print(c, n_rows_a, n_cols_b);
    printf("Matrix D:\n");
    print(d, n_cols_b, n_rows_a);
    printf("Matrix E:\n");
    print(e, n_rows_a, n_cols_a);
    printf("Matrix F:\n");
    print(f, n_rows_a, n_cols_a);

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);

    system("pause");
    return(0);
}

