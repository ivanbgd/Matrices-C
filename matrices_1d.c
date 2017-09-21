/*                  FINAL VERSION
 * Functions don't allocate arrays that they return. */

/* Matrices are represented as 1-D arrays in memory.
 * That means they are contiguous in memory, flat arrays.
 * Minimum dimension is 1, not 0, and internal dimensions must match. */

#define MATRICES_1D
#ifdef MATRICES_1D

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

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c. */
double *dot(const double *a, const int n_rows_a, const int n_cols_a, const double *b, const int n_rows_b, const int n_cols_b, double *c) {
    /* Check lengths of the input arrays */
    if (n_cols_a != n_rows_b) {
        printf("#columns A must be equal to #rows B!\n");
        system("pause");
        exit(-2);
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

/*  Adds two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *add_arrays(const double *a, const int n_a, const double *b, const int n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result[0] = *a + *b;
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        for (size_t i = 0; i < n_b; i++) {
            result[i] = *a + b[i];
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] + *b;
        }
    }
    /* Neither a nor b are scalars. */
    else {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] + b[i];
        }
    }

    return result;
}

/*  Subtracts the second array from the first one, element-wise, and puts the result
    in an array that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *subtract_arrays(const double *a, const int n_a, const double *b, const int n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result[0] = *a - *b;

    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        for (size_t i = 0; i < n_b; i++) {
            result[i] = *a - b[i];
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] - *b;
        }
    }
    /* Neither a nor b are scalars. */
    else {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] - b[i];
        }
    }
    
    return result;
}

/*  Multiplies two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *multiply_arrays(const double *a, const int n_a, const double *b, const int n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result[0] = *a * *b;
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        for (size_t i = 0; i < n_b; i++) {
            result[i] = *a * b[i];
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] * *b;
        }
    }
    /* Neither a nor b are scalars. */
    else {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] * b[i];
        }
    }

    return result;
}

/*  Multiplies two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *divide_arrays(const double *a, const int n_a, const double *b, const int n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    /* Size of result is maximum of n_a and n_b. */
    int size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result[0] = *a / *b;
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        for (size_t i = 0; i < n_b; i++) {
            result[i] = *a / b[i];
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] / *b;
        }
    }
    /* Neither a nor b are scalars. */
    else {
        for (size_t i = 0; i < n_a; i++) {
            result[i] = a[i] / b[i];
        }
    }

    return result;
}

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
    It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
    of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
    The original matrix m stays intact. */
double *transpose(const double *m, const int n_rows_m, const int n_cols_m, double *t) {
    for (size_t i = 0; i < n_rows_m; i++) {
        for (size_t j = 0; j < n_cols_m; j++) {
            t[j*n_rows_m + i] = m[i*n_cols_m + j];
        }
    }

    /* Visual validation - Prints t like m, the original */
    const int validate = 0;
    if (validate) {
        for (size_t i = 0; i < n_rows_m; i++) {
            for (size_t j = 0; j < n_cols_m; j++) {
                printf("%8.3f ", t[j*n_rows_m + i]);
            }
            printf("\n");
        }
        printf("\n");
    }

    return t;
}

/* Prints vector, or matrix. */
void print(const double *m, const int n_rows_m, const int n_cols_m) {
    for (size_t i = 0; i < n_rows_m; i++) {
        for (size_t j = 0; j < n_cols_m; j++) {
            printf("%8.3f ", m[i*n_cols_m + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void test() {
    const double one = 1.0;
    const double two = 2.0;
    const double ten = 10.0;

    const int n_rows_a = 4;
    const int n_cols_a = 3;
    const int n_rows_b = 3;
    const int n_cols_b = 2;

    const int n_rows_x = 4;
    const int n_cols_x = 3;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_b * n_cols_b * sizeof(*b));
    double *c = malloc(n_rows_a * n_cols_b * sizeof(*c));
    double *d = malloc(n_rows_a * n_cols_b * sizeof(*d));
    double *e = malloc(n_rows_a * n_cols_a * sizeof(*e));
    double *f = malloc(n_rows_a * n_cols_a * sizeof(*f));
    double *x = malloc(n_rows_x * n_cols_x * sizeof(*x));
    double *y = malloc(n_rows_x * n_cols_x * sizeof(*y));

    if (!a || !b || !c || !d || !e || !f || !x || !y) {
        printf("Couldn't allocate memory!\a\n");
        system("pause");
        exit(-1);
    }

    init_rand(a, n_rows_a, n_cols_a);
    init_rand(b, n_rows_b, n_cols_b);

    init_seq(a, n_rows_a, n_cols_a);
    init_seq(b, n_rows_b, n_cols_b);

    c = dot(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c);               // shape (4, 2)
    d = transpose(c, n_rows_a, n_cols_b, d);                                // shape (2, 4)
    e = divide_arrays(&ten, 0, a, n_rows_a * n_cols_a, e);                  // shape (4, 3)
    e = divide_arrays(a, n_rows_a * n_cols_a, &ten, 0, e);                  // shape (4, 3)
    memset(e, -1, n_rows_a * n_cols_a * sizeof(*e));
    e = divide_arrays(&one, 0, &ten, 0, e);                                 // shape (), that is, (1, 1)
    e = add_arrays(a, n_rows_a * n_cols_a, a, n_rows_a * n_cols_a, e);      // shape (4, 3)
    f = multiply_arrays(&two, 0, a, n_rows_a * n_cols_a, f);                // shape (4, 3)

    init_seq(x, n_rows_x, n_cols_x);                                        // shape (4, 3)
    transpose(x, n_rows_x, n_cols_x, y);                                    // shape (3, 4)

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

    printf("Matrix X:\n");
    print(x, n_rows_x, n_cols_x);
    printf("Matrix Y:\n");
    print(y, n_cols_x, n_rows_x);

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
    free(x);
    free(y);
}

int main(int argc, char *argv[]) {
    /* Intializes random number generator */
    time_t t;
    srand((unsigned)time(&t));

    test();

    system("pause");
    return(0);
}

#endif // MATRICES_1D
