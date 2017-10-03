//#define MATRICES_1D_V_0_1
#ifdef MATRICES_1D_V_0_1

/*           INITIAL VERSION, v0.1
 * Functions allocate arrays that they return. */

/* Matrices are represented as 1-D arrays in memory.
 * That means they are contiguous in memory, flat arrays.
 * Minimum dimension is 1, not 0, and internal dimensions must match. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Initializes vector or matrix with sequentially growing values. */
void init_seq(double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = i*n_cols_a + j;
        }
    }
}

/* Initializes vector or matrix, randomly. */
void init_rand(double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = rand() / (double)RAND_MAX;
        }
    }
}

/* Mean value of an array */
double mean(const double *arr, const unsigned length) {
    double sum = 0.;
    for (size_t i = 0; i < length; i++) {
        sum += arr[i];
    }
    return sum / length;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array. */
double *dot(const double *a, const unsigned n_rows_a, const unsigned n_cols_a,\
            const double *b, const unsigned n_rows_b, const unsigned n_cols_b) {

    /* Check lengths of the input arrays */
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
double *add_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    unsigned size = n_a > n_b ? n_a : n_b;

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
double *subtract_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    unsigned size = n_a > n_b ? n_a : n_b;

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
double *multiply_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    unsigned size = n_a > n_b ? n_a : n_b;

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
double *divide_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    unsigned size = n_a > n_b ? n_a : n_b;

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

/*  Updates an array, element-wise, by adding another array to it.
    Takes both arrays in, and returns the updated one (the first one).
    The return value (address of the first array) doesn't have to be used.
    Arrays must be of the same length, or, the second one can be a scalar.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *add_update(double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if (n_a == 0) {
        printf("'A' cannot be a scalar!\n");
        system("pause");
        exit(-2);
    }
    if ((n_a != n_b) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    /* b is scalar */
    if (n_b == 0) {
        for (size_t i = 0; i < n_a; i++) {
            a[i] += *b;
        }
    }
    /* b is array */
    else {
        for (size_t i = 0; i < n_a; i++) {
            a[i] += b[i];
        }
    }

    return a;
}

/*  Compares two arrays element-wise.
    It creates a new array with results, that it returns.
    If an element of array a is greater than a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *greater_than(const double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    unsigned size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result = malloc(1 * sizeof(*result));
        if (result != NULL) {
            result[0] = *a > *b;
        }
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        result = malloc(n_b * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_b; i++) {
                result[i] = *a > b[i];
            }
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] > *b;
            }
        }
    }
    /* Neither a nor b are scalars. */
    else {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] > b[i];
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

/*  Compares two arrays element-wise.
    It creates a new array with results, that it returns.
    If an element of array a is equal to a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *equal(const double *a, const unsigned n_a, const double *b, const unsigned n_b) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    double *result = NULL;

    /* Size of result is maximum of n_a and n_b. */
    unsigned size = n_a > n_b ? n_a : n_b;

    /* Both a and b are scalars. */
    if (size == 0) {
        result = malloc(1 * sizeof(*result));
        if (result != NULL) {
            result[0] = *a == *b;
        }
    }
    /* Only a is scalar. */
    else if (n_a == 0) {
        result = malloc(n_b * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_b; i++) {
                result[i] = *a == b[i];
            }
        }
    }
    /* Only b is scalar. */
    else if (n_b == 0) {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] == *b;
            }
        }
    }
    /* Neither a nor b are scalars. */
    else {
        result = malloc(n_a * sizeof(*result));
        if (result != NULL) {
            for (size_t i = 0; i < n_a; i++) {
                result[i] = a[i] == b[i];
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
double *transpose(const double *m, const unsigned n_rows_m, const unsigned n_cols_m) {
    double *t = malloc(n_rows_m * n_cols_m * sizeof(*t));
    if (t == NULL) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

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
void print(const double *m, const unsigned n_rows_m, const unsigned n_cols_m) {
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

    const unsigned n_rows_a = 4;
    const unsigned n_cols_a = 3;
    const unsigned n_rows_b = 3;
    const unsigned n_cols_b = 2;

    const unsigned n_rows_x = 4;
    const unsigned n_cols_x = 3;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_b * n_cols_b * sizeof(*b));
    double *c = NULL;
    double *d = NULL;
    double *e = NULL;
    double *f = NULL;
    double *g = malloc(n_rows_a * n_cols_a * sizeof(*g));
    double *h = NULL;
    double *i = NULL;
    double *x = malloc(n_rows_x * n_cols_x * sizeof(*x));
    double *y = NULL;

    if (!a || !b ) {
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
    e = divide_arrays(&one, 0, &ten, 0);                                    // shape (), that is, (1, 1)
    free(e);
    e = add_arrays(a, n_rows_a * n_cols_a, a, n_rows_a * n_cols_a);         // shape (4, 3)
    f = multiply_arrays(&two, 0, a, n_rows_a * n_cols_a);                   // shape (4, 3)
    memcpy(g, a, n_rows_a * n_cols_a * sizeof(*a));
    add_update(g, n_rows_a * n_cols_a, &ten, 0);                            // shape (4, 3)
    g = add_update(g, n_rows_a * n_cols_a, e, n_rows_a * n_cols_a);         // shape (4, 3)
    h = greater_than(&two, 0, &one, 0);                                     // shape (), that is, (1, 1)
    free(h);
    h = greater_than(&one, 0, f, n_rows_a * n_cols_a);                      // shape (4, 3)
    free(h);
    h = greater_than(f, n_rows_a * n_cols_a, &one, 0);                      // shape (4, 3)
    free(h);
    init_rand(x, n_rows_x, n_cols_x);                                       // shape (4, 3)
    h = greater_than(f, n_rows_a * n_cols_a, x, n_rows_x * n_cols_x);       // shape (4, 3)
    i = equal(e, n_rows_a * n_cols_a, f, n_rows_a * n_cols_a);              // shape (4, 3)

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
    printf("Matrix G:\n");
    print(g, n_rows_a, n_cols_a);
    printf("Matrix X:\n");
    print(x, n_rows_x, n_cols_x);
    printf("Matrix H:\n");
    print(h, n_rows_a, n_cols_a);
    printf("Matrix I:\n");
    print(i, n_rows_a, n_cols_a);

    init_seq(x, n_rows_x, n_cols_x);                                        // shape (4, 3)
    y = transpose(x, n_rows_x, n_cols_x);                                   // shape (3, 4)
    free(i);
    i = equal(x, n_rows_x * n_cols_x, &two, 0);                             // shape (4, 3)

    printf("\n");
    printf("Matrix X:\n");
    print(x, n_rows_x, n_cols_x);
    printf("Matrix Y:\n");
    print(y, n_cols_x, n_rows_x);
    printf("Matrix I:\n");
    print(i, n_rows_a, n_cols_a);

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
    free(g);
    free(h);
    free(i);
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

#endif // MATRICES_1D_V_0_1
