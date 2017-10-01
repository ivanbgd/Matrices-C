#define MATRICES_1D_OPEN_MP
#ifdef MATRICES_1D_OPEN_MP

/*                  FINAL VERSION
 * Functions don't allocate arrays that they return. */

/* Matrices are represented as 1-D arrays in memory.
 * That means they are contiguous in memory, flat arrays.
 * Minimum dimension is 1, not 0, and internal dimensions must match. */

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

 /* Initializes vector or matrix, sequentially, with indices. */
void init_seq(double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    int i = 0;

#pragma omp parallel for default(none) private(i) shared(a, n_rows_a, n_cols_a) schedule(static)
    for (i = 0; i < (int)n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = i*n_cols_a + j;
        }
    }
}

/* Initializes vector or matrix, randomly. */
void init_rand(double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    int i = 0;

    /* Schedule has to be either guided or dynamic; if it's static or runtime, the random numbers repeat. */
#pragma omp parallel for default(none) private(i) shared(a, n_rows_a, n_cols_a) schedule(dynamic)
    for (i = 0; i < (int)n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = rand() / (double)RAND_MAX;
        }
    }
}

/* Mean value of an array */
double mean(const double *arr, const unsigned length) {
    double sum = 0.;
    int i = 0;

#pragma omp parallel for default(none) private (i) shared(arr, length) reduction(+:sum) schedule(static)
    for (i = 0; i < (int)length; i++) {
        sum += arr[i];
    }
    return sum / length;
}

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
    It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
    of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
    The original matrix m stays intact. */
double *transpose(const double *m, const unsigned n_rows_m, const unsigned n_cols_m, double *t) {
    int i = 0, j = 0;
    
#pragma omp parallel for default(none) private(i, j) shared(m, n_rows_m, n_cols_m, t) schedule(static)
    for (i = 0; i < (int)n_rows_m; i++) {
        for (j = 0; j < (int)n_cols_m; j++) {
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

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is a slower version of the function. */
double *dot_simple(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
    const double *b, const unsigned n_rows_b, const unsigned n_cols_b, double *c) {

    int i = 0, j = 0, k = 0;

    /* Check lengths of the input arrays */
    if (n_cols_a != n_rows_b) {
        printf("#columns A must be equal to #rows B!\n");
        system("pause");
        exit(-2);
    }

    double *bt = malloc(n_rows_b * n_cols_b * sizeof(*b));

    bt = transpose(b, n_rows_b, n_cols_b, bt);

#pragma omp parallel for default(none) private(i, j, k) shared(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c, bt) schedule(static)
    for (i = 0; i < (int)n_rows_a; i++) {
        for (k = 0; k < (int)n_cols_b; k++) {
            double sum = 0.0;
            for (j = 0; j < (int)n_cols_a; j++) {
                sum += a[i*n_cols_a + j] * bt[k*n_rows_b + j];
            }
            c[i*n_cols_b + k] = sum;
        }
    }

    free(bt);

    return c;
}

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is a faster version of the function. */
double *dot(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
    const double *b, const unsigned n_rows_b, const unsigned n_cols_b, double *c) {

    int i = 0, j = 0, k = 0;

    /* Check lengths of the input arrays */
    if (n_cols_a != n_rows_b) {
        printf("#columns A must be equal to #rows B!\n");
        system("pause");
        exit(-2);
    }

    double *bt = malloc(n_rows_b * n_cols_b * sizeof(*b));

    bt = transpose(b, n_rows_b, n_cols_b, bt);

#pragma omp parallel for default(none) private(i, j, k) shared(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c, bt) schedule(static)
    for (i = 0; i < (int)n_rows_a; i++) {
        for (k = 0; k < (int)n_cols_b; k++) {
            double sum = 0.0;
            for (j = 0; j < (int)n_cols_a; j++) {
                sum += a[i*n_cols_a + j] * bt[k*n_rows_b + j];
            }
            c[i*n_cols_b + k] = sum;
        }
    }

    free(bt);

    return c;
}

/*  Adds two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *add_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    int i = 0;

    /* Neither a nor b are scalars. */
    if ((n_a > 0) && (n_b > 0)) {        
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] + b[i];
        }
    }
    /* Only b is scalar. */
    else if ((n_b == 0) && (n_a > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] + *b;
        }
    }
    /* Only a is scalar. */
    else if ((n_a == 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_b; i++) {
            result[i] = *a + b[i];
        }
    }
    /* Both a and b are scalars. */
    else {
        result[0] = *a + *b;
    }

    return result;
}

/*  Subtracts the second array from the first one, element-wise, and puts the result
    in an array that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *subtract_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }
    
    int i = 0;

    /* Neither a nor b are scalars. */
    if ((n_a > 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] - b[i];
        }
    }
    /* Only b is scalar. */
    else if ((n_b == 0) && (n_a > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] - *b;
        }
    }
    /* Only a is scalar. */
    else if ((n_a == 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_b; i++) {
            result[i] = *a - b[i];
        }
    }
    /* Both a and b are scalars. */
    else {
        result[0] = *a - *b;
    }

    return result;
}

/*  Multiplies two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *multiply_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    int i = 0;

    /* Neither a nor b are scalars. */
    if ((n_a > 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] * b[i];
        }
    }
    /* Only b is scalar. */
    else if ((n_b == 0) && (n_a > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] * *b;
        }
    }
    /* Only a is scalar. */
    else if ((n_a == 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_b; i++) {
            result[i] = *a * b[i];
        }
    }
    /* Both a and b are scalars. */
    else {
        result[0] = *a * *b;
    }

    return result;
}

/*  Multiplies two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *divide_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    int i = 0;

    /* Neither a nor b are scalars. */
    if ((n_a > 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] / b[i];
        }
    }
    /* Only b is scalar. */
    else if ((n_b == 0) && (n_a > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] / *b;
        }
    }
    /* Only a is scalar. */
    else if ((n_a == 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_b; i++) {
            result[i] = *a / b[i];
        }
    }
    /* Both a and b are scalars. */
    else {
        result[0] = *a / *b;
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

    int i = 0;

    /* b is scalar */
    if (n_b == 0) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            a[i] += *b;
        }
    }
    /* b is array */
    else {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            a[i] += b[i];
        }
    }

    return a;
}

/*  Compares two arrays element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    If an element of array a is greater than a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *greater_than(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }

    int i = 0;

    /* Neither a nor b are scalars. */
    if ((n_a > 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] > b[i];
        }
    }
    /* Only b is scalar. */
    else if ((n_b == 0) && (n_a > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] > *b;
        }
    }
    /* Only a is scalar. */
    else if ((n_a == 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_b; i++) {
            result[i] = *a > b[i];
        }
    }
    /* Both a and b are scalars. */
    else {
        result[0] = *a > *b;
    }

    return result;
}

/*  Compares two arrays element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    If an element of array a is equal to a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *equal(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result) {
    /* Check lengths of the input arrays */
    if ((n_a != n_b) && (n_a != 0) && (n_b != 0)) {
        printf("Length of A must be equal to length of B!\n");
        system("pause");
        exit(-2);
    }
    
    int i = 0;
    
    /* Neither a nor b are scalars. */
    if ((n_a > 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] == b[i];
        }
    }
    /* Only b is scalar. */
    else if ((n_b == 0) && (n_a > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_a; i++) {
            result[i] = a[i] == *b;
        }
    }
    /* Only a is scalar. */
    else if ((n_a == 0) && (n_b > 0)) {
    #pragma omp parallel for default(none) private(i) shared(a, n_a, b, n_b, result) schedule(static)
        for (i = 0; i < (int)n_b; i++) {
            result[i] = *a == b[i];
        }
    }
    /* Both a and b are scalars. */
    else {
        result[0] = *a == *b;
    }

    return result;
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
    double *c = malloc(n_rows_a * n_cols_b * sizeof(*c));
    double *d = malloc(n_rows_a * n_cols_b * sizeof(*d));
    double *e = malloc(n_rows_a * n_cols_a * sizeof(*e));
    double *f = malloc(n_rows_a * n_cols_a * sizeof(*f));
    double *g = malloc(n_rows_a * n_cols_a * sizeof(*g));
    double *h = malloc(n_rows_a * n_cols_a * sizeof(*h));
    double *i = malloc(n_rows_a * n_cols_a * sizeof(*i));
    double *x = malloc(n_rows_x * n_cols_x * sizeof(*x));
    double *y = malloc(n_rows_x * n_cols_x * sizeof(*y));

    if (!a || !b || !c || !d || !e || !f || !g || !x || !y) {
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
    memcpy(g, a, n_rows_a * n_cols_a * sizeof(*a));
    add_update(g, n_rows_a * n_cols_a, &ten, 0);                            // shape (4, 3)
    g = add_update(g, n_rows_a * n_cols_a, e, n_rows_a * n_cols_a);         // shape (4, 3)
    h = greater_than(&two, 0, &one, 0, h);                                  // shape (), that is, (1, 1)
    h = greater_than(&one, 0, f, n_rows_a * n_cols_a, h);                   // shape (4, 3)
    h = greater_than(f, n_rows_a * n_cols_a, &one, 0, h);                   // shape (4, 3)
    init_rand(x, n_rows_x, n_cols_x);                                       // shape (4, 3)
    h = greater_than(f, n_rows_a * n_cols_a, x, n_rows_x * n_cols_x, h);    // shape (4, 3)
    i = equal(e, n_rows_a * n_cols_a, f, n_rows_a * n_cols_a, i);           // shape (4, 3)

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
    transpose(x, n_rows_x, n_cols_x, y);                                    // shape (3, 4)
    i = equal(x, n_rows_x * n_cols_x, &two, 0, i);                          // shape (4, 3)

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

void test_speed_1() {
    /* For measuring time */
    double t0, t1;

    const unsigned scale = 2000;
    const unsigned n_rows_a = 4 * scale;
    const unsigned n_cols_a = 3 * scale;
    const unsigned n_rows_b = 4 * scale;
    const unsigned n_cols_b = 3 * scale;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_b * n_cols_b * sizeof(*b));
    double *c = malloc(n_rows_a * n_cols_b * sizeof(*c));

    if (!a || !b || !c) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    puts("");
    omp_set_num_threads(4);
    printf("omp_get_num_procs %i\n", omp_get_num_procs());
    printf("omp_get_max_threads %i\n", omp_get_max_threads());
    puts("");

    init_rand(a, n_rows_a, n_cols_a);
    init_rand(b, n_rows_b, n_cols_b);
    init_rand(c, n_rows_a, n_cols_a);

    const unsigned loop = 200u;

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        c = multiply_arrays(a, n_rows_a * n_cols_a, c, n_rows_a * n_cols_a, c);
    }
    t1 = omp_get_wtime();
    printf("multiply_arrays(): Elapsed time %.3f s\n", t1 - t0);

    if (scale == 1) {
        puts("");
        printf("Matrix A:\n");
        print(a, n_rows_a, n_cols_a);
        printf("Matrix B:\n");
        print(b, n_rows_b, n_cols_b);
        printf("Matrix C:\n");
        print(c, n_rows_a, n_cols_a);
    }

    a = multiply_arrays(c, n_rows_a * n_cols_a, b, n_rows_b * n_cols_b, a);

    free(a);
    free(b);
    free(c);
}

void test_speed_2() {
    /* For measuring time */
    double t0, t1;

    const unsigned scale = 200;
    const unsigned n_rows_a = 4 * scale;
    const unsigned n_cols_a = 3 * scale;
    const unsigned n_rows_b = 3 * scale;
    const unsigned n_cols_b = 2 * scale;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_b * n_cols_b * sizeof(*b));
    double *c = malloc(n_rows_a * n_cols_b * sizeof(*c));
    double *d = malloc(n_rows_a * n_cols_b * sizeof(*d));

    if (!a || !b || !c || !d) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    puts("");

    init_rand(a, n_rows_a, n_cols_a);
    init_rand(b, n_rows_b, n_cols_b);

    const unsigned loop = 100u;

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        c = dot_simple(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c);
    }
    t1 = omp_get_wtime();
    printf("dot_simple(): Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        d = dot(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, d);
    }
    t1 = omp_get_wtime();
    printf("dot(): Elapsed time %.3f s\n", t1 - t0);

    if (scale == 1) {
        puts("");
        printf("Matrix A:\n");
        print(a, n_rows_a, n_cols_a);
        printf("Matrix B:\n");
        print(b, n_rows_b, n_cols_b);
        printf("Matrix C:\n");
        print(c, n_rows_a, n_cols_b);
        printf("Matrix D:\n");
        print(d, n_rows_a, n_cols_b);
    }

    free(a);
    free(b);
    free(c);
    free(d);
}

int main(int argc, char *argv[]) {
    /* Intializes random number generator */
    time_t t;
    srand((unsigned)time(&t));
    srand(0);

    test();
    test_speed_2();

    system("pause");
    return(0);
}

#endif // MATRICES_1D_OPEN_MP
