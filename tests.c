#include "tests.h"

/* All functions use restricted pointers, so care should be taken
 * to make sure that arrays that they point to do not overlap, if
 * we want to modify them inside of the functions.
 * On the other hand, it's easy to change type of the pointers
 * from restricted to non-restricted versions, by using definitions
 * given at the beginning of the corresponding "matrices_1d.h"
 * header file, if necessary. */

/* Function that calls other test functions */
void test() {
    test_0();
    test_speed_4();
}

/* Simple basic tests */
void test_0() {
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

    const int rand = 0;

    printf("\nRunning \"test_0\" with %s data...\n\n", rand ? "random" : "sequential");

    if (rand) {
        init_rand(a, n_rows_a, n_cols_a);
        init_rand(b, n_rows_b, n_cols_b);
    }
    else {
        init_seq(a, n_rows_a, n_cols_a);
        init_seq(b, n_rows_b, n_cols_b);
    }    
    
    c = dot_faster(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c);        // shape (4, 2)
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

    printf("\n");
    printf("Sum of A: %f\n", sum_array(a, n_rows_a * n_cols_a));
    printf("Mean of A:        %f\n", mean(a, n_rows_a * n_cols_a));
    printf("Mean of A by sum: %f\n", sum_array(a, n_rows_a * n_cols_a) / (n_rows_a * n_cols_a));

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

    puts("");

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

/* Test speed of multiply_arrays */
void test_speed_1() {
    /* For measuring time */
    double t0, t1;

    const unsigned scale = 1500;
    const unsigned n_rows_a = 4 * scale;
    const unsigned n_cols_a = 3 * scale;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_a * n_cols_a * sizeof(*b));
    double *c = malloc(n_rows_a * n_cols_a * sizeof(*c));
    double *d = malloc(n_rows_a * n_cols_a * sizeof(*d));

    const double x = 3.7;

    if (!a || !b || !c || !d) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    const int rand = 1;

    printf("\nRunning \"test_speed_1\" with %s data...\n\n", rand ? "random" : "sequential");

    if (rand) {
        init_rand(a, n_rows_a, n_cols_a);
        init_rand(b, n_rows_a, n_cols_a);
    }
    else {
        init_seq(a, n_rows_a, n_cols_a);
        init_seq(b, n_rows_a, n_cols_a);
    }

    const unsigned loop = 100u;

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        c = multiply_arrays(a, n_rows_a * n_cols_a, b, n_rows_a * n_cols_a, c);
        //c = multiply_arrays(a, n_rows_a * n_cols_a, &x, 0, c);
    }
    t1 = omp_get_wtime();
    printf("multiply_arrays():       Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        d = multiply_arrays_tiled(a, n_rows_a * n_cols_a, b, n_rows_a * n_cols_a, d);
        //d = multiply_arrays_tiled(a, n_rows_a * n_cols_a, &x, 0, d);
    }
    t1 = omp_get_wtime();
    printf("multiply_arrays_tiled(): Elapsed time %.3f s\n", t1 - t0);

    puts("");
    printf("Compare (memcmp) multiply_arrays_tiled() with multiply_arrays(): %s\n", compare_memcmp(c, n_rows_a * n_cols_a, d, n_rows_a * n_cols_a) ? "\aDifferent" : "Same");
    printf("Compare (loop)   multiply_arrays_tiled() with multiply_arrays(): %s\n", compare(c, n_rows_a * n_cols_a, d, n_rows_a * n_cols_a) ? "\aDifferent" : "Same");

    if (scale == 1) {
        puts("");
        printf("Matrix A:\n");
        print(a, n_rows_a, n_cols_a);
        printf("Matrix B:\n");
        print(b, n_rows_a, n_cols_a);
        printf("Matrix C:\n");
        print(c, n_rows_a, n_cols_a);
        printf("Matrix D:\n");
        print(d, n_rows_a, n_cols_a);
    }

    double sum_c = sum_array(c, n_rows_a * n_cols_a);
    double sum_d = sum_array(d, n_rows_a * n_cols_a);

    puts("");
    printf("Sum of C: %f\n", sum_c);
    printf("Sum of D: %f\n", sum_d);
    printf("Compare sums of C and D: %s\n", compare_scalars(sum_c, sum_d) ? "\aDifferent" : "Same");
    puts("");

    free(a);
    free(b);
    free(c);
    free(d);
}

/* Test correctnes and speed of dot_simple, dot_simple_tiled, dot_faster and dot_faster_tiled */
void test_speed_2() {
    /* For measuring time */
    double t0, t1;

    const unsigned scale = 500;
    const unsigned n_rows_a = 4 * scale;
    const unsigned n_cols_a = 3 * scale;
    const unsigned n_rows_b = 3 * scale;
    const unsigned n_cols_b = 2 * scale;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_b * n_cols_b * sizeof(*b));
    double *c = malloc(n_rows_a * n_cols_b * sizeof(*c));
    double *d = malloc(n_rows_a * n_cols_b * sizeof(*d));
    double *e = malloc(n_rows_a * n_cols_b * sizeof(*e));
    double *f = malloc(n_rows_a * n_cols_b * sizeof(*f));

    if (!a || !b || !c || !d || !e) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    const int rand = 1;

    printf("\nRunning \"test_speed_2\" with %s data...\n\n", rand ? "random" : "sequential");

    if (rand) {
        init_rand(a, n_rows_a, n_cols_a);
        init_rand(b, n_rows_b, n_cols_b);
    }
    else {
        init_seq(a, n_rows_a, n_cols_a);
        init_seq(b, n_rows_b, n_cols_b);
    }

    const unsigned loop = 1u;

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        c = dot_simple(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c);
    }
    t1 = omp_get_wtime();
    printf("dot_simple():       Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        d = dot_simple_tiled(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, d);
    }
    t1 = omp_get_wtime();
    printf("dot_simple_tiled(): Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        e = dot_faster(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, e);
    }
    t1 = omp_get_wtime();
    printf("dot_faster():       Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        f = dot_faster_tiled(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, f);
    }
    t1 = omp_get_wtime();
    printf("dot_faster_tiled(): Elapsed time %.3f s\n", t1 - t0);

    puts("");
    printf("Compare (memcmp) dot_simple_tiled() with dot_simple(): %s\n", compare_memcmp(c, n_rows_a * n_cols_b, d, n_rows_a * n_cols_b) ? "\aDifferent" : "Same");
    printf("Compare (loop)   dot_simple_tiled() with dot_simple(): %s\n", compare(c, n_rows_a * n_cols_b, d, n_rows_a * n_cols_b) ? "\aDifferent" : "Same");
    printf("Compare (memcmp) dot_faster()       with dot_simple(): %s\n", compare_memcmp(c, n_rows_a * n_cols_b, e, n_rows_a * n_cols_b) ? "\aDifferent" : "Same");
    printf("Compare (loop)   dot_faster()       with dot_simple(): %s\n", compare(c, n_rows_a * n_cols_b, e, n_rows_a * n_cols_b) ? "\aDifferent" : "Same");
    printf("Compare (memcmp) dot_faster_tiled() with dot_simple(): %s\n", compare_memcmp(c, n_rows_a * n_cols_b, f, n_rows_a * n_cols_b) ? "\aDifferent" : "Same");
    printf("Compare (loop)   dot_faster_tiled() with dot_simple(): %s\n", compare(c, n_rows_a * n_cols_b, f, n_rows_a * n_cols_b) ? "\aDifferent" : "Same");

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
        printf("Matrix E:\n");
        print(e, n_rows_a, n_cols_b);
        printf("Matrix F:\n");
        print(f, n_rows_a, n_cols_b);
    }

    double sum_c = sum_array(c, n_rows_a * n_cols_b);
    double sum_d = sum_array(d, n_rows_a * n_cols_b);
    double sum_e = sum_array(e, n_rows_a * n_cols_b);
    double sum_f = sum_array(f, n_rows_a * n_cols_b);

    puts("");
    printf("Sum of C: %f\n", sum_c);
    printf("Sum of D: %f\n", sum_d);
    printf("Sum of E: %f\n", sum_e);
    printf("Sum of F: %f\n", sum_f);
    printf("Compare sums of C and D: %s\n", compare_scalars(sum_c, sum_d) ? "\aDifferent" : "Same");
    printf("Compare sums of C and E: %s\n", compare_scalars(sum_c, sum_e) ? "\aDifferent" : "Same");
    printf("Compare sums of C and F: %s\n", compare_scalars(sum_c, sum_f) ? "\aDifferent" : "Same");
    puts("");

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
}

/* Test correctnes and speed of transpose */
void test_speed_3() {
    /* For measuring time */
    double t0, t1;

    const unsigned scale = 500;
    const unsigned n_rows_a = 6 * scale;
    const unsigned n_cols_a = 5 * scale;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_a * n_cols_a * sizeof(*b));
    double *c = malloc(n_rows_a * n_cols_a * sizeof(*c));

    if (!(a && b && c)) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    const int rand = 1;

    printf("\nRunning \"test_speed_3\" with %s data...\n\n", rand ? "random" : "sequential");

    if (rand) {
        init_rand(a, n_rows_a, n_cols_a);
    }
    else {
        init_seq(a, n_rows_a, n_cols_a);
    }

    const unsigned loop = 100u;

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        b = transpose(a, n_rows_a, n_cols_a, b);
    }
    t1 = omp_get_wtime();
    printf("transpose():           Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        c = transpose_non_tiled(a, n_rows_a, n_cols_a, c);
    }
    t1 = omp_get_wtime();
    printf("transpose_non_tiled(): Elapsed time %.3f s\n", t1 - t0);

    if (scale == 1) {
        puts("");
        printf("Matrix A:\n");
        print(a, n_rows_a, n_cols_a);
        printf("Matrix B:\n");
        print(b, n_cols_a, n_rows_a);
        printf("Matrix C:\n");
        print(c, n_cols_a, n_rows_a);
    }

    puts("");
    printf("Compare (memcmp) transpose() with transpose_non_tiled(): %s\n", compare_memcmp(b, n_rows_a * n_cols_a, c, n_rows_a * n_cols_a) ? "\aDifferent" : "Same");
    printf("Compare (loop)   transpose() with transpose_non_tiled(): %s\n", compare(b, n_rows_a * n_cols_a, c, n_rows_a * n_cols_a) ? "\aDifferent" : "Same");

    double sum_a = sum_array(a, n_rows_a * n_cols_a);
    double sum_b = sum_array(b, n_rows_a * n_cols_a);
    double sum_c = sum_array(c, n_rows_a * n_cols_a);

    puts("");
    printf("Sum of A: %f\n", sum_a);
    printf("Sum of B: %f\n", sum_b);
    printf("Sum of C: %f\n", sum_c);
    printf("Compare sums of A and B: %s\n", compare_scalars(sum_a, sum_b) ? "\aDifferent" : "Same");
    printf("Compare sums of A and C: %s\n", compare_scalars(sum_a, sum_c) ? "\aDifferent" : "Same");
    puts("");

    free(a);
    free(b);
    free(c);
}

/* Test correctnes and speed of init and sum */
void test_speed_4() {
    /* For measuring time */
    double t0, t1;

    const unsigned scale = 1600;
    const unsigned n_rows_a = 6 * scale;
    const unsigned n_cols_a = 5 * scale;

    double *a = malloc(n_rows_a * n_cols_a * sizeof(*a));
    double *b = malloc(n_rows_a * n_cols_a * sizeof(*b));

    if (!(a && b)) {
        printf("Couldn't allocate memory!\n");
        system("pause");
        exit(-1);
    }

    int rand = 0;

    printf("\nRunning \"test_speed_4\" with %s data...\n\n", rand ? "random" : "sequential");

    const unsigned loop = 10u;

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        init_seq(a, n_rows_a, n_cols_a);
    }
    t1 = omp_get_wtime();
    printf("init_seq():       Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        init_seq_tiled(b, n_rows_a, n_cols_a);
    }
    t1 = omp_get_wtime();
    printf("init_seq_tiled(): Elapsed time %.3f s\n", t1 - t0);

    puts("");
    printf("Compare (memcmp) init_seq() with init_seq_tiled(): %s\n", compare_memcmp(a, n_rows_a * n_cols_a, b, n_rows_a * n_cols_a) ? "\aDifferent" : "Same");
    printf("Compare (loop)   init_seq() with init_seq_tiled(): %s\n", compare(a, n_rows_a * n_cols_a, b, n_rows_a * n_cols_a) ? "\aDifferent" : "Same");

    volatile double sum_a = sum_array(a, n_rows_a * n_cols_a);
    volatile double sum_b = sum_array(b, n_rows_a * n_cols_a);
    
    puts("");
    printf("Sum of A: %f\n", sum_a);
    printf("Sum of B: %f\n", sum_b);
    printf("Compare sums of A and B: %s\n", compare_scalars(sum_a, sum_b) ? "\aDifferent" : "Same");

    if (scale == 1) {
        puts("");
        printf("Matrix A:\n");
        print(a, n_rows_a, n_cols_a);
        printf("Matrix B:\n");
        print(b, n_rows_a, n_cols_a);
    }

    rand = 1;

    printf("\nRunning \"test_speed_4\" with %s data...\n\n", rand ? "random" : "sequential");

    sum_a = sum_b = 0;
    
    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        init_rand(a, n_rows_a, n_cols_a);
    }
    t1 = omp_get_wtime();
    printf("init_rand():       Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        init_rand_tiled(b, n_rows_a, n_cols_a);
    }
    t1 = omp_get_wtime();
    printf("init_rand_tiled(): Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        sum_a += sum_array(a, n_rows_a * n_cols_a);
    }
    t1 = omp_get_wtime();
    printf("sum_array():       Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    for (size_t i = 0; i < loop; i++) {
        sum_b += sum_array_tiled(a, n_rows_a * n_cols_a);
    }
    t1 = omp_get_wtime();
    printf("sum_array_tiled(): Elapsed time %.3f s\n", t1 - t0);

    puts("");
    printf("Sum of A: %f\n", sum_a / loop);
    printf("Sum of A: %f\n", sum_b / loop);
    printf("Compare sums of A and A: %s\n", compare_scalars(sum_a, sum_b) ? "\aDifferent" : "Same");
    puts("");

    free(a);
    free(b);
}
