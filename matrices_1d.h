#ifndef MATRICES_1D_INCLUDE
#define MATRICES_1D_INCLUDE

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* This can be used as type of our data. */
typedef double data_type;

/* Type of data is double, which is 8 bytes in size.
   Cache line size depends on CPU. On mine, it's 64 bytes, so it can hold 8 elements.
   Tile order is size of a dimension of a square matrix (a tile) to use
   when working with matrices, to speed operations up. So, one tile has
   TILE_ORDER * TILE_ORDER many elements. For example, if TILE_ORDER is 8,
   one tile will occupy exactly 8 cache lines.  */
#define TILE_ORDER 32

/* Tolerance for comparing doubles */
#define TOLERANCE 1.e-6

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/* Initializes vector or matrix with sequentially growing values. */
void init_seq(double *a, const unsigned n_rows_a, const unsigned n_cols_a);

/* Initializes vector or matrix, randomly. */
void init_rand(double *a, const unsigned n_rows_a, const unsigned n_cols_a);

/* Sum of an array */
double sum_array(const double *arr, const unsigned length);

/* Mean value of an array */
double mean(const double *arr, const unsigned length);

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
    It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
    of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
    The original matrix m stays intact. */
double *transpose(const double *m, const unsigned n_rows_m, const unsigned n_cols_m, double *t);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is by far the slowest version of the function, sequentially or parallely. */
double *dot_simple(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
    const double *b, const unsigned n_rows_b, const unsigned n_cols_b, double *c);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is a tiled version of the simple function, and it's much faster than it. */
double *dot_simple_tiled(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
    const double *b, const unsigned n_rows_b, const unsigned n_cols_b, double *c);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is a much faster version of the function.
 * It's the fastest one, sequential or Open MP. */
double *dot_faster(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
    const double *b, const unsigned n_rows_b, const unsigned n_cols_b, double *c);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This was supposed to be the fastest version of the function,
 * but it's similar in speed to dot_simple_tiled. */
double *dot_faster_tiled(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
    const double *b, const unsigned n_rows_b, const unsigned n_cols_b, double *c);

/*  Adds two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *add_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);

/*  Subtracts the second array from the first one, element-wise, and puts the result
    in an array that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *subtract_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);

/*  Multiplies two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it).
    Tiled version is slightly slower in Open MP, and evidently slower sequentially. */
double *multiply_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);
double *multiply_arrays_tiled(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);

/*  Divides two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *divide_arrays(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);

/*  Updates an array, element-wise, by adding another array to it.
    Takes both arrays in, and returns the updated one (the first one).
    The return value (address of the first array) doesn't have to be used.
    Arrays must be of the same length, or, the second one can be a scalar.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *add_update(double *a, const unsigned n_a, const double *b, const unsigned n_b);

/*  Compares two arrays element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    If an element of array a is greater than a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *greater_than(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);

/*  Compares two arrays element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    If an element of array a is equal to a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
double *equal(const double *a, const unsigned n_a, const double *b, const unsigned n_b, double *result);

/* Prints vector, or matrix. */
void print(const double *m, const unsigned n_rows_m, const unsigned n_cols_m);

/* Sequential function for comparing two arrays by using memcmp
   Returns 0 if contents of the arrays are the same; -1 or 1 otherwise. */
int compare_memcmp(const double *a, const unsigned n_a, const double *b, const unsigned n_b);

/* Sequential function for comparing two arrays by using a loop
   Returns 0 if contents of the arrays are the same; 1 otherwise. */
int compare(const double *a, const unsigned n_a, const double *b, const unsigned n_b);

/* Compares two scalars within a given TOLERANCE
   Returns 0 if contents of the arrays are the same; 1 otherwise. */
int compare_scalars(const double a, const double b);

#endif // ! MATRICES_1D_INCLUDE
