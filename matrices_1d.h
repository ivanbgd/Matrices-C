#ifndef MATRICES_1D_HEADER
#define MATRICES_1D_HEADER

/* All functions use restricted pointers, so care should be taken
 * to make sure that arrays that they point to do not overlap, if
 * we want to modify them inside of the functions.
 * On the other hand, it's easy to change type of the pointers
 * from restricted to non-restricted versions, by using definitions
 * given at the beginning of this header file, if necessary. */

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* This is used as type of our data. */
typedef double data_t;

/* Pointers to data */
typedef data_t * data_ptr_t;
typedef const data_t * cdata_ptr_t;
typedef data_t * __restrict data_ptr_res_t;
typedef const data_t * __restrict cdata_ptr_res_t;

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

/* Initializes vector or matrix with sequentially growing data_t values, starting from 0.
   Tiled variant is slightly slower; larger diiference is sequentially. */
void init_seq(data_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a);
void init_seq_tiled(data_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a);

/* Initializes vector or matrix, with random data_t values in the range [0, 1].
   Lot slower than init_seq(), which is expected, since it calls rand().
   Tiled variant is slightly slower; larger diiference is sequentially. */
void init_rand(data_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a);
void init_rand_tiled(data_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a);

/* Sum of an array
   Tiled and non-tiled variants are of similar speed,
   though larger diiference is sequentially. */
data_t sum_array(cdata_ptr_res_t arr, const unsigned length);
data_t sum_array_tiled(cdata_ptr_res_t arr, const unsigned length);

/* Mean value of an array */
data_t mean(cdata_ptr_res_t arr, const unsigned length);

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
    It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
    of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
    The original matrix m stays intact. */
data_ptr_res_t transpose(cdata_ptr_res_t m, const unsigned n_rows_m, const unsigned n_cols_m, data_ptr_res_t t);
data_ptr_res_t transpose_non_tiled(cdata_ptr_res_t m, const unsigned n_rows_m, const unsigned n_cols_m, data_ptr_res_t t);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is by far the slowest version of the function, sequentially or parallely. */
data_ptr_res_t dot_simple(cdata_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a, \
    cdata_ptr_res_t b, const unsigned n_rows_b, const unsigned n_cols_b, data_ptr_res_t c);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is a tiled version of the simple function, and it's much faster than it. */
data_ptr_res_t dot_simple_tiled(cdata_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a, \
    cdata_ptr_res_t b, const unsigned n_rows_b, const unsigned n_cols_b, data_ptr_res_t c);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * This is a much faster version of the function.
 * Uses more memory than the simple version, for transposing matrix b, which
 * can be a problem if the matrix is large - there might not be enough memory.
 * It's the fastest one, sequential or Open MP, if we optimize for speed. */
data_ptr_res_t dot_faster(cdata_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a, \
    cdata_ptr_res_t b, const unsigned n_rows_b, const unsigned n_cols_b, data_ptr_res_t c);

/* Dot product of two arrays, a and b, or matrix product
 * Returns an array that's passed in as the last argument, c.
 * Uses more memory than the simple version, for transposing matrix b, which
 * can be a problem if the matrix is large - there might not be enough memory.
 * This was supposed to be the fastest version of the function,
 * but it's similar in speed to dot_simple_tiled, if we optimize for speed.
 * But, if we optimize for the smallest code, this version is the fastest,
 * though, not faster than dot_faster optimized for speed, but of the same speed as it. */
data_ptr_res_t dot_faster_tiled(cdata_ptr_res_t a, const unsigned n_rows_a, const unsigned n_cols_a, \
    cdata_ptr_res_t b, const unsigned n_rows_b, const unsigned n_cols_b, data_ptr_res_t c);

/*  Adds two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
data_ptr_res_t add_arrays(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);

/*  Subtracts the second array from the first one, element-wise, and puts the result
    in an array that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
data_ptr_res_t subtract_arrays(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);

/*  Multiplies two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it).
    Tiled version is slightly slower in Open MP, and evidently slower sequentially. */
data_ptr_res_t multiply_arrays(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);
data_ptr_res_t multiply_arrays_tiled(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);

/*  Divides two arrays, element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
data_ptr_res_t divide_arrays(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);

/*  Updates an array, element-wise, by adding another array to it.
    Takes both arrays in, and returns the updated one (the first one).
    The return value (address of the first array) doesn't have to be used.
    Arrays must be of the same length, or, the second one can be a scalar.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
data_ptr_res_t add_update(data_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b);

/*  Compares two arrays element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    If an element of array a is greater than a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
data_ptr_res_t greater_than(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);

/*  Compares two arrays element-wise, and puts the result in an array
    that is passed in as the last argument, and also returns it.
    If an element of array a is equal to a corresponding element of
    array b, the resulting array will have 1.0 in that position;
    it will have 0.0 otherwise.
    Arrays must be of the same length, or, one of them, or both, can be scalars.
    Use 0 as the length of a scalar, and pass its address in (a pointer to it). */
data_ptr_res_t equal(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b, data_ptr_res_t result);

/* Prints vector, or matrix. */
void print(cdata_ptr_res_t m, const unsigned n_rows_m, const unsigned n_cols_m);

/* Sequential function for comparing two arrays by using memcmp
   Returns 0 if contents of the arrays are the same; -1 or 1 otherwise. */
int compare_memcmp(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b);

/* Sequential function for comparing two arrays by using a loop
   Returns 0 if contents of the arrays are the same; 1 otherwise. */
int compare(cdata_ptr_res_t a, const unsigned n_a, cdata_ptr_res_t b, const unsigned n_b);

/* Compares two scalars within a given TOLERANCE
   Returns 0 if contents of the arrays are the same; 1 otherwise. */
int compare_scalars(const data_t a, const data_t b);

#endif // ! MATRICES_1D_HEADER
