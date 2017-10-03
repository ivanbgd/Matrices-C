#ifndef TESTS_HEADER
#define TESTS_HEADER

#include "matrices_1d.h"

/* Function that calls other test functions */
void test();

/* Simple basic tests */
void test_0();

/* Test speed of multiply_arrays */
void test_speed_1();

/* Test correctnes and speed of dot_simple, dot_simple_tiled, dot_faster and dot_faster_tiled */
void test_speed_2();

/* Test correctnes and speed of transpose */
void test_speed_3();

/* Test speed of init and sum
   Try with different TILE_ORDER, starting with 1. */
void test_speed_4();

#endif  // TESTS_HEADER
