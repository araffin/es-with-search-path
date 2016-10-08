#include <iostream>
#include <math.h>
#include <random> // Normal distribution

#include "coco.h"
/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

void algo4(evaluate_function_t evaluate,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t max_budget);

void select_mu_best(double mu, double lambda, double** X_k, double** Z, evaluate_function_t evaluate);

void printArray(double* array, size_t n);
void printMatrix(double** array, size_t n, size_t m);
void elementProduct(double* a, double* b, double* result, size_t n);
void arraySum(double* a, double* b, double* result, size_t n);
void sumVectors(double** mat, double* result, size_t k, size_t n);
void normalMatrix(double** mat, std::normal_distribution<> N, std::mt19937 gen,
                  size_t nb_rows, size_t nb_cols);
void freeMatrix(double** mat, size_t nb_rows);
