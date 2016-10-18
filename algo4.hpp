#include <iostream>
#include <math.h>
#include <random> // Normal distribution

#include "coco.h"
/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

/**
 * A ES(mu/mu,lambda)-search-path optimizer that can be used
 * for single-objective optimization.
 *
 * @param evaluate The evaluation function used to evaluate the solutions.
 * @param dimension The number of variables (= n in our case)
 * @param number_of_objectives The number of objectives. = 1 in our case
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param max_budget The maximal number of evaluations.
 */
void algo4(evaluate_function_t evaluate,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t max_budget);
/**
 * Euclidian norm
 *
 * @param m pointer to a vector
 * @param dimension size of the vector
 * @return Euclidian norm of the vector
 */
double normE(double *m, size_t dimension);

/**
 * Destructor for handcrafted matrices
 *
 * @param mat pointer to a matrix
 * @param nb_rows Number of rows of the matrix
 */
void freeMatrix(double** mat, size_t nb_rows);
