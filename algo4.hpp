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