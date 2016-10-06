
/************************************************************************
 * WARNING
 *
 * This file is an auto-generated amalgamation. Any changes made to this
 * file will be lost when it is regenerated!
 ************************************************************************/

#line 1 "code-experiments/src/coco_random.c"
/**
 * @file coco_random.c
 * @brief Definitions of functions regarding COCO random numbers.
 *
 * @note This file contains non-C89-standard types (such as uint32_t and uint64_t), which should
 * eventually be fixed.
 */

#include <math.h>

#line 1 "code-experiments/src/coco.h"
/**
 * @file coco.h
 * @brief All public COCO functions, constants and variables are defined in this file.
 *
 * It is the authoritative reference, if any function deviates from the documented behavior it is considered
 * a bug. See the function definitions for their detailed descriptions.
 */
 
#ifndef __COCO_H__
#define __COCO_H__

#include <stddef.h>

/* Definitions of some 32 and 64-bit types (used by the random number generator) */
#ifdef _MSC_VER
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

/* Include definition for NAN among other things */
#include <math.h>
#include <float.h>
#ifndef NAN
/** @brief Definition of NAN to be used only if undefined by the included headers */
#define NAN 8.8888e88
#endif
#ifndef isnan
/** @brief Definition of isnan to be used only if undefined by the included headers */
#define isnan(x) (0)
#endif
#ifndef INFINITY
/** @brief Definition of INFINITY to be used only if undefined by the included headers */
#define INFINITY 1e22
#endif
#ifndef isinf
/** @brief Definition of isinf to be used only if undefined by the included headers */
#define isinf(x) (0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief COCO's version.
 *
 * Automatically updated by do.py.
 */
/**@{*/
static const char coco_version[32] = "1.2";
/**@}*/

/***********************************************************************************************************/
/**
 * @brief COCO's own pi constant. Simplifies the case, when the value of pi changes.
 */
/**@{*/
static const double coco_pi = 3.14159265358979323846;
static const double coco_two_pi = 2.0 * 3.14159265358979323846;
/**@}*/

/***********************************************************************************************************/

/** @brief Logging level type. */
typedef enum {
  COCO_ERROR,     /**< @brief only error messages are output */
  COCO_WARNING,   /**< @brief error and warning messages are output */
  COCO_INFO,      /**< @brief error, warning and info messages are output */
  COCO_DEBUG      /**< @brief error, warning, info and debug messages are output */
} coco_log_level_type_e;

/***********************************************************************************************************/

/** @brief Structure containing a COCO problem. */
struct coco_problem_s;

/**
 * @brief The COCO problem type.
 *
 * See coco_problem_s for more information on its fields. */
typedef struct coco_problem_s coco_problem_t;

/** @brief Structure containing a COCO suite. */
struct coco_suite_s;

/**
 * @brief The COCO suite type.
 *
 * See coco_suite_s for more information on its fields. */
typedef struct coco_suite_s coco_suite_t;

/** @brief Structure containing a COCO observer. */
struct coco_observer_s;

/**
 * @brief The COCO observer type.
 *
 * See coco_observer_s for more information on its fields. */
typedef struct coco_observer_s coco_observer_t;

/** @brief Structure containing a COCO archive. */
struct coco_archive_s;

/**
 * @brief The COCO archive type.
 *
 * See coco_archive_s for more information on its fields. */
typedef struct coco_archive_s coco_archive_t;

/** @brief Structure containing a COCO random state. */
struct coco_random_state_s;

/**
 * @brief The COCO random state type.
 *
 * See coco_random_state_s for more information on its fields. */
typedef struct coco_random_state_s coco_random_state_t;

/***********************************************************************************************************/

/**
 * @name Methods regarding COCO suite
 */
/**@{*/

/**
 * @brief Constructs a COCO suite.
 */
coco_suite_t *coco_suite(const char *suite_name, const char *suite_instance, const char *suite_options);

/**
 * @brief Frees the given suite.
 */
void coco_suite_free(coco_suite_t *suite);

/**
 * @brief Returns the next (observed) problem of the suite or NULL if there is no next problem left.
 */
coco_problem_t *coco_suite_get_next_problem(coco_suite_t *suite, coco_observer_t *observer);

/**
 * @brief Returns the problem of the suite defined by problem_index.
 */
coco_problem_t *coco_suite_get_problem(coco_suite_t *suite, const size_t problem_index);

/**
 * @brief Returns the first problem of the suite defined by function, dimension and instance numbers.
 */
coco_problem_t *coco_suite_get_problem_by_function_dimension_instance(coco_suite_t *suite,
                                                                      const size_t function,
                                                                      const size_t dimension,
                                                                      const size_t instance);

/**
 * @brief Returns the number of problems in the given suite.
 */
size_t coco_suite_get_number_of_problems(const coco_suite_t *suite);

/**
 * @brief Returns the function number in the suite in position function_idx (counting from 0).
 */
size_t coco_suite_get_function_from_function_index(const coco_suite_t *suite, const size_t function_idx);

/**
 * @brief Returns the dimension number in the suite in position dimension_idx (counting from 0).
 */
size_t coco_suite_get_dimension_from_dimension_index(const coco_suite_t *suite, const size_t dimension_idx);

/**
 * @brief Returns the instance number in the suite in position instance_idx (counting from 0).
 */
size_t coco_suite_get_instance_from_instance_index(const coco_suite_t *suite, const size_t instance_idx);
/**@}*/

/**
 * @name Encoding/decoding problem index
 *
 * General schema for encoding/decoding a problem index. Note that the index depends on the number of
 * instances a suite is defined with (it should be called a suite-instance-depending index...).
 * Also, while functions, instances and dimensions start from 1, function_idx, instance_idx and dimension_idx
 * as well as suite_dep_index start from 0!
 *
 * Showing an example with 2 dimensions (2, 3), 5 instances (6, 7, 8, 9, 10) and 2 functions (1, 2):
 *
   \verbatim
   index | instance | function | dimension
   ------+----------+----------+-----------
       0 |        6 |        1 |         2
       1 |        7 |        1 |         2
       2 |        8 |        1 |         2
       3 |        9 |        1 |         2
       4 |       10 |        1 |         2
       5 |        6 |        2 |         2
       6 |        7 |        2 |         2
       7 |        8 |        2 |         2
       8 |        9 |        2 |         2
       9 |       10 |        2 |         2
      10 |        6 |        1 |         3
      11 |        7 |        1 |         3
      12 |        8 |        1 |         3
      13 |        9 |        1 |         3
      14 |       10 |        1 |         3
      15 |        6 |        2 |         2
      16 |        7 |        2 |         3
      17 |        8 |        2 |         3
      18 |        9 |        2 |         3
      19 |       10 |        2 |         3

   index | instance_idx | function_idx | dimension_idx
   ------+--------------+--------------+---------------
       0 |            0 |            0 |             0
       1 |            1 |            0 |             0
       2 |            2 |            0 |             0
       3 |            3 |            0 |             0
       4 |            4 |            0 |             0
       5 |            0 |            1 |             0
       6 |            1 |            1 |             0
       7 |            2 |            1 |             0
       8 |            3 |            1 |             0
       9 |            4 |            1 |             0
      10 |            0 |            0 |             1
      11 |            1 |            0 |             1
      12 |            2 |            0 |             1
      13 |            3 |            0 |             1
      14 |            4 |            0 |             1
      15 |            0 |            1 |             1
      16 |            1 |            1 |             1
      17 |            2 |            1 |             1
      18 |            3 |            1 |             1
      19 |            4 |            1 |             1
   \endverbatim
 */
/**@{*/
/**
 * @brief Computes the index of the problem in the suite that corresponds to the given function, dimension
 * and instance indices.
 */
size_t coco_suite_encode_problem_index(const coco_suite_t *suite,
                                       const size_t function_idx,
                                       const size_t dimension_idx,
                                       const size_t instance_idx);

/**
 * @brief Computes the function, dimension and instance indexes of the problem with problem_index in the
 * given suite.
 */
void coco_suite_decode_problem_index(const coco_suite_t *suite,
                                     const size_t problem_index,
                                     size_t *function_idx,
                                     size_t *dimension_idx,
                                     size_t *instance_idx);
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding COCO observer
 */
/**@{*/
/**
 * @brief Constructs a COCO observer.
 */
coco_observer_t *coco_observer(const char *observer_name, const char *options);

/**
 * @brief Frees the given observer.
 */
void coco_observer_free(coco_observer_t *observer);

/**
 * @brief Adds an observer to the given problem.
 */
coco_problem_t *coco_problem_add_observer(coco_problem_t *problem, coco_observer_t *observer);

/**
 * @brief Removes an observer from the given problem.
 */
coco_problem_t *coco_problem_remove_observer(coco_problem_t *problem, coco_observer_t *observer);

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding COCO problem
 */
/**@{*/
/**
 * @brief Evaluates the problem function in point x and save the result in y.
 */
void coco_evaluate_function(coco_problem_t *problem, const double *x, double *y);

/**
 * @brief Evaluates the problem constraints in point x and save the result in y.
 */
void coco_evaluate_constraint(coco_problem_t *problem, const double *x, double *y);

/**
 * @brief Recommends a solution as the current best guesses to the problem.
 */
void coco_recommend_solution(coco_problem_t *problem, const double *x);

/**
 * @brief Frees the given problem.
 */
void coco_problem_free(coco_problem_t *problem);

/**
 * @brief Returns the name of the problem.
 */
const char *coco_problem_get_name(const coco_problem_t *problem);

/**
 * @brief Returns the ID of the problem.
 */
const char *coco_problem_get_id(const coco_problem_t *problem);

/**
 * @brief Returns the number of variables i.e. the dimension of the problem.
 */
size_t coco_problem_get_dimension(const coco_problem_t *problem);

/**
 * @brief Returns the number of objectives of the problem.
 */
size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem);

/**
 * @brief Returns the number of constraints of the problem.
 */
size_t coco_problem_get_number_of_constraints(const coco_problem_t *problem);

/**
 * @brief Returns the number of evaluations done on the problem.
 */
size_t coco_problem_get_evaluations(const coco_problem_t *problem);

/**
 * @brief Returns 1 if the final target was hit, 0 otherwise.
 */
int coco_problem_final_target_hit(const coco_problem_t *problem);

/**
 * @brief Returns the best observed value for the first objective.
 */
double coco_problem_get_best_observed_fvalue1(const coco_problem_t *problem);

/**
 * @brief Returns the target value for the first objective.
 */
double depreciated_coco_problem_get_final_target_fvalue1(const coco_problem_t *problem);

/**
 * @brief Returns a vector of size 'dimension' with lower bounds of the region of interest in
 * the decision space.
 */
const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem);

/**
 * @brief Returns a vector of size 'dimension' with upper bounds of the region of interest in
 * the decision space.
 */
const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem);

/**
 * @brief Returns the problem_index of the problem in its current suite.
 */
size_t coco_problem_get_suite_dep_index(const coco_problem_t *problem);

/**
 * @brief Returns an initial solution, i.e. a feasible variable setting, to the problem.
 */
void coco_problem_get_initial_solution(const coco_problem_t *problem, double *initial_solution);
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding random numbers
 */
/**@{*/

/**
 * @brief Creates and returns a new random number state using the given seed.
 */
coco_random_state_t *coco_random_new(uint32_t seed);

/**
 * @brief Frees all memory associated with the random state.
 */
void coco_random_free(coco_random_state_t *state);

/**
 * @brief Returns one uniform [0, 1) random value from the random number generator associated with the given
 * state.
 */
double coco_random_uniform(coco_random_state_t *state);

/**
 * @brief Generates an approximately normal random number.
 */
double coco_random_normal(coco_random_state_t *state);
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods managing memory
 */
/**@{*/
/**
 * @brief Safe memory allocation that either succeeds or triggers a coco_error.
 */
void *coco_allocate_memory(const size_t size);

/**
 * @brief Safe memory allocation for a vector of doubles that either succeeds or triggers a coco_error.
 */
double *coco_allocate_vector(const size_t size);

/**
 * @brief Frees the allocated memory.
 */
void coco_free_memory(void *data);
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding COCO messages
 */
/**@{*/
/**
 * @brief Signals a fatal error.
 */
void coco_error(const char *message, ...);

/**
 * @brief Warns about error conditions.
 */
void coco_warning(const char *message, ...);

/**
 * @brief Outputs some information.
 */
void coco_info(const char *message, ...);

/**
 * @brief Prints only the given message without any prefix and new line.
 *
 * A function similar to coco_info but producing no additional text than
 * the given message.
 *
 * The output is only produced if coco_log_level >= COCO_INFO.
 */
void coco_info_partial(const char *message, ...);

/**
 * @brief Outputs detailed information usually used for debugging.
 */
void coco_debug(const char *message, ...);

/**
 * @brief Sets the COCO log level to the given value and returns the previous value of the log level.
 */
const char *coco_set_log_level(const char *level);
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding COCO archives and log files (used when pre-processing MO data)
 */
/**@{*/

/**
 * @brief Constructs a COCO archive.
 */
coco_archive_t *coco_archive(const char *suite_name,
                             const size_t function,
                             const size_t dimension,
                             const size_t instance);
/**
 * @brief Adds a solution with objectives (y1, y2) to the archive if none of the existing solutions in the
 * archive dominates it. In this case, returns 1, otherwise the archive is not updated and the method
 * returns 0.
 */
int coco_archive_add_solution(coco_archive_t *archive, const double y1, const double y2, const char *text);

/**
 * @brief Returns the number of (non-dominated) solutions in the archive (computed first, if needed).
 */
size_t coco_archive_get_number_of_solutions(coco_archive_t *archive);

/**
 * @brief Returns the hypervolume of the archive (computed first, if needed).
 */
double coco_archive_get_hypervolume(coco_archive_t *archive);

/**
 * @brief Returns the text of the next (non-dominated) solution in the archive and "" when there are no
 * solutions left. The first two solutions are always the extreme ones.
 */
const char *coco_archive_get_next_solution_text(coco_archive_t *archive);

/**
 * @brief Frees the archive.
 */
void coco_archive_free(coco_archive_t *archive);

/**
 * @brief Feeds the solution to the bi-objective logger for logger output reconstruction purposes.
 */
int coco_logger_biobj_feed_solution(coco_problem_t *problem, const size_t evaluation, const double *y);
/**@}*/

/***********************************************************************************************************/

/**
 * @name Other useful methods
 */
/**@{*/
/**
 * @brief Removes the given directory and all its contents.
 */
int coco_remove_directory(const char *path);

/**
 * @brief Formatted string duplication.
 */
char *coco_strdupf(const char *str, ...);
/**@}*/

/***********************************************************************************************************/

#ifdef __cplusplus
}
#endif
#endif
#line 12 "code-experiments/src/coco_random.c"

#define COCO_NORMAL_POLAR /* Use polar transformation method */

#define COCO_SHORT_LAG 273
#define COCO_LONG_LAG 607

/**
 * @brief A structure containing the state of the COCO random generator.
 */
struct coco_random_state_s {
  double x[COCO_LONG_LAG];
  size_t index;
};

/**
 * @brief A lagged Fibonacci random number generator.
 *
 * This generator is nice because it is reasonably small and directly generates double values. The chosen
 * lags (607 and 273) lead to a generator with a period in excess of 2^607-1.
 */
static void coco_random_generate(coco_random_state_t *state) {
  size_t i;
  for (i = 0; i < COCO_SHORT_LAG; ++i) {
    double t = state->x[i] + state->x[i + (COCO_LONG_LAG - COCO_SHORT_LAG)];
    if (t >= 1.0)
      t -= 1.0;
    state->x[i] = t;
  }
  for (i = COCO_SHORT_LAG; i < COCO_LONG_LAG; ++i) {
    double t = state->x[i] + state->x[i - COCO_SHORT_LAG];
    if (t >= 1.0)
      t -= 1.0;
    state->x[i] = t;
  }
  state->index = 0;
}

coco_random_state_t *coco_random_new(uint32_t seed) {
  coco_random_state_t *state = (coco_random_state_t *) coco_allocate_memory(sizeof(*state));
  size_t i;
  /* Expand seed to fill initial state array. */
  for (i = 0; i < COCO_LONG_LAG; ++i) {
    /* Uses uint64_t to silence the compiler ("shift count negative or too big, undefined behavior" warning) */
    state->x[i] = ((double) seed) / (double) (((uint64_t) 1UL << 32) - 1);
    /* Advance seed based on simple RNG from TAOCP */
    seed = (uint32_t) 1812433253UL * (seed ^ (seed >> 30)) + ((uint32_t) i + 1);
  }
  state->index = 0;
  return state;
}

void coco_random_free(coco_random_state_t *state) {
  coco_free_memory(state);
}

double coco_random_uniform(coco_random_state_t *state) {
  /* If we have consumed all random numbers in our archive, it is time to run the actual generator for one
   * iteration to refill the state with 'LONG_LAG' new values. */
  if (state->index >= COCO_LONG_LAG)
    coco_random_generate(state);
  return state->x[state->index++];
}

/**
 * Instead of using the (expensive) polar method, we may cheat and abuse the central limit theorem. The sum
 * of 12 uniform random values has mean 6, variance 1 and is approximately normal. Subtract 6 and you get
 * an approximately N(0, 1) random number.
 */
double coco_random_normal(coco_random_state_t *state) {
  double normal;
#ifdef COCO_NORMAL_POLAR
  const double u1 = coco_random_uniform(state);
  const double u2 = coco_random_uniform(state);
  normal = sqrt(-2 * log(u1)) * cos(2 * coco_pi * u2);
#else
  int i;
  normal = 0.0;
  for (i = 0; i < 12; ++i) {
    normal += coco_random_uniform(state);
  }
  normal -= 6.0;
#endif
  return normal;
}

/* Be hygienic (for amalgamation) and undef lags. */
#undef COCO_SHORT_LAG
#undef COCO_LONG_LAG
#line 1 "code-experiments/src/coco_suite.c"
/**
 * @file coco_suite.c
 * @brief Definitions of functions regarding COCO suites.
 *
 * When a new suite is added, the functions coco_suite_intialize, coco_suite_get_instances_by_year and
 * coco_suite_get_problem_from_indices need to be updated.
 *
 * @see <a href="index.html">Instructions</a> on how to write new test functions and combine them into test
 * suites.
 */

#include <time.h>

#line 15 "code-experiments/src/coco_suite.c"
#line 1 "code-experiments/src/coco_internal.h"
/**
 * @file coco_internal.h
 * @brief Definitions of internal COCO structures and typedefs.
 *
 * These are used throughout the COCO code base but should not be used by any external code.
 */

#ifndef __COCO_INTERNAL__
#define __COCO_INTERNAL__

#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************************************************************/
/**
 * @brief The data free function type.
 *
 * This is a template for functions that free the contents of data (used to free the contents of data
 * fields in coco_problem, coco_suite and coco_observer).
 */
typedef void (*coco_data_free_function_t)(void *data);

/**
 * @brief The problem free function type.
 *
 * This is a template for functions that free the problem structure.
 */
typedef void (*coco_problem_free_function_t)(coco_problem_t *problem);

/**
 * @brief The initial solution function type.
 *
 * This is a template for functions that return an initial solution of the problem.
 */
typedef void (*coco_initial_solution_function_t)(const coco_problem_t *problem, double *y);

/**
 * @brief The evaluate function type.
 *
 * This is a template for functions that perform an evaluation of the problem (to evaluate the problem
 * function, the problems constraints etc.).
 */
typedef void (*coco_evaluate_function_t)(coco_problem_t *problem, const double *x, double *y);

/**
 * @brief The recommend solutions function type.
 *
 * This is a template for functions that log a recommended solution.
 */
typedef void (*coco_recommend_function_t)(coco_problem_t *problem, const double *x);

/**
 * @brief The allocate logger function type.
 *
 * This is a template for functions that allocate a logger (wrap a logger around the given problem and return
 * the wrapped problem).
 */
typedef coco_problem_t *(*coco_logger_allocate_function_t)(coco_observer_t *observer,
                                                           coco_problem_t *problem);

/**
 * @brief The free logger function type.
 *
 * This is a template for functions that free a logger.
 */
typedef void (*coco_logger_free_function_t)(void *logger);


/**
 * @brief The transformed COCO problem data type.
 *
 * This is a type of a generic structure for a transformed ("outer") coco_problem. It makes possible the
 * wrapping of problems as layers of an onion. Initialized in the coco_problem_transformed_allocate function,
 * it makes the current ("outer") transformed problem a "derived problem class", which inherits from the
 * "inner" problem, the "base class".
 *
 * From the perspective of the inner problem:
 * - data holds the meta-information to administer the inheritance
 * - data->data holds the additional fields of the derived class (the outer problem)
 * - data->inner_problem points to the inner problem (now we have a linked list)
 */
typedef struct {
  coco_problem_t *inner_problem;                  /**< @brief Pointer to the inner problem */
  void *data;                                     /**< @brief Pointer to data, which enables further
                                                  wrapping of the problem */
  coco_data_free_function_t data_free_function;   /**< @brief Function to free the contents of data */
} coco_problem_transformed_data_t;

/**
 * @brief The stacked COCO problem data type.
 *
 * This is a type of a structure used when stacking two problems (especially useful for constructing
 * multi-objective problems).
 */
typedef struct {
  coco_problem_t *problem1; /**< @brief Pointer to the first problem (objective) */
  coco_problem_t *problem2; /**< @brief Pointer to the second problem (objective) */
} coco_problem_stacked_data_t;

/**
 * @brief The option keys data type.
 *
 * This is a type of a structure used to contain a set of known option keys (used by suites and observers).
 */
typedef struct {
  size_t count;  /**< @brief Number of option keys */
  char **keys;   /**< @brief Pointer to option keys */
} coco_option_keys_t;


/***********************************************************************************************************/

/**
 * @brief The COCO problem structure.
 *
 * This is one of the main structures in COCO. It contains information about a problem to be optimized. The
 * problems can be wrapped around each other (similar to the onion layers) by means of the data field and
 * the coco_problem_transformed_data_t structure creating some kind of "object inheritance". Even the logger
 * is considered as just another coco_problem instance wrapped around the original problem.
 */
struct coco_problem_s {

  coco_initial_solution_function_t initial_solution;  /**< @brief  The function for creating an initial solution. */
  coco_evaluate_function_t evaluate_function;         /**< @brief  The function for evaluating the problem. */
  coco_evaluate_function_t evaluate_constraint;       /**< @brief  The function for evaluating the constraints. */
  coco_recommend_function_t recommend_solution;       /**< @brief  The function for recommending a solution. */
  coco_problem_free_function_t problem_free_function; /**< @brief  The function for freeing this problem. */

  size_t number_of_variables;          /**< @brief Number of variables expected by the function, i.e.
                                       problem dimension */
  size_t number_of_objectives;         /**< @brief Number of objectives. */
  size_t number_of_constraints;        /**< @brief Number of constraints. */

  double *smallest_values_of_interest; /**< @brief The lower bounds of the ROI in the decision space. */
  double *largest_values_of_interest;  /**< @brief The upper bounds of the ROI in the decision space. */

  double *best_value;                  /**< @brief Optimal (smallest) function value */
  double *nadir_value;                 /**< @brief The nadir point (defined when number_of_objectives > 1) */
  double *best_parameter;              /**< @brief Optimal decision vector (defined only when unique) */

  char *problem_name;                  /**< @brief Problem name. */
  char *problem_id;                    /**< @brief Problem ID (unique in the containing suite) */
  char *problem_type;                  /**< @brief Problem type */

  size_t evaluations;                  /**< @brief Number of evaluations performed on the problem. */

  /* Convenience fields for output generation */

  double final_target_delta[1];        /**< @brief Final target delta. */
  double best_observed_fvalue[1];      /**< @brief The best observed value so far. */
  size_t best_observed_evaluation[1];  /**< @brief The evaluation when the best value so far was achieved. */

  /* Fields depending on the containing benchmark suite */

  coco_suite_t *suite;                 /**< @brief Pointer to the containing suite (NULL if not given) */
  size_t suite_dep_index;              /**< @brief Suite-depending problem index (starting from 0) */
  size_t suite_dep_function;           /**< @brief Suite-depending function */
  size_t suite_dep_instance;           /**< @brief Suite-depending instance */

  void *data;                          /**< @brief Pointer to a data instance @see coco_problem_transformed_data_t */
};

/**
 * @brief The COCO observer structure.
 *
 * An observer observes the whole benchmark process. It is independent of suites and problems. Each time a
 * new problem of the suite is being observed, the observer initializes a new logger (wraps the observed
 * problem with the corresponding logger).
 */
struct coco_observer_s {

  int is_active;             /**< @brief Whether the observer is active (the logger will log some output). */
  char *observer_name;       /**< @brief Name of the observer for identification purposes. */
  char *result_folder;       /**< @brief Name of the result folder. */
  char *algorithm_name;      /**< @brief Name of the algorithm to be used in logger output. */
  char *algorithm_info;      /**< @brief Additional information on the algorithm to be used in logger output. */
  size_t number_target_triggers;
                             /**< @brief The number of targets between each 10**i and 10**(i+1). */
  double target_precision;   /**< @brief The minimal precision used for targets. */
  size_t number_evaluation_triggers;
                             /**< @brief The number of triggers between each 10**i and 10**(i+1) evaluation number. */
  char *base_evaluation_triggers;
                             /**< @brief The "base evaluations" used to evaluations that trigger logging. */
  int precision_x;           /**< @brief Output precision for decision variables. */
  int precision_f;           /**< @brief Output precision for function values. */
  void *data;                /**< @brief Void pointer that can be used to point to data specific to an observer. */

  coco_data_free_function_t data_free_function;             /**< @brief  The function for freeing this observer. */
  coco_logger_allocate_function_t logger_allocate_function; /**< @brief  The function for allocating the logger. */
  coco_logger_free_function_t logger_free_function;         /**< @brief  The function for freeing the logger. */
};

/**
 * @brief The COCO suite structure.
 *
 * A suite is a collection of problems constructed by a Cartesian product of the suite's optimization
 * functions, dimensions and instances. The functions and dimensions are fixed for a suite with some name,
 * while the instances are defined dynamically. The suite can be filtered - only the chosen functions,
 * dimensions and instances will be taken into account when iterating through the suite.
 */
struct coco_suite_s {

  char *suite_name;                /**< @brief Name of the suite. */

  size_t number_of_dimensions;     /**< @brief Number of dimensions contained in the suite. */
  size_t *dimensions;              /**< @brief The dimensions contained in the suite. */

  size_t number_of_functions;      /**< @brief Number of functions contained in the suite. */
  size_t *functions;               /**< @brief The functions contained in the suite. */

  size_t number_of_instances;      /**< @brief Number of instances contained in the suite. */
  char *default_instances;         /**< @brief The instances contained in the suite by default. */
  size_t *instances;               /**< @brief The instances contained in the suite. */

  coco_problem_t *current_problem; /**< @brief Pointer to the currently tackled problem. */
  long current_dimension_idx;      /**< @brief The dimension index of the currently tackled problem. */
  long current_function_idx;       /**< @brief The function index of the currently tackled problem. */
  long current_instance_idx;       /**< @brief The instance index of the currently tackled problem. */

  void *data;                      /**< @brief Void pointer that can be used to point to data specific to a suite. */

  coco_data_free_function_t data_free_function; /**< @brief The function for freeing this suite. */

};

#ifdef __cplusplus
}
#endif
#endif

#line 16 "code-experiments/src/coco_suite.c"
#line 1 "code-experiments/src/coco_utilities.c"
/**
 * @file coco_utilities.c
 * @brief Definitions of miscellaneous functions used throughout the COCO framework.
 */

#line 1 "code-experiments/src/coco_platform.h"
/**
 * @file coco_platform.h
 * @brief Automatic platform-dependent configuration of the COCO framework.
 *
 * Some platforms and standard conforming compilers require extra defines or includes to provide some
 * functionality.
 *
 * Because most feature defines need to be set before the first system header is included and we do not
 * know when a system header is included for the first time in the amalgamation, all internal files
 * that need these definitions should include this file before any system headers.
 */

#ifndef __COCO_PLATFORM__ 
#define __COCO_PLATFORM__

#include <stddef.h>

/* Definitions of COCO_PATH_MAX, coco_path_separator, HAVE_GFA and HAVE_STAT heavily used by functions in
 * coco_utilities.c */
#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW64__) || defined(__CYGWIN__)
#include <windows.h>
static const char *coco_path_separator = "\\";
#define COCO_PATH_MAX MAX_PATH
#define HAVE_GFA 1
#elif defined(__gnu_linux__)
#include <sys/stat.h>
#include <sys/types.h>
#include <linux/limits.h>
static const char *coco_path_separator = "/";
#define HAVE_STAT 1
#define COCO_PATH_MAX PATH_MAX
#elif defined(__APPLE__)
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/syslimits.h>
static const char *coco_path_separator = "/";
#define HAVE_STAT 1
#define COCO_PATH_MAX PATH_MAX
#elif defined(__FreeBSD__)
#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>
static const char *coco_path_separator = "/";
#define HAVE_STAT 1
#define COCO_PATH_MAX PATH_MAX
#elif (defined(__sun) || defined(sun)) && (defined(__SVR4) || defined(__svr4__))
/* Solaris */
#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>
static const char *coco_path_separator = "/";
#define HAVE_STAT 1
#define COCO_PATH_MAX PATH_MAX
#else
#error Unknown platform
#endif
#if !defined(COCO_PATH_MAX)
#error COCO_PATH_MAX undefined
#endif

/* Definitions needed for creating and removing directories */
/* Separately handle the special case of Microsoft Visual Studio 2008 with x86_64-w64-mingw32-gcc */
#if _MSC_VER
#include <direct.h>
#elif defined(__MINGW32__) || defined(__MINGW64__)
#include <dirent.h>
#else
#include <dirent.h>

#ifdef __cplusplus
extern "C" {
#endif

/* To silence the compiler (implicit-function-declaration warning). */
/** @cond */
int rmdir(const char *pathname);
int unlink(const char *file_name);
int mkdir(const char *pathname, mode_t mode);
/** @endcond */
#endif

/* Definition of the S_IRWXU constant needed to set file permissions */
#if defined(HAVE_GFA)
#define S_IRWXU 0700
#endif

/* To silence the Visual Studio compiler (C4996 warnings in the python build). */
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

#ifdef __cplusplus
}
#endif

#endif
#line 7 "code-experiments/src/coco_utilities.c"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>

#line 16 "code-experiments/src/coco_utilities.c"
#line 17 "code-experiments/src/coco_utilities.c"
#line 1 "code-experiments/src/coco_string.c"
/**
 * @file coco_string.c
 * @brief Definitions of functions that manipulate strings.
 */

#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#line 11 "code-experiments/src/coco_string.c"

static size_t *coco_allocate_vector_size_t(const size_t number_of_elements);

/**
 * @brief Creates a duplicate copy of string and returns a pointer to it.
 *
 * The caller is responsible for freeing the allocated memory using coco_free_memory().
 */
static char *coco_strdup(const char *string) {
  size_t len;
  char *duplicate;
  if (string == NULL)
    return NULL;
  len = strlen(string);
  duplicate = (char *) coco_allocate_memory(len + 1);
  memcpy(duplicate, string, len + 1);
  return duplicate;
}

/**
 * @brief The length of the buffer used in the coco_vstrdupf function.
 *
 * @note This should be handled differently!
 */
#define COCO_VSTRDUPF_BUFLEN 444

/**
 * @brief Formatted string duplication, with va_list arguments.
 */
static char *coco_vstrdupf(const char *str, va_list args) {
  static char buf[COCO_VSTRDUPF_BUFLEN];
  long written;
  /* apparently args can only be used once, therefore
   * len = vsnprintf(NULL, 0, str, args) to find out the
   * length does not work. Therefore we use a buffer
   * which limits the max length. Longer strings should
   * never appear anyway, so this is rather a non-issue. */

#if 0
  written = vsnprintf(buf, COCO_VSTRDUPF_BUFLEN - 2, str, args);
  if (written < 0)
  coco_error("coco_vstrdupf(): vsnprintf failed on '%s'", str);
#else /* less safe alternative, if vsnprintf is not available */
  assert(strlen(str) < COCO_VSTRDUPF_BUFLEN / 2 - 2);
  if (strlen(str) >= COCO_VSTRDUPF_BUFLEN / 2 - 2)
    coco_error("coco_vstrdupf(): string is too long");
  written = vsprintf(buf, str, args);
  if (written < 0)
    coco_error("coco_vstrdupf(): vsprintf failed on '%s'", str);
#endif
  if (written > COCO_VSTRDUPF_BUFLEN - 3)
    coco_error("coco_vstrdupf(): A suspiciously long string is tried to being duplicated '%s'", buf);
  return coco_strdup(buf);
}

#undef COCO_VSTRDUPF_BUFLEN

/**
 * Optional arguments are used like in sprintf.
 */
char *coco_strdupf(const char *str, ...) {
  va_list args;
  char *s;

  va_start(args, str);
  s = coco_vstrdupf(str, args);
  va_end(args);
  return s;
}

/**
 * @brief Returns a concatenate copy of string1 + string2.
 *
 * The caller is responsible for freeing the allocated memory using coco_free_memory().
 */
static char *coco_strconcat(const char *s1, const char *s2) {
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  char *s = (char *) coco_allocate_memory(len1 + len2 + 1);

  memcpy(s, s1, len1);
  memcpy(&s[len1], s2, len2 + 1);
  return s;
}

/**
 * @brief Returns the first index where seq occurs in base and -1 if it doesn't.
 *
 * @note If there is an equivalent standard C function, this can/should be removed.
 */
static long coco_strfind(const char *base, const char *seq) {
  const size_t strlen_seq = strlen(seq);
  const size_t last_first_idx = strlen(base) - strlen(seq);
  size_t i, j;

  if (strlen(base) < strlen(seq))
    return -1;

  for (i = 0; i <= last_first_idx; ++i) {
    if (base[i] == seq[0]) {
      for (j = 0; j < strlen_seq; ++j) {
        if (base[i + j] != seq[j])
          break;
      }
      if (j == strlen_seq) {
        if (i > 1e9)
          coco_error("coco_strfind(): strange values observed i=%lu, j=%lu, strlen(base)=%lu",
          		(unsigned long) i, (unsigned long) j, (unsigned long) strlen(base));
        return (long) i;
      }
    }
  }
  return -1;
}

/**
 * @brief Splits a string based on the given delimiter.
 *
 * Returns a pointer to the resulting substrings with NULL as the last one.
 * The caller is responsible for freeing the allocated memory using coco_free_memory().
 */
static char **coco_string_split(const char *string, const char delimiter) {

  char **result;
  char *str_copy, *ptr, *token;
  char str_delimiter[2];
  size_t i;
  size_t count = 1;

  str_copy = coco_strdup(string);

  /* Counts the parts between delimiters */
  ptr = str_copy;
  while (*ptr != '\0') {
    if (*ptr == delimiter) {
      count++;
    }
    ptr++;
  }
  /* Makes room for an empty string that will be appended at the end */
  count++;

  result = (char **) coco_allocate_memory(count * sizeof(char *));

  /* Iterates through tokens
   * NOTE: strtok() ignores multiple delimiters, therefore the final number of detected substrings might be
   * lower than the count. This is OK. */
  i = 0;
  /* A char* delimiter needs to be used, otherwise strtok() can surprise */
  str_delimiter[0] = delimiter;
  str_delimiter[1] = '\0';
  token = strtok(str_copy, str_delimiter);
  while (token)
  {
      assert(i < count);
      *(result + i++) = coco_strdup(token);
      token = strtok(NULL, str_delimiter);
  }
  *(result + i) = NULL;

  coco_free_memory(str_copy);

  return result;
}

/**
 * @brief Creates and returns a string with removed characters between from and to.
 *
 * If you wish to remove characters from the beginning of the string, set from to "".
 * If you wish to remove characters until the end of the string, set to to "".
 *
 * The caller is responsible for freeing the allocated memory using coco_free_memory().
 */
static char *coco_remove_from_string(const char *string, const char *from, const char *to) {

  char *result, *start, *stop;

  result = coco_strdup(string);

  if (strcmp(from, "") == 0) {
    /* Remove from the start */
    start = result;
  } else
    start = strstr(result, from);

  if (strcmp(to, "") == 0) {
    /* Remove until the end */
    stop = result + strlen(result);
  } else
    stop = strstr(result, to);

  if ((start == NULL) || (stop == NULL) || (stop < start)) {
    coco_error("coco_remove_from_string(): failed to remove characters between %s and %s from string %s",
        from, to, string);
    return NULL; /* Never reached */
  }

  memmove(start, stop, strlen(stop) + 1);

  return result;
}


/**
 * @brief Returns the numbers defined by the ranges.
 *
 * Reads ranges from a string of positive ranges separated by commas. For example: "-3,5-6,8-". Returns the
 * numbers that are defined by the ranges if min and max are used as their extremes. If the ranges with open
 * beginning/end are not allowed, use 0 as min/max. The returned string has an appended 0 to mark its end.
 * A maximum of max_count values is returned. If there is a problem with one of the ranges, the parsing stops
 * and the current result is returned. The memory of the returned object needs to be freed by the caller.
 */
static size_t *coco_string_parse_ranges(const char *string,
                                        const size_t min,
                                        const size_t max,
                                        const char *name,
                                        const size_t max_count) {

  char *ptr, *dash = NULL;
  char **ranges, **numbers;
  size_t i, j, count;
  size_t num[2];

  size_t *result;
  size_t i_result = 0;

  char *str = coco_strdup(string);

  /* Check for empty string */
  if ((str == NULL) || (strlen(str) == 0)) {
    coco_warning("coco_string_parse_ranges(): cannot parse empty ranges");
    coco_free_memory(str);
    return NULL;
  }

  ptr = str;
  /* Check for disallowed characters */
  while (*ptr != '\0') {
    if ((*ptr != '-') && (*ptr != ',') && !isdigit((unsigned char )*ptr)) {
      coco_warning("coco_string_parse_ranges(): problem parsing '%s' - cannot parse ranges with '%c'", str,
          *ptr);
      coco_free_memory(str);
      return NULL;
    } else
      ptr++;
  }

  /* Check for incorrect boundaries */
  if ((max > 0) && (min > max)) {
    coco_warning("coco_string_parse_ranges(): incorrect boundaries");
    coco_free_memory(str);
    return NULL;
  }

  result = coco_allocate_vector_size_t(max_count + 1);

  /* Split string to ranges w.r.t commas */
  ranges = coco_string_split(str, ',');
  coco_free_memory(str);

  if (ranges) {
    /* Go over the current range */
    for (i = 0; *(ranges + i); i++) {

      ptr = *(ranges + i);
      /* Count the number of '-' */
      count = 0;
      while (*ptr != '\0') {
        if (*ptr == '-') {
          if (count == 0)
            /* Remember the position of the first '-' */
            dash = ptr;
          count++;
        }
        ptr++;
      }
      /* Point again to the start of the range */
      ptr = *(ranges + i);

      /* Check for incorrect number of '-' */
      if (count > 1) {
        coco_warning("coco_string_parse_ranges(): problem parsing '%s' - too many '-'s", string);
        /* Cleanup */
        for (j = i; *(ranges + j); j++)
          coco_free_memory(*(ranges + j));
        coco_free_memory(ranges);
        if (i_result == 0) {
          coco_free_memory(result);
          return NULL;
        }
        result[i_result] = 0;
        return result;
      } else if (count == 0) {
        /* Range is in the format: n (no range) */
        num[0] = (size_t) strtol(ptr, NULL, 10);
        num[1] = num[0];
      } else {
        /* Range is in one of the following formats: n-m / -n / n- / - */

        /* Split current range to numbers w.r.t '-' */
        numbers = coco_string_split(ptr, '-');
        j = 0;
        if (numbers) {
          /* Read the numbers */
          for (j = 0; *(numbers + j); j++) {
            assert(j < 2);
            num[j] = (size_t) strtol(*(numbers + j), NULL, 10);
            coco_free_memory(*(numbers + j));
          }
        }
        coco_free_memory(numbers);

        if (j == 0) {
          /* Range is in the format - (open ends) */
          if ((min == 0) || (max == 0)) {
            coco_warning("coco_string_parse_ranges(): '%s' ranges cannot have an open ends; some ranges ignored", name);
            /* Cleanup */
            for (j = i; *(ranges + j); j++)
              coco_free_memory(*(ranges + j));
            coco_free_memory(ranges);
            if (i_result == 0) {
              coco_free_memory(result);
              return NULL;
            }
            result[i_result] = 0;
            return result;
          }
          num[0] = min;
          num[1] = max;
        } else if (j == 1) {
          if (dash - *(ranges + i) == 0) {
            /* Range is in the format -n */
            if (min == 0) {
              coco_warning("coco_string_parse_ranges(): '%s' ranges cannot have an open beginning; some ranges ignored", name);
              /* Cleanup */
              for (j = i; *(ranges + j); j++)
                coco_free_memory(*(ranges + j));
              coco_free_memory(ranges);
              if (i_result == 0) {
                coco_free_memory(result);
                return NULL;
              }
              result[i_result] = 0;
              return result;
            }
            num[1] = num[0];
            num[0] = min;
          } else {
            /* Range is in the format n- */
            if (max == 0) {
              coco_warning("coco_string_parse_ranges(): '%s' ranges cannot have an open end; some ranges ignored", name);
              /* Cleanup */
              for (j = i; *(ranges + j); j++)
                coco_free_memory(*(ranges + j));
              coco_free_memory(ranges);
              if (i_result == 0) {
                coco_free_memory(result);
                return NULL;
              }
              result[i_result] = 0;
              return result;
            }
            num[1] = max;
          }
        }
        /* if (j == 2), range is in the format n-m and there is nothing to do */
      }

      /* Make sure the boundaries are taken into account */
      if ((min > 0) && (num[0] < min)) {
        num[0] = min;
        coco_warning("coco_string_parse_ranges(): '%s' ranges adjusted to be >= %lu", name,
        		(unsigned long) min);
      }
      if ((max > 0) && (num[1] > max)) {
        num[1] = max;
        coco_warning("coco_string_parse_ranges(): '%s' ranges adjusted to be <= %lu", name,
        		(unsigned long) max);
      }
      if (num[0] > num[1]) {
        coco_warning("coco_string_parse_ranges(): '%s' ranges not within boundaries; some ranges ignored", name);
        /* Cleanup */
        for (j = i; *(ranges + j); j++)
          coco_free_memory(*(ranges + j));
        coco_free_memory(ranges);
        if (i_result == 0) {
          coco_free_memory(result);
          return NULL;
        }
        result[i_result] = 0;
        return result;
      }

      /* Write in result */
      for (j = num[0]; j <= num[1]; j++) {
        if (i_result > max_count - 1)
          break;
        result[i_result++] = j;
      }

      coco_free_memory(*(ranges + i));
      *(ranges + i) = NULL;
    }
  }

  coco_free_memory(ranges);

  if (i_result == 0) {
    coco_free_memory(result);
    return NULL;
  }

  result[i_result] = 0;
  return result;
}

/**
 * @brief Trims the given string (removes any leading and trailing spaces).
 *
 * If the string contains any leading spaces, the contents are shifted so that if it was dynamically
 * allocated, it can be still freed on the returned pointer.
 */
static char *coco_string_trim(char *string) {
	size_t len = 0;
	char *frontp = string;
	char *endp = NULL;

	if (string == NULL) {
		return NULL;
	}
	if (string[0] == '\0') {
		return string;
	}

	len = strlen(string);
	endp = string + len;

	/* Move the front and back pointers to address the first non-whitespace characters from each end. */
	while (isspace((unsigned char) *frontp)) {
		++frontp;
	}
	if (endp != frontp) {
		while (isspace((unsigned char) *(--endp)) && endp != frontp) {
		}
	}

	if (string + len - 1 != endp)
		*(endp + 1) = '\0';
	else if (frontp != string && endp == frontp)
		*string = '\0';

	/* Shift the string. Note the reuse of endp to mean the front of the string buffer now. */
	endp = string;
	if (frontp != string) {
		while (*frontp) {
			*endp++ = *frontp++;
		}
		*endp = '\0';
	}

	return string;
}
#line 18 "code-experiments/src/coco_utilities.c"

/***********************************************************************************************************/

/**
 * @brief Initializes the logging level to COCO_INFO.
 */
static coco_log_level_type_e coco_log_level = COCO_INFO;

/**
 * @param log_level Denotes the level of information given to the user through the standard output and
 * error streams. Can take on the values:
 * - "error" (only error messages are output),
 * - "warning" (only error and warning messages are output),
 * - "info" (only error, warning and info messages are output) and
 * - "debug" (all messages are output).
 * - "" does not set a new value
 * The default value is info.
 *
 * @return The previous coco_log_level value as an immutable string.
 */
const char *coco_set_log_level(const char *log_level) {

  coco_log_level_type_e previous_log_level = coco_log_level;

  if (strcmp(log_level, "error") == 0)
    coco_log_level = COCO_ERROR;
  else if (strcmp(log_level, "warning") == 0)
    coco_log_level = COCO_WARNING;
  else if (strcmp(log_level, "info") == 0)
    coco_log_level = COCO_INFO;
  else if (strcmp(log_level, "debug") == 0)
    coco_log_level = COCO_DEBUG;
  else if (strcmp(log_level, "") == 0) {
    /* Do nothing */
  } else {
    coco_warning("coco_set_log_level(): unknown level %s", log_level);
  }

  if (previous_log_level == COCO_ERROR)
    return "error";
  else if (previous_log_level == COCO_WARNING)
    return "warning";
  else if (previous_log_level == COCO_INFO)
    return "info";
  else if (previous_log_level == COCO_DEBUG)
    return "debug";
  else {
    coco_error("coco_set_log_level(): unknown previous log level");
    return "";
  }
}

/***********************************************************************************************************/

/**
 * @name Methods regarding file, directory and path manipulations
 */
/**@{*/
/**
 * @brief Creates a platform-dependent path from the given strings.
 *
 * @note The last argument must be NULL.
 * @note The first parameter must be able to accommodate path_max_length characters and the length
 * of the joined path must not exceed path_max_length characters.
 * @note Should work cross-platform.
 *
 * Usage examples:
 * - coco_join_path(base_path, 100, folder1, folder2, folder3, NULL) creates base_path/folder1/folder2/folder3
 * - coco_join_path(base_path, 100, folder1, file_name, NULL) creates base_path/folder1/file_name
 * @param path The base path; it's also where the joined path is stored to.
 * @param path_max_length The maximum length of the path.
 * @param ... Additional strings, must end with NULL
 */
static void coco_join_path(char *path, const size_t path_max_length, ...) {
  const size_t path_separator_length = strlen(coco_path_separator);
  va_list args;
  char *path_component;
  size_t path_length = strlen(path);

  va_start(args, path_max_length);
  while (NULL != (path_component = va_arg(args, char *))) {
    size_t component_length = strlen(path_component);
    if (path_length + path_separator_length + component_length >= path_max_length) {
      coco_error("coco_join_path() failed because the ${path} is too short.");
      return; /* never reached */
    }
    /* Both should be safe because of the above check. */
    if (strlen(path) > 0)
      strncat(path, coco_path_separator, path_max_length - strlen(path) - 1);
    strncat(path, path_component, path_max_length - strlen(path) - 1);
  }
  va_end(args);
}

/**
 * @brief Checks if the given directory exists.
 *
 * @note Should work cross-platform.
 *
 * @param path The given path.
 *
 * @return 1 if the path exists and corresponds to a directory and 0 otherwise.
 */
static int coco_directory_exists(const char *path) {
  int res;
#if defined(HAVE_GFA)
  DWORD dwAttrib = GetFileAttributesA(path);
  res = (dwAttrib != INVALID_FILE_ATTRIBUTES && (dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
#elif defined(HAVE_STAT)
  struct stat buf;
  res = (!stat(path, &buf) && S_ISDIR(buf.st_mode));
#else
#error Ooops
#endif
  return res;
}

/**
 * @brief Checks if the given file exists.
 *
 * @note Should work cross-platform.
 *
 * @param path The given path.
 *
 * @return 1 if the path exists and corresponds to a file and 0 otherwise.
 */
static int coco_file_exists(const char *path) {
  int res;
#if defined(HAVE_GFA)
  DWORD dwAttrib = GetFileAttributesA(path);
  res = (dwAttrib != INVALID_FILE_ATTRIBUTES) && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY);
#elif defined(HAVE_STAT)
  struct stat buf;
  res = (!stat(path, &buf) && !S_ISDIR(buf.st_mode));
#else
#error Ooops
#endif
  return res;
}

/**
 * @brief Calls the right mkdir() method (depending on the platform).
 *
 * @param path The directory path.
 *
 * @return 0 on successful completion, and -1 on error.
 */
static int coco_mkdir(const char *path) {
#if _MSC_VER
  return _mkdir(path);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  return mkdir(path);
#else
  return mkdir(path, S_IRWXU);
#endif
}

/**
 * @brief Creates a directory with full privileges for the user.
 *
 * @note Should work cross-platform.
 *
 * @param path The directory path.
 */
static void coco_create_directory(const char *path) {
  char *tmp = NULL;
  char *p;
  size_t len = strlen(path);
  char path_sep = coco_path_separator[0];

  /* Nothing to do if the path exists. */
  if (coco_directory_exists(path))
    return;

  tmp = coco_strdup(path);
  /* Remove possible trailing slash */
  if (tmp[len - 1] == path_sep)
    tmp[len - 1] = 0;
  for (p = tmp + 1; *p; p++) {
    if (*p == path_sep) {
      *p = 0;
      if (!coco_directory_exists(tmp)) {
        if (0 != coco_mkdir(tmp))
          coco_error("coco_create_path(): failed creating %s", tmp);
      }
      *p = path_sep;
    }
  }
  if (0 != coco_mkdir(tmp))
    coco_error("coco_create_path(): failed creating %s", tmp);
  coco_free_memory(tmp);
  return;
}

/* Commented to silence the compiler (unused function warning) */
#if 0
/**
 * @brief Creates a unique file name from the given file_name.
 *
 * If the file_name does not yet exit, it is left as is, otherwise it is changed(!) by prepending a number
 * to it. If filename.ext already exists, 01-filename.ext will be tried. If this one exists as well,
 * 02-filename.ext will be tried, and so on. If 99-filename.ext exists as well, the function throws
 * an error.
 */
static void coco_create_unique_filename(char **file_name) {

  int counter = 1;
  char *new_file_name;

  /* Do not change the file_name if it does not yet exist */
  if (!coco_file_exists(*file_name)) {
    return;
  }

  while (counter < 99) {

    new_file_name = coco_strdupf("%02d-%s", counter, *file_name);

    if (!coco_file_exists(new_file_name)) {
      coco_free_memory(*file_name);
      *file_name = new_file_name;
      return;
    } else {
      counter++;
      coco_free_memory(new_file_name);
    }

  }

  coco_free_memory(new_file_name);
  coco_error("coco_create_unique_filename(): could not create a unique file name");
  return; /* Never reached */
}
#endif

/**
 * @brief Creates a unique directory from the given path.
 *
 * If the given path does not yet exit, it is left as is, otherwise it is changed(!) by appending a number
 * to it. If path already exists, path-01 will be tried. If this one exists as well, path-02 will be tried,
 * and so on. If path-99 exists as well, the function throws an error.
 */
static void coco_create_unique_directory(char **path) {

  int counter = 1;
  char *new_path;

  /* Create the path if it does not yet exist */
  if (!coco_directory_exists(*path)) {
    coco_create_directory(*path);
    return;
  }

  while (counter < 999) {

    new_path = coco_strdupf("%s-%03d", *path, counter);

    if (!coco_directory_exists(new_path)) {
      coco_free_memory(*path);
      *path = new_path;
      coco_create_directory(*path);
      return;
    } else {
      counter++;
      coco_free_memory(new_path);
    }

  }

  coco_error("coco_create_unique_path(): could not create a unique path with name %s", *path);
  return; /* Never reached */
}

/**
 * The method should work across different platforms/compilers.
 *
 * @path The path to the directory
 *
 * @return 0 on successful completion, and -1 on error.
 */
int coco_remove_directory(const char *path) {
#if _MSC_VER
  WIN32_FIND_DATA find_data_file;
  HANDLE find_handle = NULL;
  char *buf;
  int r = -1;
  int r2 = -1;

  buf = coco_strdupf("%s\\*.*", path);
  /* Nothing to do if the folder does not exist */
  if ((find_handle = FindFirstFile(buf, &find_data_file)) == INVALID_HANDLE_VALUE) {
    coco_free_memory(buf);
    return 0;
  }
  coco_free_memory(buf);

  do {
    r = 0;

    /* Skip the names "." and ".." as we don't want to recurse on them */
    if (strcmp(find_data_file.cFileName, ".") != 0 && strcmp(find_data_file.cFileName, "..") != 0) {
      /* Build the new path using the argument path the file/folder name we just found */
      buf = coco_strdupf("%s\\%s", path, find_data_file.cFileName);

      if (find_data_file.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
        /* Buf is a directory, recurse on it */
        r2 = coco_remove_directory(buf);
      } else {
        /* Buf is a file, delete it */
        /* Careful, DeleteFile returns 0 if it fails and nonzero otherwise! */
        r2 = -(DeleteFile(buf) == 0);
      }

      coco_free_memory(buf);
    }

    r = r2;

  }while (FindNextFile(find_handle, &find_data_file)); /* Find the next file */

  FindClose(find_handle);

  if (!r) {
    /* Path is an empty directory, delete it */
    /* Careful, RemoveDirectory returns 0 if it fails and nonzero otherwise! */
    r = -(RemoveDirectory(path) == 0);
  }

  return r;
#else
  DIR *d = opendir(path);
  int r = -1;
  int r2 = -1;
  char *buf;

  /* Nothing to do if the folder does not exist */
  if (!coco_directory_exists(path))
    return 0;

  if (d) {
    struct dirent *p;

    r = 0;

    while (!r && (p = readdir(d))) {

      /* Skip the names "." and ".." as we don't want to recurse on them */
      if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, "..")) {
        continue;
      }

      buf = coco_strdupf("%s/%s", path, p->d_name);
      if (buf) {
        if (coco_directory_exists(buf)) {
          /* Buf is a directory, recurse on it */
          r2 = coco_remove_directory(buf);
        } else {
          /* Buf is a file, delete it */
          r2 = unlink(buf);
        }
      }
      coco_free_memory(buf);

      r = r2;
    }

    closedir(d);
  }

  if (!r) {
    /* Path is an empty directory, delete it */
    r = rmdir(path);
  }

  return r;
#endif
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding memory allocations
 */
/**@{*/
double *coco_allocate_vector(const size_t number_of_elements) {
  const size_t block_size = number_of_elements * sizeof(double);
  return (double *) coco_allocate_memory(block_size);
}

/**
 * @brief Allocates memory for a vector and sets all its elements to value.
 */
static double *coco_allocate_vector_with_value(const size_t number_of_elements, double value) {
  const size_t block_size = number_of_elements * sizeof(double);
  double *vector = (double *) coco_allocate_memory(block_size);
  size_t i;

  for (i = 0; i < number_of_elements; i++)
  	vector[i] = value;

  return vector;
}

/**
 * @brief Safe memory allocation for a vector with size_t elements that either succeeds or triggers a
 * coco_error.
 */
static size_t *coco_allocate_vector_size_t(const size_t number_of_elements) {
  const size_t block_size = number_of_elements * sizeof(size_t);
  return (size_t *) coco_allocate_memory(block_size);
}

static char *coco_allocate_string(const size_t number_of_elements) {
  const size_t block_size = number_of_elements * sizeof(char);
  return (char *) coco_allocate_memory(block_size);
}

static double *coco_duplicate_vector(const double *src, const size_t number_of_elements) {
  size_t i;
  double *dst;

  assert(src != NULL);
  assert(number_of_elements > 0);

  dst = coco_allocate_vector(number_of_elements);
  for (i = 0; i < number_of_elements; ++i) {
    dst[i] = src[i];
  }
  return dst;
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding string options
 */
/**@{*/

/**
 * @brief Allocates an option keys structure holding the given number of option keys.
 */
static coco_option_keys_t *coco_option_keys_allocate(const size_t count, const char **keys) {

  size_t i;
  coco_option_keys_t *option_keys;

  if ((count == 0) || (keys == NULL))
    return NULL;

  option_keys = (coco_option_keys_t *) coco_allocate_memory(sizeof(*option_keys));

  option_keys->keys = (char **) coco_allocate_memory(count * sizeof(char *));
  for (i = 0; i < count; i++) {
    assert(keys[i]);
    option_keys->keys[i] = coco_strdup(keys[i]);
  }
  option_keys->count = count;

  return option_keys;
}

/**
 * @brief Frees the given option keys structure.
 */
static void coco_option_keys_free(coco_option_keys_t *option_keys) {

  size_t i;

  if (option_keys) {
    for (i = 0; i < option_keys->count; i++) {
      coco_free_memory(option_keys->keys[i]);
    }
    coco_free_memory(option_keys->keys);
    coco_free_memory(option_keys);
  }
}

/**
 * @brief Returns redundant option keys (the ones present in given_option_keys but not in known_option_keys).
 */
static coco_option_keys_t *coco_option_keys_get_redundant(const coco_option_keys_t *known_option_keys,
                                                          const coco_option_keys_t *given_option_keys) {

  size_t i, j, count = 0;
  int found;
  char **redundant_keys;
  coco_option_keys_t *redundant_option_keys;

  assert(known_option_keys != NULL);
  assert(given_option_keys != NULL);

  /* Find the redundant keys */
  redundant_keys = (char **) coco_allocate_memory(given_option_keys->count * sizeof(char *));
  for (i = 0; i < given_option_keys->count; i++) {
    found = 0;
    for (j = 0; j < known_option_keys->count; j++) {
      if (strcmp(given_option_keys->keys[i], known_option_keys->keys[j]) == 0) {
        found = 1;
        break;
      }
    }
    if (!found) {
      redundant_keys[count++] = coco_strdup(given_option_keys->keys[i]);
    }
  }
  redundant_option_keys = coco_option_keys_allocate(count, (const char**) redundant_keys);

  /* Free memory */
  for (i = 0; i < count; i++) {
    coco_free_memory(redundant_keys[i]);
  }
  coco_free_memory(redundant_keys);

  return redundant_option_keys;
}

/**
 * @brief Adds additional option keys to the given basic option keys (changes the basic keys).
 */
static void coco_option_keys_add(coco_option_keys_t **basic_option_keys,
                                 const coco_option_keys_t *additional_option_keys) {

  size_t i, j;
  size_t new_count;
  char **new_keys;
  coco_option_keys_t *new_option_keys;

  assert(*basic_option_keys != NULL);
  if (additional_option_keys == NULL)
    return;

  /* Construct the union of both keys */
  new_count = (*basic_option_keys)->count + additional_option_keys->count;
  new_keys = (char **) coco_allocate_memory(new_count * sizeof(char *));
  for (i = 0; i < (*basic_option_keys)->count; i++) {
    new_keys[i] = coco_strdup((*basic_option_keys)->keys[i]);
  }
  for (j = 0; j < additional_option_keys->count; j++) {
    new_keys[(*basic_option_keys)->count + j] = coco_strdup(additional_option_keys->keys[j]);
  }
  new_option_keys = coco_option_keys_allocate(new_count, (const char**) new_keys);

  /* Free the old basic keys */
  coco_option_keys_free(*basic_option_keys);
  *basic_option_keys = new_option_keys;
  for (i = 0; i < new_count; i++) {
    coco_free_memory(new_keys[i]);
  }
  coco_free_memory(new_keys);
}

/**
 * @brief Creates an instance of option keys from the given string of options containing keys and values
 * separated by colons.
 *
 * @note Relies heavily on the "key: value" format and might fail if the number of colons doesn't match the
 * number of keys.
 */
static coco_option_keys_t *coco_option_keys(const char *option_string) {

  size_t i;
  char **keys;
  coco_option_keys_t *option_keys = NULL;
  char *string_to_parse, *key;

  /* Check for empty string */
  if ((option_string == NULL) || (strlen(option_string) == 0)) {
	    return NULL;
  }

  /* Split the options w.r.t ':' */
  keys = coco_string_split(option_string, ':');

  if (keys) {
    /* Keys now contain something like this: "values_of_previous_key this_key" except for the first, which
     * contains only the key and the last, which contains only the previous values */
    for (i = 0; *(keys + i); i++) {
      string_to_parse = coco_strdup(*(keys + i));

      /* Remove any leading and trailing spaces */
      string_to_parse = coco_string_trim(string_to_parse);

      /* Stop if this is the last substring (contains a value and no key) */
      if ((i > 0) && (*(keys + i + 1) == NULL)) {
        coco_free_memory(string_to_parse);
        break;
      }

      /* Disregard everything before the last space */
      key = strrchr(string_to_parse, ' ');
      if ((key == NULL) || (i == 0)) {
        /* No spaces left (or this is the first key), everything is the key */
        key = string_to_parse;
      } else {
        /* Move to the start of the key (one char after the space) */
        key++;
      }

      /* Put the key in keys */
      coco_free_memory(*(keys + i));
      *(keys + i) = coco_strdup(key);
      coco_free_memory(string_to_parse);
    }

    option_keys = coco_option_keys_allocate(i, (const char**) keys);

    /* Free the keys */
    for (i = 0; *(keys + i); i++) {
      coco_free_memory(*(keys + i));
    }
    coco_free_memory(keys);
  }

  return option_keys;
}

/**
 * @brief Creates and returns a string containing the info_string and all keys from option_keys.
 *
 * Can be used to output information about the given option_keys.
 */
static char *coco_option_keys_get_output_string(const coco_option_keys_t *option_keys,
                                                const char *info_string) {
  size_t i;
  char *string = NULL, *new_string;

  if ((option_keys != NULL) && (option_keys->count > 0)) {

    string = coco_strdup(info_string);
    for (i = 0; i < option_keys->count; i++) {
      new_string = coco_strdupf("%s %s\n", string, option_keys->keys[i]);
      coco_free_memory(string);
      string = new_string;
    }
  }

  return string;
}

/**
 * @brief Parses options in the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - value needs to be a single string (no spaces allowed)
 *
 * @return The number of successful assignments.
 */
static int coco_options_read(const char *options, const char *name, const char *format, void *pointer) {

  long i1, i2;

  if ((!options) || (strlen(options) == 0))
    return 0;

  i1 = coco_strfind(options, name);
  if (i1 < 0)
    return 0;
  i2 = i1 + coco_strfind(&options[i1], ":") + 1;

  /* Remove trailing whitespaces */
  while (isspace((unsigned char) options[i2]))
    i2++;

  if (i2 <= i1){
    return 0;
  }

  return sscanf(&options[i2], format, pointer);
}

/**
 * @brief Reads an integer from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be an integer
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_int(const char *options, const char *name, int *pointer) {
  return coco_options_read(options, name, " %i", pointer);
}

/**
 * @brief Reads a size_t from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be a size_t
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_size_t(const char *options, const char *name, size_t *pointer) {
  return coco_options_read(options, name, "%lu", pointer);
}

/**
 * @brief Reads a double value from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be a double
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_double(const char *options, const char *name, double *pointer) {
  return coco_options_read(options, name, "%lf", pointer);
}

/**
 * @brief Reads a string from options using the form "name1: value1 name2: value2".
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 * - the value corresponding to the given name needs to be a string - either a single word or multiple words
 * in double quotes
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_string(const char *options, const char *name, char *pointer) {

  long i1, i2;

  if ((!options) || (strlen(options) == 0))
    return 0;

  i1 = coco_strfind(options, name);
  if (i1 < 0)
    return 0;
  i2 = i1 + coco_strfind(&options[i1], ":") + 1;

  /* Remove trailing white spaces */
  while (isspace((unsigned char) options[i2]))
    i2++;

  if (i2 <= i1){
    return 0;
  }

  if (options[i2] == '\"') {
    /* The value starts with a quote: read everything between two quotes into a string */
    return sscanf(&options[i2], "\"%[^\"]\"", pointer);
  } else
    return sscanf(&options[i2], "%s", pointer);
}

/**
 * @brief Reads (possibly delimited) values from options using the form "name1: value1,value2,value3 name2: value4",
 * i.e. reads all characters from the corresponding name up to the next whitespace or end of string.
 *
 * Formatting requirements:
 * - name and value need to be separated by a colon (spaces are optional)
 *
 * @return The number of successful assignments.
 */
static int coco_options_read_values(const char *options, const char *name, char *pointer) {

  long i1, i2;
  int i;

  if ((!options) || (strlen(options) == 0))
    return 0;

  i1 = coco_strfind(options, name);
  if (i1 < 0)
    return 0;
  i2 = i1 + coco_strfind(&options[i1], ":") + 1;

  /* Remove trailing white spaces */
  while (isspace((unsigned char) options[i2]))
    i2++;

  if (i2 <= i1) {
    return 0;
  }

  i = 0;
  while (!isspace((unsigned char) options[i2 + i]) && (options[i2 + i] != '\0')) {
    pointer[i] = options[i2 + i];
    i++;
  }
  pointer[i] = '\0';
  return i;
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods implementing functions on double values not contained in C89 standard
 */
/**@{*/

/**
 * @brief Rounds the given double to the nearest integer.
 */
static double coco_double_round(const double number) {
  return floor(number + 0.5);
}

/**
 * @brief Returns the maximum of a and b.
 */
static double coco_double_max(const double a, const double b) {
  if (a >= b) {
    return a;
  } else {
    return b;
  }
}

/**
 * @brief Returns the minimum of a and b.
 */
static double coco_double_min(const double a, const double b) {
  if (a <= b) {
    return a;
  } else {
    return b;
  }
}

/**
 * @brief Performs a "safer" double to size_t conversion.
 */
static size_t coco_double_to_size_t(const double number) {
  return (size_t) coco_double_round(number);
}

/**
 * @brief  Returns 1 if |a - b| < precision and 0 otherwise.
 */
static int coco_double_almost_equal(const double a, const double b, const double precision) {
  return (fabs(a - b) < precision);
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods handling NAN and INFINITY
 */
/**@{*/

/**
 * @brief Returns 1 if x is NAN and 0 otherwise.
 */
static int coco_is_nan(const double x) {
  return (isnan(x) || (x != x) || !(x == x) || ((x >= NAN / (1 + 1e-9)) && (x <= NAN * (1 + 1e-9))));
}

/**
 * @brief Returns 1 if the input vector of dimension dim contains any NAN values and 0 otherwise.
 */
static int coco_vector_contains_nan(const double *x, const size_t dim) {
	size_t i;
	for (i = 0; i < dim; i++) {
		if (coco_is_nan(x[i]))
		  return 1;
	}
	return 0;
}

/**
 * @brief Sets all dim values of y to NAN.
 */
static void coco_vector_set_to_nan(double *y, const size_t dim) {
	size_t i;
	for (i = 0; i < dim; i++) {
		y[i] = NAN;
	}
}

/**
 * @brief Returns 1 if x is INFINITY and 0 otherwise.
 */
static int coco_is_inf(const double x) {
	if (coco_is_nan(x))
		return 0;
	return (isinf(x) || (x <= -INFINITY) || (x >= INFINITY));
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Miscellaneous methods
 */
/**@{*/

/**
 * @brief Returns the current time as a string.
 *
 * The caller is responsible for freeing the allocated memory using coco_free_memory().
 */
static char *coco_current_time_get_string(void) {
  time_t timer;
  char *time_string = coco_allocate_string(30);
  struct tm* tm_info;
  time(&timer);
  tm_info = localtime(&timer);
  assert(tm_info != NULL);
  strftime(time_string, 30, "%d.%m.%y %H:%M:%S", tm_info);
  return time_string;
}

/**
 * @brief Returns the number of positive numbers pointed to by numbers (the count stops when the first
 * 0 is encountered of max_count numbers have been read).
 *
 * If there are more than max_count numbers, a coco_error is raised. The name argument is used
 * only to provide more informative output in case of any problems.
 */
static size_t coco_count_numbers(const size_t *numbers, const size_t max_count, const char *name) {

  size_t count = 0;
  while ((count < max_count) && (numbers[count] != 0)) {
    count++;
  }
  if (count == max_count) {
    coco_error("coco_count_numbers(): over %lu numbers in %s", (unsigned long) max_count, name);
    return 0; /* Never reached*/
  }

  return count;
}

/**@}*/

/***********************************************************************************************************/
#line 17 "code-experiments/src/coco_suite.c"

#line 1 "code-experiments/src/suite_bbob.c"
/**
 * @file suite_bbob.c
 * @brief Implementation of the bbob suite containing 24 noiseless single-objective functions in 6
 * dimensions.
 */

#line 8 "code-experiments/src/suite_bbob.c"

#line 1 "code-experiments/src/f_attractive_sector.c"
/**
 * @file f_attractive_sector.c
 * @brief Implementation of the attractive sector function and problem.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/coco_problem.c"
/**
 * @file coco_problem.c
 * @brief Definitions of functions regarding COCO problems.
 */

#include <float.h>
#line 8 "code-experiments/src/coco_problem.c"
#line 9 "code-experiments/src/coco_problem.c"

#line 11 "code-experiments/src/coco_problem.c"

/***********************************************************************************************************/

/**
 * @name Methods regarding the basic COCO problem
 */
/**@{*/
/**
 * Evaluates the problem function, increases the number of evaluations and updates the best observed value
 * and the best observed evaluation number.
 *
 * @note Both x and y must point to correctly sized allocated memory regions.
 *
 * @param problem The given COCO problem.
 * @param x The decision vector.
 * @param y The objective vector that is the result of the evaluation (in single-objective problems only the
 * first vector item is being set).
 */
void coco_evaluate_function(coco_problem_t *problem, const double *x, double *y) {
  /* implements a safer version of problem->evaluate(problem, x, y) */
	size_t i, j;
	assert(problem != NULL);
  assert(problem->evaluate_function != NULL);

  /* Set objective vector to INFINITY if the decision vector contains any INFINITY values */
	for (i = 0; i < coco_problem_get_dimension(problem); i++) {
		if (coco_is_inf(x[i])) {
			for (j = 0; j < coco_problem_get_number_of_objectives(problem); j++) {
				y[j] = fabs(x[i]);
			}
	  	return;
		}
  }

  /* Set objective vector to NAN if the decision vector contains any NAN values */
  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  problem->evaluate_function(problem, x, y);
  problem->evaluations++; /* each derived class has its own counter, only the most outer will be visible */

  /* A little bit of bookkeeping */
  if (y[0] < problem->best_observed_fvalue[0]) {
    problem->best_observed_fvalue[0] = y[0];
    problem->best_observed_evaluation[0] = problem->evaluations;
  }

}

/**
 * @note None of the problems implement this function yet!
 * @note Both x and y must point to correctly sized allocated memory regions.
 *
 * @param problem The given COCO problem.
 * @param x The decision vector.
 * @param y The vector of constraints that is the result of the evaluation.
 */
void coco_evaluate_constraint(coco_problem_t *problem, const double *x, double *y) {
  /* implements a safer version of problem->evaluate(problem, x, y) */
  assert(problem != NULL);
  if (problem->evaluate_constraint == NULL) {
    coco_error("coco_evaluate_constraint(): No constraint function implemented for problem %s",
        problem->problem_id);
  }
  problem->evaluate_constraint(problem, x, y);
}

/**
 * Evaluates and logs the given solution (as the coco_evaluate_function), but does not return the evaluated
 * value.
 *
 * @note None of the observers implements this function yet!
 * @note x must point to a correctly sized allocated memory region.

 * @param problem The given COCO problem.
 * @param x The decision vector.
 */
void coco_recommend_solution(coco_problem_t *problem, const double *x) {
  assert(problem != NULL);
  if (problem->recommend_solution == NULL) {
    coco_error("coco_recommend_solutions(): No recommend solution function implemented for problem %s",
        problem->problem_id);
  }
  problem->recommend_solution(problem, x);
}

/***********************************************************************************************************/

/**
 * @brief Allocates a new coco_problem_t for the given number of variables, number of objectives and
 * number of constraints.
 */
static coco_problem_t *coco_problem_allocate(const size_t number_of_variables,
                                             const size_t number_of_objectives,
                                             const size_t number_of_constraints) {
  coco_problem_t *problem;
  problem = (coco_problem_t *) coco_allocate_memory(sizeof(*problem));
  /* Initialize fields to sane/safe defaults */
  problem->initial_solution = NULL;
  problem->evaluate_function = NULL;
  problem->evaluate_constraint = NULL;
  problem->recommend_solution = NULL;
  problem->problem_free_function = NULL;
  problem->number_of_variables = number_of_variables;
  problem->number_of_objectives = number_of_objectives;
  problem->number_of_constraints = number_of_constraints;
  problem->smallest_values_of_interest = coco_allocate_vector(number_of_variables);
  problem->largest_values_of_interest = coco_allocate_vector(number_of_variables);
  problem->best_parameter = coco_allocate_vector(number_of_variables);
  problem->best_value = coco_allocate_vector(number_of_objectives);
  if (number_of_objectives > 1)
    problem->nadir_value = coco_allocate_vector(number_of_objectives);
  else
    problem->nadir_value = NULL;
  problem->problem_name = NULL;
  problem->problem_id = NULL;
  problem->problem_type = NULL;
  problem->evaluations = 0;
  problem->final_target_delta[0] = 1e-8; /* in case to be modified by the benchmark */
  problem->best_observed_fvalue[0] = DBL_MAX;
  problem->best_observed_evaluation[0] = 0;
  problem->suite = NULL; /* To be initialized in the coco_suite_get_problem_from_indices() function */
  problem->suite_dep_index = 0;
  problem->suite_dep_function = 0;
  problem->suite_dep_instance = 0;
  problem->data = NULL;
  return problem;
}

/**
 * @brief Creates a duplicate of the 'other' problem for all fields except for data, which points to NULL.
 */
static coco_problem_t *coco_problem_duplicate(const coco_problem_t *other) {
  size_t i;
  coco_problem_t *problem;
  problem = coco_problem_allocate(other->number_of_variables, other->number_of_objectives,
      other->number_of_constraints);

  problem->initial_solution = other->initial_solution;
  problem->evaluate_function = other->evaluate_function;
  problem->evaluate_constraint = other->evaluate_constraint;
  problem->recommend_solution = other->recommend_solution;
  problem->problem_free_function = other->problem_free_function;

  for (i = 0; i < problem->number_of_variables; ++i) {
    problem->smallest_values_of_interest[i] = other->smallest_values_of_interest[i];
    problem->largest_values_of_interest[i] = other->largest_values_of_interest[i];
    if (other->best_parameter)
      problem->best_parameter[i] = other->best_parameter[i];
  }

  if (other->best_value)
    for (i = 0; i < problem->number_of_objectives; ++i) {
      problem->best_value[i] = other->best_value[i];
    }

  if (other->nadir_value)
    for (i = 0; i < problem->number_of_objectives; ++i) {
      problem->nadir_value[i] = other->nadir_value[i];
    }

  problem->problem_name = coco_strdup(other->problem_name);
  problem->problem_id = coco_strdup(other->problem_id);
  problem->problem_type = coco_strdup(other->problem_type);

  problem->evaluations = other->evaluations;
  problem->final_target_delta[0] = other->final_target_delta[0];
  problem->best_observed_fvalue[0] = other->best_observed_fvalue[0];
  problem->best_observed_evaluation[0] = other->best_observed_evaluation[0];

  problem->suite = other->suite;
  problem->suite_dep_index = other->suite_dep_index;
  problem->suite_dep_function = other->suite_dep_function;
  problem->suite_dep_instance = other->suite_dep_instance;

  problem->data = NULL;

  return problem;
}

/**
 * @brief Allocates a problem using scalar values for smallest_value_of_interest, largest_value_of_interest
 * and best_parameter.
 */
static coco_problem_t *coco_problem_allocate_from_scalars(const char *problem_name,
                                                          coco_evaluate_function_t evaluate_function,
                                                          coco_problem_free_function_t problem_free_function,
                                                          const size_t number_of_variables,
                                                          const double smallest_value_of_interest,
                                                          const double largest_value_of_interest,
                                                          const double best_parameter) {
  size_t i;
  coco_problem_t *problem = coco_problem_allocate(number_of_variables, 1, 0);

  problem->problem_name = coco_strdup(problem_name);
  problem->number_of_variables = number_of_variables;
  problem->number_of_objectives = 1;
  problem->number_of_constraints = 0;
  problem->evaluate_function = evaluate_function;
  problem->problem_free_function = problem_free_function;

  for (i = 0; i < number_of_variables; ++i) {
    problem->smallest_values_of_interest[i] = smallest_value_of_interest;
    problem->largest_values_of_interest[i] = largest_value_of_interest;
    problem->best_parameter[i] = best_parameter;
  }
  return problem;
}

void coco_problem_free(coco_problem_t *problem) {
  assert(problem != NULL);
  if (problem->problem_free_function != NULL) {
    problem->problem_free_function(problem);
  } else {
    /* Best guess at freeing all relevant structures */
    if (problem->smallest_values_of_interest != NULL)
      coco_free_memory(problem->smallest_values_of_interest);
    if (problem->largest_values_of_interest != NULL)
      coco_free_memory(problem->largest_values_of_interest);
    if (problem->best_parameter != NULL)
      coco_free_memory(problem->best_parameter);
    if (problem->best_value != NULL)
      coco_free_memory(problem->best_value);
    if (problem->nadir_value != NULL)
      coco_free_memory(problem->nadir_value);
    if (problem->problem_name != NULL)
      coco_free_memory(problem->problem_name);
    if (problem->problem_id != NULL)
      coco_free_memory(problem->problem_id);
    if (problem->problem_type != NULL)
      coco_free_memory(problem->problem_type);
    if (problem->data != NULL)
      coco_free_memory(problem->data);
    problem->smallest_values_of_interest = NULL;
    problem->largest_values_of_interest = NULL;
    problem->best_parameter = NULL;
    problem->best_value = NULL;
    problem->nadir_value = NULL;
    problem->suite = NULL;
    problem->data = NULL;
    coco_free_memory(problem);
  }
}

/***********************************************************************************************************/

/**
 * @brief Checks whether the given string is in the right format to be a problem_id.
 *
 * No non-alphanumeric characters besides '-', '_' and '.' are allowed.
 */
static int coco_problem_id_is_fine(const char *id, ...) {
  va_list args;
  const int reject = 0;
  const int accept = 1;
  const char *cp;
  char *s;
  int result = accept;

  va_start(args, id);
  s = coco_vstrdupf(id, args);
  va_end(args);
  for (cp = s; *cp != '\0'; ++cp) {
    if (('A' <= *cp) && (*cp <= 'Z'))
      continue;
    if (('a' <= *cp) && (*cp <= 'z'))
      continue;
    if ((*cp == '_') || (*cp == '-'))
      continue;
    if (('0' <= *cp) && (*cp <= '9'))
      continue;
    result = reject;
  }
  coco_free_memory(s);
  return result;
}

/**
 * @brief Sets the problem_id using formatted printing (as in printf).
 *
 * Takes care of memory (de-)allocation and verifies that the problem_id is in the correct format.
 */
static void coco_problem_set_id(coco_problem_t *problem, const char *id, ...) {
  va_list args;

  va_start(args, id);
  if (problem->problem_id != NULL)
    coco_free_memory(problem->problem_id);
  problem->problem_id = coco_vstrdupf(id, args);
  va_end(args);
  if (!coco_problem_id_is_fine(problem->problem_id)) {
    coco_error("Problem id should only contain standard chars, not like '%s'", problem->problem_id);
  }
}

/**
 * @brief Sets the problem_name using formatted printing (as in printf).
 *
 * Takes care of memory (de-)allocation.
 */
static void coco_problem_set_name(coco_problem_t *problem, const char *name, ...) {
  va_list args;

  va_start(args, name);
  if (problem->problem_name != NULL)
    coco_free_memory(problem->problem_name);
  problem->problem_name = coco_vstrdupf(name, args);
  va_end(args);
}

/**
 * @brief Sets the problem_type using formatted printing (as in printf).
 *
 * Takes care of memory (de-)allocation.
 */
static void coco_problem_set_type(coco_problem_t *problem, const char *type, ...) {
  va_list args;

  va_start(args, type);
  if (problem->problem_type != NULL)
    coco_free_memory(problem->problem_type);
  problem->problem_type = coco_vstrdupf(type, args);
  va_end(args);
}

size_t coco_problem_get_evaluations(const coco_problem_t *problem) {
  assert(problem != NULL);
  return problem->evaluations;
}

/**
 * @brief Returns 1 if the best parameter is not (close to) zero and 0 otherwise.
 */
static int coco_problem_best_parameter_not_zero(const coco_problem_t *problem) {
	size_t i = 0;
	int best_is_zero = 1;

	if (coco_vector_contains_nan(problem->best_parameter, problem->number_of_variables))
		return 1;

	while (i < problem->number_of_variables && best_is_zero) {
	      best_is_zero = coco_double_almost_equal(problem->best_parameter[i], 0, 1e-9);
	      i++;
	  }

	return !best_is_zero;
}

/**
 * @note Can be used to prevent unnecessary burning of CPU time.
 */
int coco_problem_final_target_hit(const coco_problem_t *problem) {
  assert(problem != NULL);
  if (coco_problem_get_number_of_objectives(problem) != 1 ||
      coco_problem_get_evaluations(problem) < 1) 
    return 0;
  if (problem->best_value == NULL)
    return 0;
  return problem->best_observed_fvalue[0] <= problem->best_value[0] + problem->final_target_delta[0] ?
    1 : 0;
}

/**
 * @note Tentative...
 */
double coco_problem_get_best_observed_fvalue1(const coco_problem_t *problem) {
  assert(problem != NULL);
  return problem->best_observed_fvalue[0];
}

/**
 * @note This function breaks the black-box property: the returned  value is not
 * meant to be used by the optimization algorithm.
 */
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->best_value != NULL);
  assert(problem->final_target_delta != NULL);
  return problem->best_value[0] + problem->final_target_delta[0];
}

/**
 * @note Do not modify the returned string! If you free the problem, the returned pointer becomes invalid.
 * When in doubt, use coco_strdup() on the returned value.
 */
const char *coco_problem_get_name(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->problem_name != NULL);
  return problem->problem_name;
}

/**
 * The ID is guaranteed to contain only characters in the set [a-z0-9_-]. It should therefore be safe to use
 * it to construct filenames or other identifiers.
 *
 * Each problem ID should be unique within each benchmark suite.
 *
 * @note Do not modify the returned string! If you free the problem, the returned pointer becomes invalid.
 * When in doubt, use coco_strdup() on the returned value.
 */
const char *coco_problem_get_id(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->problem_id != NULL);
  return problem->problem_id;
}

const char *coco_problem_get_type(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->problem_type != NULL);
  return problem->problem_type;
}

size_t coco_problem_get_dimension(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->number_of_variables > 0);
  return problem->number_of_variables;
}

size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->number_of_objectives > 0);
  return problem->number_of_objectives;
}

size_t coco_problem_get_number_of_constraints(const coco_problem_t *problem) {
  assert(problem != NULL);
  return problem->number_of_constraints;
}

const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->smallest_values_of_interest != NULL);
  return problem->smallest_values_of_interest;
}

const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->largest_values_of_interest != NULL);
  return problem->largest_values_of_interest;
}

/**
 * If a special method for setting an initial solution to the problem does not exist, the center of the
 * problem's region of interest is the initial solution.
 * @param problem The given COCO problem.
 * @param initial_solution The pointer to the initial solution being set by this method.
 */
void coco_problem_get_initial_solution(const coco_problem_t *problem, double *initial_solution) {
  assert(problem != NULL);
  if (problem->initial_solution != NULL) {
    problem->initial_solution(problem, initial_solution);
  } else {
    size_t i;
    assert(problem->smallest_values_of_interest != NULL);
    assert(problem->largest_values_of_interest != NULL);
    for (i = 0; i < problem->number_of_variables; ++i)
      initial_solution[i] = problem->smallest_values_of_interest[i] + 0.5
          * (problem->largest_values_of_interest[i] - problem->smallest_values_of_interest[i]);
  }
}

static coco_suite_t *coco_problem_get_suite(const coco_problem_t *problem) {
  assert(problem != NULL);
  return problem->suite;
}

static void coco_problem_set_suite(coco_problem_t *problem, const coco_suite_t *suite) {
  assert(problem != NULL);
  problem->suite = suite;
}

size_t coco_problem_get_suite_dep_index(const coco_problem_t *problem) {
  assert(problem != NULL);
  return problem->suite_dep_index;
}

static size_t coco_problem_get_suite_dep_function(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->suite_dep_function > 0);
  return problem->suite_dep_function;
}

static size_t coco_problem_get_suite_dep_instance(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->suite_dep_instance > 0);
  return problem->suite_dep_instance;
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding the transformed COCO problem
 */
/**@{*/

/**
 * @brief Returns the data of the transformed problem.
 */
static void *coco_problem_transformed_get_data(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->data != NULL);
  assert(((coco_problem_transformed_data_t *) problem->data)->data != NULL);

  return ((coco_problem_transformed_data_t *) problem->data)->data;
}

/**
 * @brief Returns the inner problem of the transformed problem.
 */
static coco_problem_t *coco_problem_transformed_get_inner_problem(const coco_problem_t *problem) {
  assert(problem != NULL);
  assert(problem->data != NULL);
  assert(((coco_problem_transformed_data_t *) problem->data)->inner_problem != NULL);

  return ((coco_problem_transformed_data_t *) problem->data)->inner_problem;
}

/**
 * @brief Calls the coco_evaluate_function function on the inner problem.
 */
static void coco_problem_transformed_evaluate_function(coco_problem_t *problem, const double *x, double *y) {
  coco_problem_transformed_data_t *data;
  assert(problem != NULL);
  assert(problem->data != NULL);
  data = (coco_problem_transformed_data_t *) problem->data;
  assert(data->inner_problem != NULL);

  coco_evaluate_function(data->inner_problem, x, y);
}

/**
 * @brief Calls the coco_evaluate_constraint function on the inner problem.
 */
static void coco_problem_transformed_evaluate_constraint(coco_problem_t *problem, const double *x, double *y) {
  coco_problem_transformed_data_t *data;
  assert(problem != NULL);
  assert(problem->data != NULL);
  data = (coco_problem_transformed_data_t *) problem->data;
  assert(data->inner_problem != NULL);

  coco_evaluate_constraint(data->inner_problem, x, y);
}

/**
 * @brief Calls the coco_recommend_solution function on the inner problem.
 */
static void coco_problem_transformed_recommend_solution(coco_problem_t *problem, const double *x) {
  coco_problem_transformed_data_t *data;
  assert(problem != NULL);
  assert(problem->data != NULL);
  data = (coco_problem_transformed_data_t *) problem->data;
  assert(data->inner_problem != NULL);

  coco_recommend_solution(data->inner_problem, x);
}

/**
 * @brief Frees only the data of the transformed problem leaving the inner problem intact.
 *
 * @note If there is no other pointer to the inner problem, access to it will be lost.
 */
static void coco_problem_transformed_free_data(coco_problem_t *problem) {
  coco_problem_transformed_data_t *data;

  assert(problem != NULL);
  assert(problem->data != NULL);
  data = (coco_problem_transformed_data_t *) problem->data;

  if (data->data != NULL) {
    if (data->data_free_function != NULL) {
      data->data_free_function(data->data);
      data->data_free_function = NULL;
    }
    coco_free_memory(data->data);
    data->data = NULL;
  }
  /* Let the generic free problem code deal with the rest of the fields. For this we clear the free_problem
   * function pointer and recall the generic function. */
  problem->problem_free_function = NULL;
  coco_problem_free(problem);
}

/**
 * @brief Frees the transformed problem.
 */
static void coco_problem_transformed_free(coco_problem_t *problem) {
  coco_problem_transformed_data_t *data;

  assert(problem != NULL);
  assert(problem->data != NULL);
  data = (coco_problem_transformed_data_t *) problem->data;
  assert(data->inner_problem != NULL);
  if (data->inner_problem != NULL) {
    coco_problem_free(data->inner_problem);
    data->inner_problem = NULL;
  }
  coco_problem_transformed_free_data(problem);
}

/**
 * @brief Allocates a transformed problem that wraps the inner_problem.
 *
 * By default all methods will dispatch to the inner_problem. A prefix is prepended to the problem name
 * in order to reflect the transformation somewhere.
 */
static coco_problem_t *coco_problem_transformed_allocate(coco_problem_t *inner_problem,
                                                         void *user_data,
                                                         coco_data_free_function_t data_free_function,
                                                         const char *name_prefix) {
  coco_problem_transformed_data_t *problem;
  coco_problem_t *inner_copy;
  char *old_name = coco_strdup(inner_problem->problem_name);

  problem = (coco_problem_transformed_data_t *) coco_allocate_memory(sizeof(*problem));
  problem->inner_problem = inner_problem;
  problem->data = user_data;
  problem->data_free_function = data_free_function;

  inner_copy = coco_problem_duplicate(inner_problem);
  inner_copy->evaluate_function = coco_problem_transformed_evaluate_function;
  inner_copy->evaluate_constraint = coco_problem_transformed_evaluate_constraint;
  inner_copy->recommend_solution = coco_problem_transformed_recommend_solution;
  inner_copy->problem_free_function = coco_problem_transformed_free;
  inner_copy->data = problem;

  coco_problem_set_name(inner_copy, "%s(%s)", name_prefix, old_name);
  coco_free_memory(old_name);

  return inner_copy;
}
/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding the stacked COCO problem
 */
/**@{*/

/**
 * @brief Calls the coco_evaluate_function function on the underlying problems.
 */
static void coco_problem_stacked_evaluate_function(coco_problem_t *problem, const double *x, double *y) {
  coco_problem_stacked_data_t* data = (coco_problem_stacked_data_t *) problem->data;

  assert(
      coco_problem_get_number_of_objectives(problem)
          == coco_problem_get_number_of_objectives(data->problem1)
              + coco_problem_get_number_of_objectives(data->problem2));

  coco_evaluate_function(data->problem1, x, &y[0]);
  coco_evaluate_function(data->problem2, x, &y[coco_problem_get_number_of_objectives(data->problem1)]);
}

/**
 * @brief Calls the coco_evaluate_constraint function on the underlying problems.
 */
static void coco_problem_stacked_evaluate_constraint(coco_problem_t *problem, const double *x, double *y) {
  coco_problem_stacked_data_t* data = (coco_problem_stacked_data_t*) problem->data;

  assert(
      coco_problem_get_number_of_constraints(problem)
          == coco_problem_get_number_of_constraints(data->problem1)
              + coco_problem_get_number_of_constraints(data->problem2));

  if (coco_problem_get_number_of_constraints(data->problem1) > 0)
    coco_evaluate_constraint(data->problem1, x, y);
  if (coco_problem_get_number_of_constraints(data->problem2) > 0)
    coco_evaluate_constraint(data->problem2, x, &y[coco_problem_get_number_of_constraints(data->problem1)]);
}

/* TODO: Missing coco_problem_stacked_recommend_solution function! */

/**
 * @brief Frees the stacked problem.
 */
static void coco_problem_stacked_free(coco_problem_t *problem) {
  coco_problem_stacked_data_t *data;

  assert(problem != NULL);
  assert(problem->data != NULL);
  data = (coco_problem_stacked_data_t*) problem->data;

  if (data->problem1 != NULL) {
    coco_problem_free(data->problem1);
    data->problem1 = NULL;
  }
  if (data->problem2 != NULL) {
    coco_problem_free(data->problem2);
    data->problem2 = NULL;
  }
  /* Let the generic free problem code deal with the rest of the fields. For this we clear the free_problem
   * function pointer and recall the generic function. */
  problem->problem_free_function = NULL;
  coco_problem_free(problem);
}

/**
 * @brief Allocates a problem constructed by stacking two COCO problems.
 * 
 * This is particularly useful for generating multi-objective problems, e.g. a bi-objective problem from two
 * single-objective problems. The stacked problem must behave like a normal COCO problem accepting the same
 * input. The region of interest in the decision space is defined by parameters smallest_values_of_interest
 * and largest_values_of_interest, which are two arrays of size equal to the dimensionality of both problems.
 *
 * @note Regions of interest in the decision space must either agree or at least one of them must be NULL.
 * @note Best parameter becomes somewhat meaningless, but the nadir value make sense now.
 */
static coco_problem_t *coco_problem_stacked_allocate(coco_problem_t *problem1,
																										 coco_problem_t *problem2,
																										 const double *smallest_values_of_interest,
																										 const double *largest_values_of_interest) {

  const size_t number_of_variables = coco_problem_get_dimension(problem1);
  const size_t number_of_objectives = coco_problem_get_number_of_objectives(problem1)
      + coco_problem_get_number_of_objectives(problem2);
  const size_t number_of_constraints = coco_problem_get_number_of_constraints(problem1)
      + coco_problem_get_number_of_constraints(problem2);
  size_t i;
  char *s;
  coco_problem_stacked_data_t *data;
  coco_problem_t *problem; /* the new coco problem */

  assert(coco_problem_get_dimension(problem1) == coco_problem_get_dimension(problem2));

  problem = coco_problem_allocate(number_of_variables, number_of_objectives, number_of_constraints);

  s = coco_strconcat(coco_problem_get_id(problem1), "__");
  problem->problem_id = coco_strconcat(s, coco_problem_get_id(problem2));
  coco_free_memory(s);
  s = coco_strconcat(coco_problem_get_name(problem1), " + ");
  problem->problem_name = coco_strconcat(s, coco_problem_get_name(problem2));
  coco_free_memory(s);

  problem->evaluate_function = coco_problem_stacked_evaluate_function;
  if (number_of_constraints > 0)
    problem->evaluate_constraint = coco_problem_stacked_evaluate_constraint;

	assert(smallest_values_of_interest);
	assert(largest_values_of_interest);
  for (i = 0; i < number_of_variables; ++i) {
    problem->smallest_values_of_interest[i] = smallest_values_of_interest[i];
    problem->largest_values_of_interest[i] = largest_values_of_interest[i];
  }

	if (problem->best_parameter) /* logger_bbob doesn't work then anymore */
		coco_free_memory(problem->best_parameter);
	problem->best_parameter = NULL;

  /* Compute the ideal and nadir values */
  assert(problem->best_value);
  assert(problem->nadir_value);
  problem->best_value[0] = problem1->best_value[0];
  problem->best_value[1] = problem2->best_value[0];
  coco_evaluate_function(problem1, problem2->best_parameter, &problem->nadir_value[0]);
  coco_evaluate_function(problem2, problem1->best_parameter, &problem->nadir_value[1]);

  /* setup data holder */
  data = (coco_problem_stacked_data_t *) coco_allocate_memory(sizeof(*data));
  data->problem1 = problem1;
  data->problem2 = problem2;

  problem->data = data;
  problem->problem_free_function = coco_problem_stacked_free;

  return problem;
}
/**@}*/

/***********************************************************************************************************/
#line 11 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/suite_bbob_legacy_code.c"
/**
 * @file suite_bbob_legacy_code.c
 * @brief Legacy code from BBOB2009 required to replicate the 2009 functions.
 *
 * All of this code should only be used by the suite_bbob2009 functions to provide compatibility to the
 * legacy code. New test beds should strive to use the new COCO facilities for random number generation etc.
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#line 13 "code-experiments/src/suite_bbob_legacy_code.c"

/** @brief Maximal dimension used in BBOB2009. */
#define SUITE_BBOB2009_MAX_DIM 40

/** @brief Computes the minimum of the two values. */
static double bbob2009_fmin(double a, double b) {
  return (a < b) ? a : b;
}

/** @brief Computes the maximum of the two values. */
static double bbob2009_fmax(double a, double b) {
  return (a > b) ? a : b;
}

/** @brief Rounds the given value. */
static double bbob2009_round(double x) {
  return floor(x + 0.5);
}

/**
 * @brief Allocates a n by m matrix structured as an array of pointers to double arrays.
 */
static double **bbob2009_allocate_matrix(const size_t n, const size_t m) {
  double **matrix = NULL;
  size_t i;
  matrix = (double **) coco_allocate_memory(sizeof(double *) * n);
  for (i = 0; i < n; ++i) {
    matrix[i] = coco_allocate_vector(m);
  }
  return matrix;
}

/**
 * @brief Frees the matrix structured as an array of pointers to double arrays.
 */
static void bbob2009_free_matrix(double **matrix, const size_t n) {
  size_t i;
  for (i = 0; i < n; ++i) {
    if (matrix[i] != NULL) {
      coco_free_memory(matrix[i]);
      matrix[i] = NULL;
    }
  }
  coco_free_memory(matrix);
}

/**
 * @brief Generates N uniform random numbers using inseed as the seed and stores them in r.
 */
static void bbob2009_unif(double *r, size_t N, long inseed) {
  /* generates N uniform numbers with starting seed */
  long aktseed;
  long tmp;
  long rgrand[32];
  long aktrand;
  long i;

  if (inseed < 0)
    inseed = -inseed;
  if (inseed < 1)
    inseed = 1;
  aktseed = inseed;
  for (i = 39; i >= 0; i--) {
    tmp = (int) floor((double) aktseed / (double) 127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if (aktseed < 0)
      aktseed = aktseed + 2147483647;
    if (i < 32)
      rgrand[i] = aktseed;
  }
  aktrand = rgrand[0];
  for (i = 0; i < N; i++) {
    tmp = (int) floor((double) aktseed / (double) 127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if (aktseed < 0)
      aktseed = aktseed + 2147483647;
    tmp = (int) floor((double) aktrand / (double) 67108865);
    aktrand = rgrand[tmp];
    rgrand[tmp] = aktseed;
    r[i] = (double) aktrand / 2.147483647e9;
    if (r[i] == 0.) {
      r[i] = 1e-99;
    }
  }
  return;
}

/**
 * @brief Converts from packed matrix storage to an array of array of double representation.
 */
static double **bbob2009_reshape(double **B, double *vector, const size_t m, const size_t n) {
  size_t i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      B[i][j] = vector[j * m + i];
    }
  }
  return B;
}

/**
 * @brief Generates N Gaussian random numbers using the given seed and stores them in g.
 */
static void bbob2009_gauss(double *g, const size_t N, const long seed) {
  size_t i;
  double uniftmp[6000];
  assert(2 * N < 6000);
  bbob2009_unif(uniftmp, 2 * N, seed);

  for (i = 0; i < N; i++) {
    g[i] = sqrt(-2 * log(uniftmp[i])) * cos(2 * coco_pi * uniftmp[N + i]);
    if (g[i] == 0.)
      g[i] = 1e-99;
  }
  return;
}

/**
 * @brief Computes a DIM by DIM rotation matrix based on seed and stores it in B.
 */
static void bbob2009_compute_rotation(double **B, const long seed, const size_t DIM) {
  /* To ensure temporary data fits into gvec */
  double prod;
  double gvect[2000];
  long i, j, k; /* Loop over pairs of column vectors. */

  assert(DIM * DIM < 2000);

  bbob2009_gauss(gvect, DIM * DIM, seed);
  bbob2009_reshape(B, gvect, DIM, DIM);
  /*1st coordinate is row, 2nd is column.*/

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < i; j++) {
      prod = 0;
      for (k = 0; k < DIM; k++)
        prod += B[k][i] * B[k][j];
      for (k = 0; k < DIM; k++)
        B[k][i] -= prod * B[k][j];
    }
    prod = 0;
    for (k = 0; k < DIM; k++)
      prod += B[k][i] * B[k][i];
    for (k = 0; k < DIM; k++)
      B[k][i] /= sqrt(prod);
  }
}

static void bbob2009_copy_rotation_matrix(double **rot, double *M, double *b, const size_t DIM) {
  size_t row, column;
  double *current_row;

  for (row = 0; row < DIM; ++row) {
    current_row = M + row * DIM;
    for (column = 0; column < DIM; ++column) {
      current_row[column] = rot[row][column];
    }
    b[row] = 0.0;
  }
}

/**
 * @brief Randomly computes the location of the global optimum.
 */
static void bbob2009_compute_xopt(double *xopt, const long seed, const size_t DIM) {
  long i;
  bbob2009_unif(xopt, DIM, seed);
  for (i = 0; i < DIM; i++) {
    xopt[i] = 8 * floor(1e4 * xopt[i]) / 1e4 - 4;
    if (xopt[i] == 0.0)
      xopt[i] = -1e-5;
  }
}

/**
 * @brief Randomly chooses the objective offset for the given function and instance.
 */
static double bbob2009_compute_fopt(const size_t function, const size_t instance) {
  long rseed, rrseed;
  double gval, gval2;

  if (function == 4)
    rseed = 3;
  else if (function == 18)
    rseed = 17;
  else if (function == 101 || function == 102 || function == 103 || function == 107
      || function == 108 || function == 109)
    rseed = 1;
  else if (function == 104 || function == 105 || function == 106 || function == 110
      || function == 111 || function == 112)
    rseed = 8;
  else if (function == 113 || function == 114 || function == 115)
    rseed = 7;
  else if (function == 116 || function == 117 || function == 118)
    rseed = 10;
  else if (function == 119 || function == 120 || function == 121)
    rseed = 14;
  else if (function == 122 || function == 123 || function == 124)
    rseed = 17;
  else if (function == 125 || function == 126 || function == 127)
    rseed = 19;
  else if (function == 128 || function == 129 || function == 130)
    rseed = 21;
  else
    rseed = (long) function;

  rrseed = rseed + (long) (10000 * instance);
  bbob2009_gauss(&gval, 1, rrseed);
  bbob2009_gauss(&gval2, 1, rrseed + 1);
  return bbob2009_fmin(1000., bbob2009_fmax(-1000., bbob2009_round(100. * 100. * gval / gval2) / 100.));
}
#line 12 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/transform_obj_oscillate.c"
/**
 * @file transform_obj_oscillate.c
 * @brief Implementation of oscillating the objective value.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/transform_obj_oscillate.c"
#line 11 "code-experiments/src/transform_obj_oscillate.c"

/**
 * @brief Evaluates the transformation.
 */
static void transform_obj_oscillate_evaluate(coco_problem_t *problem, const double *x, double *y) {
  static const double factor = 0.1;
  size_t i;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }
  coco_evaluate_function(coco_problem_transformed_get_inner_problem(problem), x, y);

  for (i = 0; i < problem->number_of_objectives; i++) {
      if (y[i] != 0) {
          double log_y;
          log_y = log(fabs(y[i])) / factor;
          if (y[i] > 0) {
              y[i] = pow(exp(log_y + 0.49 * (sin(log_y) + sin(0.79 * log_y))), factor);
          } else {
              y[i] = -pow(exp(log_y + 0.49 * (sin(0.55 * log_y) + sin(0.31 * log_y))), factor);
          }
      }
  }
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_oscillate(coco_problem_t *inner_problem) {
  coco_problem_t *problem;
  problem = coco_problem_transformed_allocate(inner_problem, NULL, NULL, "transform_obj_oscillate");
  problem->evaluate_function = transform_obj_oscillate_evaluate;
  /* Compute best value */
  /* Maybe not the most efficient solution */
  transform_obj_oscillate_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}
#line 13 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/transform_obj_power.c"
/**
 * @file transform_obj_power.c
 * @brief Implementation of raising the objective value to the power of a given exponent.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/transform_obj_power.c"
#line 11 "code-experiments/src/transform_obj_power.c"

/**
 * @brief Data type for transform_obj_power.
 */
typedef struct {
  double exponent;
} transform_obj_power_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_obj_power_evaluate(coco_problem_t *problem, const double *x, double *y) {
  transform_obj_power_data_t *data;
  size_t i;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_obj_power_data_t *) coco_problem_transformed_get_data(problem);
  coco_evaluate_function(coco_problem_transformed_get_inner_problem(problem), x, y);

  for (i = 0; i < problem->number_of_objectives; i++) {
      y[i] = pow(y[i], data->exponent);
  }
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_power(coco_problem_t *inner_problem, const double exponent) {
  transform_obj_power_data_t *data;
  coco_problem_t *problem;

  data = (transform_obj_power_data_t *) coco_allocate_memory(sizeof(*data));
  data->exponent = exponent;

  problem = coco_problem_transformed_allocate(inner_problem, data, NULL, "transform_obj_power");
  problem->evaluate_function = transform_obj_power_evaluate;
  /* Compute best value */
  transform_obj_power_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}
#line 14 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/transform_obj_shift.c"
/**
 * @file transform_obj_shift.c
 * @brief Implementation of shifting the objective value by the given offset.
 */

#include <assert.h>

#line 9 "code-experiments/src/transform_obj_shift.c"
#line 10 "code-experiments/src/transform_obj_shift.c"

/**
 * @brief Data type for transform_obj_shift.
 */
typedef struct {
  double offset;
} transform_obj_shift_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_obj_shift_evaluate(coco_problem_t *problem, const double *x, double *y) {
  transform_obj_shift_data_t *data;
  size_t i;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_obj_shift_data_t *) coco_problem_transformed_get_data(problem);
  coco_evaluate_function(coco_problem_transformed_get_inner_problem(problem), x, y);

  for (i = 0; i < problem->number_of_objectives; i++) {
      y[i] += data->offset;
  }
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_shift(coco_problem_t *inner_problem, const double offset) {
  coco_problem_t *problem;
  transform_obj_shift_data_t *data;
  size_t i;
  data = (transform_obj_shift_data_t *) coco_allocate_memory(sizeof(*data));
  data->offset = offset;

  problem = coco_problem_transformed_allocate(inner_problem, data, NULL, "transform_obj_shift");
  problem->evaluate_function = transform_obj_shift_evaluate;
  for (i = 0; i < problem->number_of_objectives; i++) {
      problem->best_value[0] += offset;
  }
  return problem;
}
#line 15 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/transform_vars_affine.c"
/**
 * @file transform_vars_affine.c
 * @brief Implementation of performing an affine transformation on decision values.
 *
 * x |-> Mx + b <br>
 * The matrix M is stored in row-major format.
 */

#include <assert.h>

#line 12 "code-experiments/src/transform_vars_affine.c"
#line 13 "code-experiments/src/transform_vars_affine.c"

/**
 * @brief Data type for transform_vars_affine.
 */
typedef struct {
  double *M, *b, *x;
} transform_vars_affine_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_affine_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i, j;
  transform_vars_affine_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_affine_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  for (i = 0; i < inner_problem->number_of_variables; ++i) {
    /* data->M has problem->number_of_variables columns and inner_problem->number_of_variables rows. */
    const double *current_row = data->M + i * problem->number_of_variables;
    data->x[i] = data->b[i];
    for (j = 0; j < problem->number_of_variables; ++j) {
      data->x[i] += x[j] * current_row[j];
    }
  }
  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_affine_free(void *thing) {
  transform_vars_affine_data_t *data = (transform_vars_affine_data_t *) thing;
  coco_free_memory(data->M);
  coco_free_memory(data->b);
  coco_free_memory(data->x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_affine(coco_problem_t *inner_problem,
                                             const double *M,
                                             const double *b,
                                             const size_t number_of_variables) {
  /*
   * TODO:
   * - Calculate new smallest/largest values of interest?
   * - Resize bounds vectors if input and output dimensions do not match
   */

  coco_problem_t *problem;
  transform_vars_affine_data_t *data;
  size_t entries_in_M;

  entries_in_M = inner_problem->number_of_variables * number_of_variables;
  data = (transform_vars_affine_data_t *) coco_allocate_memory(sizeof(*data));
  data->M = coco_duplicate_vector(M, entries_in_M);
  data->b = coco_duplicate_vector(b, inner_problem->number_of_variables);
  data->x = coco_allocate_vector(inner_problem->number_of_variables);
  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_affine_free, "transform_vars_affine");
  problem->evaluate_function = transform_vars_affine_evaluate;
  if (coco_problem_best_parameter_not_zero(inner_problem)) {
    coco_debug("transform_vars_affine(): 'best_parameter' not updated, set to NAN");
    coco_vector_set_to_nan(inner_problem->best_parameter, inner_problem->number_of_variables);
  }
  return problem;
}
#line 16 "code-experiments/src/f_attractive_sector.c"
#line 1 "code-experiments/src/transform_vars_shift.c"
/**
 * @file transform_vars_shift.c
 * @brief Implementation of shifting all decision values by an offset.
 */

#include <assert.h>

#line 9 "code-experiments/src/transform_vars_shift.c"
#line 10 "code-experiments/src/transform_vars_shift.c"

/**
 * @brief Data type for transform_vars_shift.
 */
typedef struct {
  double *offset;
  double *shifted_x;
  coco_problem_free_function_t old_free_problem;
} transform_vars_shift_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_shift_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_shift_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_shift_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  for (i = 0; i < problem->number_of_variables; ++i) {
    data->shifted_x[i] = x[i] - data->offset[i];
  }
  coco_evaluate_function(inner_problem, data->shifted_x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_shift_free(void *thing) {
  transform_vars_shift_data_t *data = (transform_vars_shift_data_t *) thing;
  coco_free_memory(data->shifted_x);
  coco_free_memory(data->offset);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_shift(coco_problem_t *inner_problem,
                                            const double *offset,
                                            const int shift_bounds) {
  transform_vars_shift_data_t *data;
  coco_problem_t *problem;
  size_t i;
  if (shift_bounds)
    coco_error("shift_bounds not implemented.");

  data = (transform_vars_shift_data_t *) coco_allocate_memory(sizeof(*data));
  data->offset = coco_duplicate_vector(offset, inner_problem->number_of_variables);
  data->shifted_x = coco_allocate_vector(inner_problem->number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_shift_free, "transform_vars_shift");
  problem->evaluate_function = transform_vars_shift_evaluate;
  /* Compute best parameter */
  for (i = 0; i < problem->number_of_variables; i++) {
      problem->best_parameter[i] += data->offset[i];
  }
  return problem;
}
#line 17 "code-experiments/src/f_attractive_sector.c"

/**
 * @brief Data type for the attractive sector problem.
 */
typedef struct {
  double *xopt;
} f_attractive_sector_data_t;

/**
 * @brief Implements the attractive sector function without connections to any COCO structures.
 */
static double f_attractive_sector_raw(const double *x,
                                      const size_t number_of_variables,
                                      f_attractive_sector_data_t *data) {
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    if (data->xopt[i] * x[i] > 0.0) {
      result += 100.0 * 100.0 * x[i] * x[i];
    } else {
      result += x[i] * x[i];
    }
  }
  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_attractive_sector_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_attractive_sector_raw(x, problem->number_of_variables, (f_attractive_sector_data_t *) problem->data);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the attractive sector data object.
 */
static void f_attractive_sector_free(coco_problem_t *problem) {
  f_attractive_sector_data_t *data;
  data = (f_attractive_sector_data_t *) problem->data;
  coco_free_memory(data->xopt);
  problem->problem_free_function = NULL;
  coco_problem_free(problem);
}

/**
 * @brief Allocates the basic attractive sector problem.
 */
static coco_problem_t *f_attractive_sector_allocate(const size_t number_of_variables, const double *xopt) {

  f_attractive_sector_data_t *data;
  coco_problem_t *problem = coco_problem_allocate_from_scalars("attractive sector function",
      f_attractive_sector_evaluate, f_attractive_sector_free, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "attractive_sector", number_of_variables);

  data = (f_attractive_sector_data_t *) coco_allocate_memory(sizeof(*data));
  data->xopt = coco_duplicate_vector(xopt, number_of_variables);
  problem->data = data;

  /* Compute best solution */
  f_attractive_sector_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB attractive sector problem.
 */
static coco_problem_t *f_attractive_sector_bbob_problem_allocate(const size_t function,
                                                                 const size_t dimension,
                                                                 const size_t instance,
                                                                 const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j, k;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, **rot2;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  /* Compute affine transformation M from two rotation matrices */
  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  rot2 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2, rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      for (k = 0; k < dimension; ++k) {
        double exponent = 1.0 * (int) k / ((double) (long) dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(sqrt(10.0), exponent) * rot2[k][j];
      }
    }
  }
  bbob2009_free_matrix(rot1, dimension);
  bbob2009_free_matrix(rot2, dimension);

  problem = f_attractive_sector_allocate(dimension, xopt);
  problem = transform_obj_oscillate(problem);
  problem = transform_obj_power(problem, 0.9);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "2-moderate");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 10 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_bent_cigar.c"
/**
 * @file f_bent_cigar.c
 * @brief Implementation of the bent cigar function and problem.
 */

#include <stdio.h>
#include <assert.h>

#line 10 "code-experiments/src/f_bent_cigar.c"
#line 11 "code-experiments/src/f_bent_cigar.c"
#line 12 "code-experiments/src/f_bent_cigar.c"
#line 13 "code-experiments/src/f_bent_cigar.c"
#line 14 "code-experiments/src/f_bent_cigar.c"
#line 1 "code-experiments/src/transform_vars_asymmetric.c"
/**
 * @file transform_vars_asymmetric.c
 * @brief Implementation of performing an asymmetric transformation on decision values.
 */

#include <math.h>
#include <assert.h>

#line 10 "code-experiments/src/transform_vars_asymmetric.c"
#line 11 "code-experiments/src/transform_vars_asymmetric.c"

/**
 * @brief Data type for transform_vars_asymmetric.
 */
typedef struct {
  double *x;
  double beta;
} transform_vars_asymmetric_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_asymmetric_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  double exponent;
  transform_vars_asymmetric_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_asymmetric_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  for (i = 0; i < problem->number_of_variables; ++i) {
    if (x[i] > 0.0) {
      exponent = 1.0
          + (data->beta * (double) (long) i) / ((double) (long) problem->number_of_variables - 1.0) * sqrt(x[i]);
      data->x[i] = pow(x[i], exponent);
    } else {
      data->x[i] = x[i];
    }
  }
  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

static void transform_vars_asymmetric_free(void *thing) {
  transform_vars_asymmetric_data_t *data = (transform_vars_asymmetric_data_t *) thing;
  coco_free_memory(data->x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_asymmetric(coco_problem_t *inner_problem, const double beta) {
  transform_vars_asymmetric_data_t *data;
  coco_problem_t *problem;

  data = (transform_vars_asymmetric_data_t *) coco_allocate_memory(sizeof(*data));
  data->x = coco_allocate_vector(inner_problem->number_of_variables);
  data->beta = beta;
  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_asymmetric_free, "transform_vars_asymmetric");
  problem->evaluate_function = transform_vars_asymmetric_evaluate;
  if (coco_problem_best_parameter_not_zero(inner_problem)) {
    coco_warning("transform_vars_asymmetric(): 'best_parameter' not updated, set to NAN");
    coco_vector_set_to_nan(inner_problem->best_parameter, inner_problem->number_of_variables);
  }
  return problem;
}
#line 15 "code-experiments/src/f_bent_cigar.c"
#line 16 "code-experiments/src/f_bent_cigar.c"

/**
 * @brief Implements the bent cigar function without connections to any COCO structures.
 */
static double f_bent_cigar_raw(const double *x, const size_t number_of_variables) {

  static const double condition = 1.0e6;
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = x[0] * x[0];
  for (i = 1; i < number_of_variables; ++i) {
    result += condition * x[i] * x[i];
  }
  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_bent_cigar_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_bent_cigar_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic bent cigar problem.
 */
static coco_problem_t *f_bent_cigar_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("bent cigar function",
      f_bent_cigar_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "bent_cigar", number_of_variables);

  /* Compute best solution */
  f_bent_cigar_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB bent cigar problem.
 */
static coco_problem_t *f_bent_cigar_bbob_problem_allocate(const size_t function,
                                                          const size_t dimension,
                                                          const size_t instance,
                                                          const long rseed,
                                                          const char *problem_id_template,
                                                          const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double **rot1;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed + 1000000, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  bbob2009_free_matrix(rot1, dimension);

  problem = f_bent_cigar_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_asymmetric(problem, 0.5);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 11 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_bueche_rastrigin.c"
/**
 * @file f_bueche_rastrigin.c
 * @brief Implementation of the Bueche-Rastrigin function and problem.
 */

#include <math.h>
#include <assert.h>

#line 10 "code-experiments/src/f_bueche_rastrigin.c"
#line 11 "code-experiments/src/f_bueche_rastrigin.c"
#line 12 "code-experiments/src/f_bueche_rastrigin.c"
#line 1 "code-experiments/src/transform_vars_brs.c"
/**
 * @file transform_vars_brs.c
 * @brief Implementation of the ominous 's_i scaling' of the BBOB Bueche-Rastrigin problem.
 */

#include <math.h>
#include <assert.h>

#line 10 "code-experiments/src/transform_vars_brs.c"
#line 11 "code-experiments/src/transform_vars_brs.c"

/**
 * @brief Data type for transform_vars_brs.
 */
typedef struct {
  double *x;
} transform_vars_brs_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_brs_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  double factor;
  transform_vars_brs_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_brs_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  for (i = 0; i < problem->number_of_variables; ++i) {
    /* Function documentation says we should compute 10^(0.5 *
     * (i-1)/(D-1)). Instead we compute the equivalent
     * sqrt(10)^((i-1)/(D-1)) just like the legacy code.
     */
    factor = pow(sqrt(10.0), (double) (long) i / ((double) (long) problem->number_of_variables - 1.0));
    /* Documentation specifies odd indices and starts indexing
     * from 1, we use all even indices since C starts indexing
     * with 0.
     */
    if (x[i] > 0.0 && i % 2 == 0) {
      factor *= 10.0;
    }
    data->x[i] = factor * x[i];
  }
  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_brs_free(void *thing) {
  transform_vars_brs_data_t *data = (transform_vars_brs_data_t *) thing;
  coco_free_memory(data->x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_brs(coco_problem_t *inner_problem) {
  transform_vars_brs_data_t *data;
  coco_problem_t *problem;

  data = (transform_vars_brs_data_t *) coco_allocate_memory(sizeof(*data));
  data->x = coco_allocate_vector(inner_problem->number_of_variables);
  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_brs_free, "transform_vars_brs");
  problem->evaluate_function = transform_vars_brs_evaluate;

  if (coco_problem_best_parameter_not_zero(inner_problem)) {
    coco_warning("transform_vars_brs(): 'best_parameter' not updated, set to NAN");
    coco_vector_set_to_nan(inner_problem->best_parameter, inner_problem->number_of_variables);
  }
  return problem;
}
#line 13 "code-experiments/src/f_bueche_rastrigin.c"
#line 1 "code-experiments/src/transform_vars_oscillate.c"
/**
 * @file transform_vars_oscillate.c
 * @brief Implementation of oscillating the decision values.
 */

#include <math.h>
#include <assert.h>

#line 10 "code-experiments/src/transform_vars_oscillate.c"
#line 11 "code-experiments/src/transform_vars_oscillate.c"

/**
 * @brief Data type for transform_vars_oscillate.
 */
typedef struct {
  double *oscillated_x;
} transform_vars_oscillate_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_oscillate_evaluate(coco_problem_t *problem, const double *x, double *y) {
  static const double alpha = 0.1;
  double tmp, base, *oscillated_x;
  size_t i;
  transform_vars_oscillate_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_oscillate_data_t *) coco_problem_transformed_get_data(problem);
  oscillated_x = data->oscillated_x; /* short cut to make code more readable */
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  for (i = 0; i < problem->number_of_variables; ++i) {
    if (x[i] > 0.0) {
      tmp = log(x[i]) / alpha;
      base = exp(tmp + 0.49 * (sin(tmp) + sin(0.79 * tmp)));
      oscillated_x[i] = pow(base, alpha);
    } else if (x[i] < 0.0) {
      tmp = log(-x[i]) / alpha;
      base = exp(tmp + 0.49 * (sin(0.55 * tmp) + sin(0.31 * tmp)));
      oscillated_x[i] = -pow(base, alpha);
    } else {
      oscillated_x[i] = 0.0;
    }
  }
  coco_evaluate_function(inner_problem, oscillated_x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_oscillate_free(void *thing) {
  transform_vars_oscillate_data_t *data = (transform_vars_oscillate_data_t *) thing;
  coco_free_memory(data->oscillated_x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_oscillate(coco_problem_t *inner_problem) {
  transform_vars_oscillate_data_t *data;
  coco_problem_t *problem;
  data = (transform_vars_oscillate_data_t *) coco_allocate_memory(sizeof(*data));
  data->oscillated_x = coco_allocate_vector(inner_problem->number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_oscillate_free, "transform_vars_oscillate");
  problem->evaluate_function = transform_vars_oscillate_evaluate;
  return problem;
}
#line 14 "code-experiments/src/f_bueche_rastrigin.c"
#line 15 "code-experiments/src/f_bueche_rastrigin.c"
#line 16 "code-experiments/src/f_bueche_rastrigin.c"
#line 1 "code-experiments/src/transform_obj_penalize.c"
/**
 * @file transform_obj_penalize.c
 * @brief Implementation of adding a penalty to the objective value for solutions outside of the ROI in the
 * decision space.
 */

#include <assert.h>

#line 10 "code-experiments/src/transform_obj_penalize.c"
#line 11 "code-experiments/src/transform_obj_penalize.c"

/**
 * @brief Data type for transform_obj_penalize.
 */
typedef struct {
  double factor;
} transform_obj_penalize_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_obj_penalize_evaluate(coco_problem_t *problem, const double *x, double *y) {
  transform_obj_penalize_data_t *data = (transform_obj_penalize_data_t *) coco_problem_transformed_get_data(problem);
  const double *lower_bounds = problem->smallest_values_of_interest;
  const double *upper_bounds = problem->largest_values_of_interest;
  double penalty = 0.0;
  size_t i;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  for (i = 0; i < problem->number_of_variables; ++i) {
    const double c1 = x[i] - upper_bounds[i];
    const double c2 = lower_bounds[i] - x[i];
    assert(lower_bounds[i] < upper_bounds[i]);
    if (c1 > 0.0) {
      penalty += c1 * c1;
    } else if (c2 > 0.0) {
      penalty += c2 * c2;
    }
  }
  assert(coco_problem_transformed_get_inner_problem(problem) != NULL);
  coco_evaluate_function(coco_problem_transformed_get_inner_problem(problem), x, y);

  for (i = 0; i < problem->number_of_objectives; ++i) {
    y[i] += data->factor * penalty;
  }
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_penalize(coco_problem_t *inner_problem, const double factor) {
  coco_problem_t *problem;
  transform_obj_penalize_data_t *data;
  assert(inner_problem != NULL);

  data = (transform_obj_penalize_data_t *) coco_allocate_memory(sizeof(*data));
  data->factor = factor;
  problem = coco_problem_transformed_allocate(inner_problem, data, NULL, "transform_obj_penalize");
  problem->evaluate_function = transform_obj_penalize_evaluate;
  /* No need to update the best value as the best parameter is feasible */
  return problem;
}
#line 17 "code-experiments/src/f_bueche_rastrigin.c"

/**
 * @brief Implements the Bueche-Rastrigin function without connections to any COCO structures.
 */
static double f_bueche_rastrigin_raw(const double *x, const size_t number_of_variables) {

  double tmp = 0., tmp2 = 0.;
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    tmp += cos(2 * coco_pi * x[i]);
    tmp2 += x[i] * x[i];
  }
  result = 10.0 * ((double) (long) number_of_variables - tmp) + tmp2 + 0;
  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_bueche_rastrigin_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_bueche_rastrigin_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Bueche-Rastrigin problem.
 */
static coco_problem_t *f_bueche_rastrigin_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Bueche-Rastrigin function",
      f_bueche_rastrigin_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "bueche-rastrigin", number_of_variables);

  /* Compute best solution */
  f_bueche_rastrigin_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Bueche-Rastrigin problem.
 */
static coco_problem_t *f_bueche_rastrigin_bbob_problem_allocate(const size_t function,
                                                                const size_t dimension,
                                                                const size_t instance,
                                                                const long rseed,
                                                                const char *problem_id_template,
                                                                const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  const double penalty_factor = 100.0;
  size_t i;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  /* OME: This step is in the legacy C code but _not_ in the function description. */
  for (i = 0; i < dimension; i += 2) {
    xopt[i] = fabs(xopt[i]);
  }

  problem = f_bueche_rastrigin_allocate(dimension);
  problem = transform_vars_brs(problem);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_obj_penalize(problem, penalty_factor);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "1-separable");

  coco_free_memory(xopt);
  return problem;
}
#line 12 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_different_powers.c"
/**
 * @file f_different_powers.c
 * @brief Implementation of the different powers function and problem.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/f_different_powers.c"
#line 11 "code-experiments/src/f_different_powers.c"
#line 12 "code-experiments/src/f_different_powers.c"
#line 13 "code-experiments/src/f_different_powers.c"
#line 14 "code-experiments/src/f_different_powers.c"
#line 15 "code-experiments/src/f_different_powers.c"

/**
 * @brief Implements the different powers function without connections to any COCO structures.
 */
static double f_different_powers_raw(const double *x, const size_t number_of_variables) {

  size_t i;
  double sum = 0.0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double exponent = 2.0 + (4.0 * (double) (long) i) / ((double) (long) number_of_variables - 1.0);
    sum += pow(fabs(x[i]), exponent);
  }
  result = sqrt(sum);

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_different_powers_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_different_powers_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic different powers problem.
 */
static coco_problem_t *f_different_powers_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("different powers function",
      f_different_powers_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "different_powers", number_of_variables);

  /* Compute best solution */
  f_different_powers_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB different powers problem.
 */
static coco_problem_t *f_different_powers_bbob_problem_allocate(const size_t function,
                                                                const size_t dimension,
                                                                const size_t instance,
                                                                const long rseed,
                                                                const char *problem_id_template,
                                                                const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double **rot1;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  bbob2009_free_matrix(rot1, dimension);

  problem = f_different_powers_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 13 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_discus.c"
/**
 * @file f_discus.c
 * @brief Implementation of the discus function and problem.
 */

#include <assert.h>

#line 9 "code-experiments/src/f_discus.c"
#line 10 "code-experiments/src/f_discus.c"
#line 11 "code-experiments/src/f_discus.c"
#line 12 "code-experiments/src/f_discus.c"
#line 13 "code-experiments/src/f_discus.c"
#line 14 "code-experiments/src/f_discus.c"
#line 15 "code-experiments/src/f_discus.c"

/**
 * @brief Implements the discus function without connections to any COCO structures.
 */
static double f_discus_raw(const double *x, const size_t number_of_variables) {

  static const double condition = 1.0e6;
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = condition * x[0] * x[0];
  for (i = 1; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_discus_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_discus_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic discus problem.
 */
static coco_problem_t *f_discus_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("discus function",
      f_discus_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "discus", number_of_variables);

  /* Compute best solution */
  f_discus_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB discus problem.
 */
static coco_problem_t *f_discus_bbob_problem_allocate(const size_t function,
                                                      const size_t dimension,
                                                      const size_t instance,
                                                      const long rseed,
                                                      const char *problem_id_template,
                                                      const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double **rot1;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  bbob2009_free_matrix(rot1, dimension);

  problem = f_discus_allocate(dimension);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 14 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_ellipsoid.c"
/**
 * @file f_ellipsoid.c
 * @brief Implementation of the ellipsoid function and problem.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#line 11 "code-experiments/src/f_ellipsoid.c"
#line 12 "code-experiments/src/f_ellipsoid.c"
#line 13 "code-experiments/src/f_ellipsoid.c"
#line 14 "code-experiments/src/f_ellipsoid.c"
#line 15 "code-experiments/src/f_ellipsoid.c"
#line 16 "code-experiments/src/f_ellipsoid.c"
#line 1 "code-experiments/src/transform_vars_permblockdiag.c"
/**
 * @file transform_vars_permblockdiag.c
 */

#include <assert.h>

#line 8 "code-experiments/src/transform_vars_permblockdiag.c"
#line 9 "code-experiments/src/transform_vars_permblockdiag.c"
#line 1 "code-experiments/src/large_scale_transformations.c"
#include <stdio.h>
#include <assert.h>
#line 4 "code-experiments/src/large_scale_transformations.c"

#line 6 "code-experiments/src/large_scale_transformations.c"
#line 7 "code-experiments/src/large_scale_transformations.c"

#include <time.h> /*tmp*/

/* TODO: Document this file in doxygen style! */

static double *ls_random_data;/* global variable used to generate the random permutations */

/**
 * ls_allocate_blockmatrix(n, m, bs):
 *
 * Allocate a ${n} by ${m} block matrix of nb_blocks block sizes block_sizes structured as an array of pointers
 * to double arrays.
 * each row constains only the block_sizes[i] possibly non-zero elements
 */
static double **ls_allocate_blockmatrix(const size_t n, const size_t* block_sizes, const size_t nb_blocks) {
  double **matrix = NULL;
  size_t current_blocksize;
  size_t next_bs_change;
  size_t idx_blocksize;
  size_t i;
  size_t sum_block_sizes;
  
  sum_block_sizes = 0;
  for (i = 0; i < nb_blocks; i++){
    sum_block_sizes += block_sizes[i];
  }
  assert(sum_block_sizes == n);
  
  matrix = (double **) coco_allocate_memory(sizeof(double *) * n);
  idx_blocksize = 0;
  next_bs_change = block_sizes[idx_blocksize];
  
  for (i = 0; i < n; ++i) {
    if (i >= next_bs_change) {
      idx_blocksize++;
      next_bs_change += block_sizes[idx_blocksize];
    }
    current_blocksize=block_sizes[idx_blocksize];
    matrix[i] = coco_allocate_vector(current_blocksize);
    
  }
  return matrix;
}


/*
 * frees a block diagonal matrix (same as a matrix but in case of change, easier to update separatly from free_matrix)
 */
static void ls_free_block_matrix(double **matrix, const size_t n) {
  size_t i;
  for (i = 0; i < n; ++i) {
    if (matrix[i] != NULL) {
      coco_free_memory(matrix[i]);
      matrix[i] = NULL;
    }
  }
  coco_free_memory(matrix);
}



/**
 * ls_compute_blockrotation(B, seed, DIM):
 *
 * Compute a ${DIM}x${DIM} block-diagonal matrix based on ${seed} and block_sizes and stores it in ${B}.
 * B is a 2D vector with DIM lines and each line has blocksize(line) elements (the zeros are not stored)
 */
static void ls_compute_blockrotation(double **B, long seed, size_t n, size_t *block_sizes, size_t nb_blocks) {
  double prod;
  /*double *gvect;*/
  double **current_block;
  size_t i, j, k; /* Loop over pairs of column vectors. */
  size_t idx_block, current_blocksize,cumsum_prev_block_sizes, sum_block_sizes;
  size_t nb_entries;
  coco_random_state_t *rng = coco_random_new((uint32_t) seed);
  
  nb_entries = 0;
  sum_block_sizes = 0;
  for (i = 0; i < nb_blocks; i++){
    sum_block_sizes += block_sizes[i];
    nb_entries += block_sizes[i] * block_sizes[i];
  }
  assert(sum_block_sizes == n);
  
  cumsum_prev_block_sizes = 0;/* shift in rows to account for the previous blocks */
  for (idx_block = 0; idx_block < nb_blocks; idx_block++) {
    current_blocksize = block_sizes[idx_block];
    current_block = bbob2009_allocate_matrix(current_blocksize, current_blocksize);
    for (i = 0; i < current_blocksize; i++) {
      for (j = 0; j < current_blocksize; j++) {
        current_block[i][j] = coco_random_normal(rng);
      }
    }
    
    for (i = 0; i < current_blocksize; i++) {
      for (j = 0; j < i; j++) {
        prod = 0;
        for (k = 0; k < current_blocksize; k++){
          prod += current_block[k][i] * current_block[k][j];
        }
        for (k = 0; k < current_blocksize; k++){
          current_block[k][i] -= prod * current_block[k][j];
        }
      }
      prod = 0;
      for (k = 0; k < current_blocksize; k++){
        prod += current_block[k][i] * current_block[k][i];
      }
      for (k = 0; k < current_blocksize; k++){
        current_block[k][i] /= sqrt(prod);
      }
    }
    
    /* now fill the block matrix*/
    for (i = 0 ; i < current_blocksize; i++) {
      for (j = 0; j < current_blocksize; j++) {
        B[i + cumsum_prev_block_sizes][j]=current_block[i][j];
      }
    }
    
    cumsum_prev_block_sizes+=current_blocksize;
    /*current_gvect_pos += current_blocksize * current_blocksize;*/
    ls_free_block_matrix(current_block, current_blocksize);
  }
  /*coco_free_memory(gvect);*/
  coco_random_free(rng);
}

/*
 * makes a copy of a block_matrix
 */
static double **ls_copy_block_matrix(const double *const *B, const size_t dimension, const size_t *block_sizes, const size_t nb_blocks) {
  double **dest;
  size_t i, j, idx_blocksize, current_blocksize, next_bs_change;
  
  dest = ls_allocate_blockmatrix(dimension, block_sizes, nb_blocks);
  idx_blocksize = 0;
  current_blocksize = block_sizes[idx_blocksize];
  next_bs_change = block_sizes[idx_blocksize];
  assert(nb_blocks != 0); /*tmp*//*to silence warning*/
  for (i = 0; i < dimension; i++) {
    if (i >= next_bs_change) {
      idx_blocksize++;
      next_bs_change += block_sizes[idx_blocksize];
    }
    current_blocksize=block_sizes[idx_blocksize];
    for (j = 0; j < current_blocksize; j++) {
      dest[i][j] = B[i][j];
    }
  }
  return dest;
}

/**
 * Comparison function used for sorting.
 * In our case, it serves as a random permutation generator
 */
static int f_compare_doubles_for_random_permutation(const void *a, const void *b) {
  double temp = ls_random_data[*(const size_t *) a] - ls_random_data[*(const size_t *) b];
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

/*
 * generates a random, uniformly sampled, permutation and puts it in P
 */
static void ls_compute_random_permutation(size_t *P, long seed, size_t n) {
  long i;
  coco_random_state_t *rng = coco_random_new((uint32_t) seed);
  ls_random_data = coco_allocate_vector(n);
  for (i = 0; i < n; i++){
    P[i] = (size_t) i;
    ls_random_data[i] = coco_random_uniform(rng);
  }
  qsort(P, n, sizeof(size_t), f_compare_doubles_for_random_permutation);
  coco_random_free(rng);
}


/*
 * returns a uniformly distributed integer between lower_bound and upper_bound using seed.
 */
long ls_rand_int(long lower_bound, long upper_bound, coco_random_state_t *rng){
  long range;
  range = upper_bound - lower_bound + 1;
  return ((long)(coco_random_uniform(rng) * (double) range)) + lower_bound;
}



/*
 * generates a random permutation resulting from nb_swaps truncated uniform swaps of range swap_range
 * missing paramteters: dynamic_not_static pool, seems empirically irrelevant
 * for now so dynamic is implemented (simple since no need for tracking indices
 * if swap_range is the largest possible size_t value ( (size_t) -1 ), a random uniform permutation is generated
 */
static void ls_compute_truncated_uniform_swap_permutation(size_t *P, long seed, size_t n, size_t nb_swaps, size_t swap_range) {
  long i, idx_swap;
  size_t lower_bound, upper_bound, first_swap_var, second_swap_var, tmp;
  size_t *idx_order;
  coco_random_state_t *rng = coco_random_new((uint32_t) seed);

  ls_random_data = coco_allocate_vector(n);
  idx_order = coco_allocate_vector_size_t(n);
  for (i = 0; i < n; i++){
    P[i] = (size_t) i;
    idx_order[i] = (size_t) i;
    ls_random_data[i] = coco_random_uniform(rng);
  }
  
  if (swap_range > 0) {
    /*sort the random data in random_data and arange idx_order accordingly*/
    /*did not use ls_compute_random_permutation to only use the seed once*/
    qsort(idx_order, n, sizeof(size_t), f_compare_doubles_for_random_permutation);
    for (idx_swap = 0; idx_swap < nb_swaps; idx_swap++) {
      first_swap_var = idx_order[idx_swap];
      if (first_swap_var < swap_range) {
        lower_bound = 0;
      }
      else{
        lower_bound = first_swap_var - swap_range;
      }
      if (first_swap_var + swap_range > n - 1) {
        upper_bound = n - 1;
      }
      else{
        upper_bound = first_swap_var + swap_range;
      }

      second_swap_var = (size_t) ls_rand_int((long) lower_bound, (long) upper_bound, rng);
      while (first_swap_var == second_swap_var) {
        second_swap_var = (size_t) ls_rand_int((long) lower_bound, (long) upper_bound, rng);
      }
      /* swap*/
      tmp = P[first_swap_var];
      P[first_swap_var] = P[second_swap_var];
      P[second_swap_var] = tmp;
    }
  } else {
    if ( swap_range == (size_t) -1) {
      /* generate random permutation instead */
      ls_compute_random_permutation(P, seed, n);
    }
    
  }
  coco_random_free(rng);
}



/*
 * duplicates a size_t vector
 */
size_t *coco_duplicate_size_t_vector(const size_t *src, const size_t number_of_elements) {
  size_t i;
  size_t *dst;
  
  assert(src != NULL);
  assert(number_of_elements > 0);
  
  dst = coco_allocate_vector_size_t(number_of_elements);
  for (i = 0; i < number_of_elements; ++i) {
    dst[i] = src[i];
  }
  return dst;
}


/*
 * returns the list of block_sizes and sets nb_blocks to its correct value
 * TODO: update with chosen parameter setting
 */
size_t *ls_get_block_sizes(size_t *nb_blocks, size_t dimension){
  size_t *block_sizes;
  size_t block_size;
  int i;
  
  block_size = coco_double_to_size_t(bbob2009_fmin((double)dimension / 4, 100));
  *nb_blocks = dimension / block_size + ((dimension % block_size) > 0);
  block_sizes = coco_allocate_vector_size_t(*nb_blocks);
  for (i = 0; i < *nb_blocks - 1; i++) {
    block_sizes[i] = block_size;
  }
  block_sizes[*nb_blocks - 1] = dimension - (*nb_blocks - 1) * block_size; /*add rest*/
  return block_sizes;
}


/*
 * return the swap_range corresponding to the problem
 * TODO: update with chosen parameter setting
 */
size_t ls_get_swap_range(size_t dimension){
  return dimension / 3;
}


/*
 * return the number of swaps corresponding to the problem
 * TODO: update with chosen parameter setting
 */
size_t ls_get_nb_swaps(size_t dimension){
  return dimension;
}




#line 10 "code-experiments/src/transform_vars_permblockdiag.c"

/**
 * @brief Data type for transform_vars_permblockdiag.
 */
typedef struct {
  double **B;
  double *x;
  size_t *P1; /*permutation matrices, P1 for the columns of B and P2 for its rows*/
  size_t *P2;
  size_t *block_sizes;
  size_t nb_blocks;
  size_t *block_size_map; /* maps rows to blocksizes, keep until better way is found */
  size_t *first_non_zero_map; /* maps a row to the index of its first non zero element */
} transform_vars_permblockdiag_t;

static void transform_vars_permblockdiag_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i, j, current_blocksize, first_non_zero_ind;
  transform_vars_permblockdiag_t *data;
  coco_problem_t *inner_problem;
  
  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_permblockdiag_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);
  
  for (i = 0; i < inner_problem->number_of_variables; ++i) {
    current_blocksize = data->block_size_map[data->P2[i]];/*the block_size is that of the permuted line*/
    first_non_zero_ind = data->first_non_zero_map[data->P2[i]];
    data->x[i] = 0;
    /*compute data->x[i] = < B[P2[i]] , x[P1] >  */
    for (j = first_non_zero_ind; j < first_non_zero_ind + current_blocksize; ++j) {/*blocksize[P2[i]]*/
      data->x[i] += data->B[data->P2[i]][j - first_non_zero_ind] * x[data->P1[j]];/*all B lines start at 0*/
    }
    if (data->x[i] > 100 || data->x[i] < -100 || 1) {
    }
    
  }
  
  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

static void transform_vars_permblockdiag_free(void *thing) {
  transform_vars_permblockdiag_t *data = (transform_vars_permblockdiag_t *) thing;
  coco_free_memory(data->B);
  coco_free_memory(data->P1);
  coco_free_memory(data->P2);
  coco_free_memory(data->block_sizes);
  coco_free_memory(data->x);
  coco_free_memory(data->block_size_map);
}

/*
 * Apply a double permuted orthogonal block-diagonal transfromation matrix to the search space
 *
 *
 * The matrix M is stored in row-major format.
 */
static coco_problem_t *transform_vars_permblockdiag(coco_problem_t *inner_problem,
                                                    const double * const *B,
                                                    const size_t *P1,
                                                    const size_t *P2,
                                                    const size_t number_of_variables,
                                                    const size_t *block_sizes,
                                                    const size_t nb_blocks) {
  coco_problem_t *problem;
  transform_vars_permblockdiag_t *data;
  size_t entries_in_M, idx_blocksize, next_bs_change, current_blocksize;
  int i;
  entries_in_M = 0;
  assert(number_of_variables > 0);/*tmp*/
  for (i = 0; i < nb_blocks; i++) {
    entries_in_M += block_sizes[i] * block_sizes[i];
  }
  data = (transform_vars_permblockdiag_t *) coco_allocate_memory(sizeof(*data));
  data->B = ls_copy_block_matrix(B, number_of_variables, block_sizes, nb_blocks);
  data->x = coco_allocate_vector(inner_problem->number_of_variables);
  data->P1 = coco_duplicate_size_t_vector(P1, inner_problem->number_of_variables);
  data->P2 = coco_duplicate_size_t_vector(P2, inner_problem->number_of_variables);
  data->block_sizes = coco_duplicate_size_t_vector(block_sizes, nb_blocks);
  data->nb_blocks = nb_blocks;
  data->block_size_map = coco_allocate_vector_size_t(number_of_variables);
  data->first_non_zero_map = coco_allocate_vector_size_t(number_of_variables);
  
  idx_blocksize = 0;
  next_bs_change = block_sizes[idx_blocksize];
  for (i = 0; i < number_of_variables; i++) {
    if (i >= next_bs_change) {
      idx_blocksize++;
      next_bs_change += block_sizes[idx_blocksize];
    }
    current_blocksize=block_sizes[idx_blocksize];
    data->block_size_map[i] = current_blocksize;
    data->first_non_zero_map[i] = next_bs_change - current_blocksize;/* next_bs_change serves also as a cumsum for blocksizes*/
  }
  
  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_permblockdiag_free, "transform_vars_permblockdiag");
  problem->evaluate_function = transform_vars_permblockdiag_evaluate;
  return problem;
}


#line 17 "code-experiments/src/f_ellipsoid.c"

/**
 * @brief Implements the ellipsoid function without connections to any COCO structures.
 */
static double f_ellipsoid_raw(const double *x, const size_t number_of_variables) {

  static const double condition = 1.0e6;
  size_t i = 0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = x[i] * x[i];
  for (i = 1; i < number_of_variables; ++i) {
    const double exponent = 1.0 * (double) (long) i / ((double) (long) number_of_variables - 1.0);
    result += pow(condition, exponent) * x[i] * x[i];
  }

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_ellipsoid_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_ellipsoid_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic ellipsoid problem.
 */
static coco_problem_t *f_ellipsoid_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("ellipsoid function",
      f_ellipsoid_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "ellipsoid", number_of_variables);

  /* Compute best solution */
  f_ellipsoid_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB ellipsoid problem.
 */
static coco_problem_t *f_ellipsoid_bbob_problem_allocate(const size_t function,
                                                         const size_t dimension,
                                                         const size_t instance,
                                                         const long rseed,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  problem = f_ellipsoid_allocate(dimension);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "1-separable");

  coco_free_memory(xopt);
  return problem;
}

/**
 * @brief Creates the BBOB rotated ellipsoid problem.
 */
static coco_problem_t *f_ellipsoid_rotated_bbob_problem_allocate(const size_t function,
                                                                 const size_t dimension,
                                                                 const size_t instance,
                                                                 const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double **rot1;

  xopt = coco_allocate_vector(dimension);
  bbob2009_compute_xopt(xopt, rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  bbob2009_free_matrix(rot1, dimension);

  problem = f_ellipsoid_allocate(dimension);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}

static coco_problem_t *f_ellipsoid_permblockdiag_bbob_problem_allocate(const size_t function,
                                                                       const size_t dimension,
                                                                       const size_t instance,
                                                                       const long rseed,
                                                                       const char *problem_id_template,
                                                                       const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  double **B;
  const double *const *B_copy;
  size_t *P1 = coco_allocate_vector_size_t(dimension);
  size_t *P2 = coco_allocate_vector_size_t(dimension);
  size_t *block_sizes;
  size_t nb_blocks;
  size_t swap_range;
  size_t nb_swaps;
  
  block_sizes = ls_get_block_sizes(&nb_blocks, dimension);
  swap_range = ls_get_swap_range(dimension);
  nb_swaps = ls_get_nb_swaps(dimension);

  /*printf("f:%zu  n:%zu  i:%zu  bs:[%zu,...,%zu,%zu]  sR:%zu\n", function, dimension, instance, block_sizes[0], block_sizes[0],block_sizes[nb_blocks-1], swap_range);*/
  
  xopt = coco_allocate_vector(dimension);
  bbob2009_compute_xopt(xopt, rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  
  B = ls_allocate_blockmatrix(dimension, block_sizes, nb_blocks);
  B_copy = (const double *const *)B;/*TODO: silences the warning, not sure if it prevents the modification of B at all levels*/

  ls_compute_blockrotation(B, rseed + 1000000, dimension, block_sizes, nb_blocks);
  ls_compute_truncated_uniform_swap_permutation(P1, rseed + 2000000, dimension, nb_swaps, swap_range);
  ls_compute_truncated_uniform_swap_permutation(P2, rseed + 3000000, dimension, nb_swaps, swap_range);

  
  problem = f_ellipsoid_allocate(dimension);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_permblockdiag(problem, B_copy, P1, P2, dimension, block_sizes, nb_blocks);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "large_scale_block_rotated");/*TODO: no large scale prefix*/

  ls_free_block_matrix(B, dimension);
  coco_free_memory(P1);
  coco_free_memory(P2);
  coco_free_memory(block_sizes);
  
  return problem;
}




#line 15 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_gallagher.c"
/**
 * @file f_gallagher.c
 * @brief Implementation of the Gallagher function and problem.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/f_gallagher.c"
#line 11 "code-experiments/src/f_gallagher.c"
#line 12 "code-experiments/src/f_gallagher.c"
#line 13 "code-experiments/src/f_gallagher.c"
#line 14 "code-experiments/src/f_gallagher.c"

/**
 * @brief A random permutation type for the Gallagher problem.
 *
 * Needed to create a random permutation that is compatible with the old code.
 */
typedef struct {
  double value;
  size_t index;
} f_gallagher_permutation_t;

/**
 * @brief Data type for the Gallagher problem.
 */
typedef struct {
  long rseed;
  double *xopt;
  double **rotation, **x_local, **arr_scales;
  size_t number_of_peaks;
  double *peak_values;
  coco_problem_free_function_t old_free_problem;
} f_gallagher_data_t;

/**
 * Comparison function used for sorting.
 */
static int f_gallagher_compare_doubles(const void *a, const void *b) {
  double temp = (*(f_gallagher_permutation_t *) a).value - (*(f_gallagher_permutation_t *) b).value;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

/**
 * @brief Implements the Gallagher function without connections to any COCO structures.
 */
static double f_gallagher_raw(const double *x, const size_t number_of_variables, f_gallagher_data_t *data) {
  size_t i, j; /* Loop over dim */
  double *tmx;
  double a = 0.1;
  double tmp2, f = 0., f_add, tmp, f_pen = 0., f_true = 0.;
  double fac;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  fac = -0.5 / (double) number_of_variables;

  /* Boundary handling */
  for (i = 0; i < number_of_variables; ++i) {
    tmp = fabs(x[i]) - 5.;
    if (tmp > 0.) {
      f_pen += tmp * tmp;
    }
  }
  f_add = f_pen;
  /* Transformation in search space */
  /* TODO: this should rather be done in f_gallagher */
  tmx = coco_allocate_vector(number_of_variables);
  for (i = 0; i < number_of_variables; i++) {
    tmx[i] = 0;
    for (j = 0; j < number_of_variables; ++j) {
      tmx[i] += data->rotation[i][j] * x[j];
    }
  }
  /* Computation core*/
  for (i = 0; i < data->number_of_peaks; ++i) {
    tmp2 = 0.;
    for (j = 0; j < number_of_variables; ++j) {
      tmp = (tmx[j] - data->x_local[j][i]);
      tmp2 += data->arr_scales[i][j] * tmp * tmp;
    }
    tmp2 = data->peak_values[i] * exp(fac * tmp2);
    f = coco_double_max(f, tmp2);
  }

  f = 10. - f;
  if (f > 0) {
    f_true = log(f) / a;
    f_true = pow(exp(f_true + 0.49 * (sin(f_true) + sin(0.79 * f_true))), a);
  } else if (f < 0) {
    f_true = log(-f) / a;
    f_true = -pow(exp(f_true + 0.49 * (sin(0.55 * f_true) + sin(0.31 * f_true))), a);
  } else
    f_true = f;

  f_true *= f_true;
  f_true += f_add;
  result = f_true;
  coco_free_memory(tmx);
  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_gallagher_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_gallagher_raw(x, problem->number_of_variables, (f_gallagher_data_t *) problem->data);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the Gallagher data object.
 */
static void f_gallagher_free(coco_problem_t *problem) {
  f_gallagher_data_t *data;
  data = (f_gallagher_data_t *) problem->data;
  coco_free_memory(data->xopt);
  coco_free_memory(data->peak_values);
  bbob2009_free_matrix(data->rotation, problem->number_of_variables);
  bbob2009_free_matrix(data->x_local, problem->number_of_variables);
  bbob2009_free_matrix(data->arr_scales, data->number_of_peaks);
  problem->problem_free_function = NULL;
  coco_problem_free(problem);
}

/**
 * @brief Creates the BBOB Gallagher problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_t *f_gallagher_bbob_problem_allocate(const size_t function,
                                                         const size_t dimension,
                                                         const size_t instance,
                                                         const long rseed,
                                                         const size_t number_of_peaks,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {

  f_gallagher_data_t *data;
  /* problem_name and best_parameter will be overwritten below */
  coco_problem_t *problem = coco_problem_allocate_from_scalars("Gallagher function",
      f_gallagher_evaluate, f_gallagher_free, dimension, -5.0, 5.0, 0.0);

  const size_t peaks_21 = 21;
  const size_t peaks_101 = 101;

  double fopt;
  size_t i, j, k;
  double maxcondition = 1000.;
  /* maxcondition1 satisfies the old code and the doc but seems wrong in that it is, with very high
   * probability, not the largest condition level!!! */
  double maxcondition1 = 1000.;
  double *arrCondition;
  double fitvalues[2] = { 1.1, 9.1 };
  /* Parameters for generating local optima. In the old code, they are different in f21 and f22 */
  double b, c;
  /* Random permutation */
  f_gallagher_permutation_t *rperm;
  double *random_numbers;

  data = (f_gallagher_data_t *) coco_allocate_memory(sizeof(*data));
  /* Allocate temporary storage and space for the rotation matrices */
  data->number_of_peaks = number_of_peaks;
  data->xopt = coco_allocate_vector(dimension);
  data->rotation = bbob2009_allocate_matrix(dimension, dimension);
  data->x_local = bbob2009_allocate_matrix(dimension, number_of_peaks);
  data->arr_scales = bbob2009_allocate_matrix(number_of_peaks, dimension);

  if (number_of_peaks == peaks_101) {
    maxcondition1 = sqrt(maxcondition1);
    b = 10.;
    c = 5.;
  } else if (number_of_peaks == peaks_21) {
    b = 9.8;
    c = 4.9;
  } else {
    coco_error("f_gallagher_bbob_problem_allocate(): '%lu' is a non-supported number of peaks",
    		(unsigned long) number_of_peaks);
  }
  data->rseed = rseed;
  bbob2009_compute_rotation(data->rotation, rseed, dimension);

  /* Initialize all the data of the inner problem */
  random_numbers = coco_allocate_vector(number_of_peaks * dimension); /* This is large enough for all cases below */
  bbob2009_unif(random_numbers, number_of_peaks - 1, data->rseed);
  rperm = (f_gallagher_permutation_t *) coco_allocate_memory(sizeof(*rperm) * (number_of_peaks - 1));
  for (i = 0; i < number_of_peaks - 1; ++i) {
    rperm[i].value = random_numbers[i];
    rperm[i].index = i;
  }
  qsort(rperm, number_of_peaks - 1, sizeof(*rperm), f_gallagher_compare_doubles);

  /* Random permutation */
  arrCondition = coco_allocate_vector(number_of_peaks);
  arrCondition[0] = maxcondition1;
  data->peak_values = coco_allocate_vector(number_of_peaks);
  data->peak_values[0] = 10;
  for (i = 1; i < number_of_peaks; ++i) {
    arrCondition[i] = pow(maxcondition, (double) (rperm[i - 1].index) / ((double) (number_of_peaks - 2)));
    data->peak_values[i] = (double) (i - 1) / (double) (number_of_peaks - 2) * (fitvalues[1] - fitvalues[0])
        + fitvalues[0];
  }
  coco_free_memory(rperm);

  rperm = (f_gallagher_permutation_t *) coco_allocate_memory(sizeof(*rperm) * dimension);
  for (i = 0; i < number_of_peaks; ++i) {
    bbob2009_unif(random_numbers, dimension, data->rseed + (long) (1000 * i));
    for (j = 0; j < dimension; ++j) {
      rperm[j].value = random_numbers[j];
      rperm[j].index = j;
    }
    qsort(rperm, dimension, sizeof(*rperm), f_gallagher_compare_doubles);
    for (j = 0; j < dimension; ++j) {
      data->arr_scales[i][j] = pow(arrCondition[i],
          ((double) rperm[j].index) / ((double) (dimension - 1)) - 0.5);
    }
  }
  coco_free_memory(rperm);

  bbob2009_unif(random_numbers, dimension * number_of_peaks, data->rseed);
  for (i = 0; i < dimension; ++i) {
    data->xopt[i] = 0.8 * (b * random_numbers[i] - c);
    problem->best_parameter[i] = 0.8 * (b * random_numbers[i] - c);
    for (j = 0; j < number_of_peaks; ++j) {
      data->x_local[i][j] = 0.;
      for (k = 0; k < dimension; ++k) {
        data->x_local[i][j] += data->rotation[i][k] * (b * random_numbers[j * dimension + k] - c);
      }
      if (j == 0) {
        data->x_local[i][j] *= 0.8;
      }
    }
  }
  coco_free_memory(arrCondition);
  coco_free_memory(random_numbers);

  problem->data = data;

  /* Compute best solution */
  f_gallagher_evaluate(problem, problem->best_parameter, problem->best_value);

  fopt = bbob2009_compute_fopt(function, instance);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "5-weakly-structured");

  return problem;
}

#line 16 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_griewank_rosenbrock.c"
/**
 * @file f_griewank_rosenbrock.c
 * @brief Implementation of the Griewank-Rosenbrock function and problem.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#line 11 "code-experiments/src/f_griewank_rosenbrock.c"
#line 12 "code-experiments/src/f_griewank_rosenbrock.c"
#line 13 "code-experiments/src/f_griewank_rosenbrock.c"
#line 14 "code-experiments/src/f_griewank_rosenbrock.c"
#line 15 "code-experiments/src/f_griewank_rosenbrock.c"
#line 16 "code-experiments/src/f_griewank_rosenbrock.c"

/**
 * @brief Implements the Griewank-Rosenbrock function without connections to any COCO structures.
 */
static double f_griewank_rosenbrock_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;
  double tmp = 0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  /* Computation core */
  result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double c1 = x[i] * x[i] - x[i + 1];
    const double c2 = 1.0 - x[i];
    tmp = 100.0 * c1 * c1 + c2 * c2;
    result += tmp / 4000. - cos(tmp);
  }
  result = 10. + 10. * result / (double) (number_of_variables - 1);

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_griewank_rosenbrock_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_griewank_rosenbrock_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Griewank-Rosenbrock problem.
 */
static coco_problem_t *f_griewank_rosenbrock_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Griewank Rosenbrock function",
      f_griewank_rosenbrock_evaluate, NULL, number_of_variables, -5.0, 5.0, 1);
  coco_problem_set_id(problem, "%s_d%02lu", "griewank_rosenbrock", number_of_variables);

  /* Compute best solution */
  f_griewank_rosenbrock_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Griewank-Rosenbrock problem.
 */
static coco_problem_t *f_griewank_rosenbrock_bbob_problem_allocate(const size_t function,
                                                                   const size_t dimension,
                                                                   const size_t instance,
                                                                   const long rseed,
                                                                   const char *problem_id_template,
                                                                   const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = NULL;
  size_t i, j;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *shift = coco_allocate_vector(dimension);
  double scales, **rot1;

  fopt = bbob2009_compute_fopt(function, instance);
  for (i = 0; i < dimension; ++i) {
    shift[i] = -0.5;
  }

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed, dimension);
  scales = coco_double_max(1., sqrt((double) dimension) / 8.);
  for (i = 0; i < dimension; ++i) {
    for (j = 0; j < dimension; ++j) {
      rot1[i][j] *= scales;
    }
  }

  problem = f_griewank_rosenbrock_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_shift(problem, shift, 0);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  problem = transform_vars_affine(problem, M, b, dimension);

  bbob2009_free_matrix(rot1, dimension);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(shift);
  return problem;
}
#line 17 "code-experiments/src/suite_bbob.c"
#line 18 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_katsuura.c"
/**
 * @file f_katsuura.c
 * @brief Implementation of the Katsuura function and problem.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#line 11 "code-experiments/src/f_katsuura.c"
#line 12 "code-experiments/src/f_katsuura.c"
#line 13 "code-experiments/src/f_katsuura.c"
#line 14 "code-experiments/src/f_katsuura.c"
#line 15 "code-experiments/src/f_katsuura.c"
#line 16 "code-experiments/src/f_katsuura.c"
#line 17 "code-experiments/src/f_katsuura.c"
#line 18 "code-experiments/src/f_katsuura.c"

/**
 * @brief Implements the Katsuura function without connections to any COCO structures.
 */
static double f_katsuura_raw(const double *x, const size_t number_of_variables) {

  size_t i, j;
  double tmp, tmp2;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  /* Computation core */
  result = 1.0;
  for (i = 0; i < number_of_variables; ++i) {
    tmp = 0;
    for (j = 1; j < 33; ++j) {
      tmp2 = pow(2., (double) j);
      tmp += fabs(tmp2 * x[i] - coco_double_round(tmp2 * x[i])) / tmp2;
    }
    tmp = 1.0 + ((double) (long) i + 1) * tmp;
    result *= tmp;
  }
  result = 10. / ((double) number_of_variables) / ((double) number_of_variables)
      * (-1. + pow(result, 10. / pow((double) number_of_variables, 1.2)));

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_katsuura_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_katsuura_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Katsuura problem.
 */
static coco_problem_t *f_katsuura_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Katsuura function",
      f_katsuura_evaluate, NULL, number_of_variables, -5.0, 5.0, 1);
  coco_problem_set_id(problem, "%s_d%02lu", "katsuura", number_of_variables);

  /* Compute best solution */
  f_katsuura_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Katsuura problem.
 */
static coco_problem_t *f_katsuura_bbob_problem_allocate(const size_t function,
                                                        const size_t dimension,
                                                        const size_t instance,
                                                        const long rseed,
                                                        const char *problem_id_template,
                                                        const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j, k;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, **rot2;

  const double penalty_factor = 1.0;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  rot2 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2, rseed, dimension);

  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      for (k = 0; k < dimension; ++k) {
        double exponent = 1.0 * (int) k / ((double) (long) dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(sqrt(100), exponent) * rot2[k][j];
      }
    }
  }

  problem = f_katsuura_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_penalize(problem, penalty_factor);

  bbob2009_free_matrix(rot1, dimension);
  bbob2009_free_matrix(rot2, dimension);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "5-weakly-structured");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 19 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_linear_slope.c"
/**
 * @file f_linear_slope.c
 * @brief Implementation of the linear slope function and problem.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#line 11 "code-experiments/src/f_linear_slope.c"
#line 12 "code-experiments/src/f_linear_slope.c"
#line 13 "code-experiments/src/f_linear_slope.c"
#line 14 "code-experiments/src/f_linear_slope.c"

/**
 * @brief Implements the linear slope function without connections to any COCO structures.
 */
static double f_linear_slope_raw(const double *x,
                                 const size_t number_of_variables,
                                 const double *best_parameter) {

  static const double alpha = 100.0;
  size_t i;
  double result = 0.0;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double base, exponent, si;

    base = sqrt(alpha);
    exponent = (double) (long) i / ((double) (long) number_of_variables - 1);
    if (best_parameter[i] > 0.0) {
      si = pow(base, exponent);
    } else {
      si = -pow(base, exponent);
    }
    /* boundary handling */
    if (x[i] * best_parameter[i] < 25.0) {
      result += 5.0 * fabs(si) - si * x[i];
    } else {
      result += 5.0 * fabs(si) - si * best_parameter[i];
    }
  }

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_linear_slope_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_linear_slope_raw(x, problem->number_of_variables, problem->best_parameter);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic linear slope problem.
 */
static coco_problem_t *f_linear_slope_allocate(const size_t number_of_variables, const double *best_parameter) {

  size_t i;
  /* best_parameter will be overwritten below */
  coco_problem_t *problem = coco_problem_allocate_from_scalars("linear slope function",
      f_linear_slope_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "linear_slope", number_of_variables);

  /* Compute best solution */
  for (i = 0; i < number_of_variables; ++i) {
    if (best_parameter[i] < 0.0) {
      problem->best_parameter[i] = problem->smallest_values_of_interest[i];
    } else {
      problem->best_parameter[i] = problem->largest_values_of_interest[i];
    }
  }
  f_linear_slope_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB linear slope problem.
 */
static coco_problem_t *f_linear_slope_bbob_problem_allocate(const size_t function,
                                                            const size_t dimension,
                                                            const size_t instance,
                                                            const long rseed,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  xopt = coco_allocate_vector(dimension);
  bbob2009_compute_xopt(xopt, rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);

  problem = f_linear_slope_allocate(dimension, xopt);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "1-separable");

  coco_free_memory(xopt);
  return problem;
}
#line 20 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_lunacek_bi_rastrigin.c"
/**
 * @file f_lunacek_bi_rastrigin.c
 * @brief Implementation of the Lunacek bi-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/f_lunacek_bi_rastrigin.c"
#line 11 "code-experiments/src/f_lunacek_bi_rastrigin.c"
#line 12 "code-experiments/src/f_lunacek_bi_rastrigin.c"
#line 13 "code-experiments/src/f_lunacek_bi_rastrigin.c"

/**
 * @brief Data type for the Lunacek bi-Rastrigin problem.
 */
typedef struct {
  double *x_hat, *z;
  double *xopt, fopt;
  double **rot1, **rot2;
  long rseed;
  coco_problem_free_function_t old_free_problem;
} f_lunacek_bi_rastrigin_data_t;

/**
 * @brief Implements the Lunacek bi-Rastrigin function without connections to any COCO structures.
 */
static double f_lunacek_bi_rastrigin_raw(const double *x,
                                         const size_t number_of_variables,
                                         f_lunacek_bi_rastrigin_data_t *data) {
  double result;
  static const double condition = 100.;
  size_t i, j;
  double penalty = 0.0;
  static const double mu0 = 2.5;
  static const double d = 1.;
  const double s = 1. - 0.5 / (sqrt((double) (number_of_variables + 20)) - 4.1);
  const double mu1 = -sqrt((mu0 * mu0 - d) / s);
  double *tmpvect, sum1 = 0., sum2 = 0., sum3 = 0.;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double tmp;
    tmp = fabs(x[i]) - 5.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  /* x_hat */
  for (i = 0; i < number_of_variables; ++i) {
    data->x_hat[i] = 2. * x[i];
    if (data->xopt[i] < 0.) {
      data->x_hat[i] *= -1.;
    }
  }

  tmpvect = coco_allocate_vector(number_of_variables);
  /* affine transformation */
  for (i = 0; i < number_of_variables; ++i) {
    double c1;
    tmpvect[i] = 0.0;
    c1 = pow(sqrt(condition), ((double) i) / (double) (number_of_variables - 1));
    for (j = 0; j < number_of_variables; ++j) {
      tmpvect[i] += c1 * data->rot2[i][j] * (data->x_hat[j] - mu0);
    }
  }
  for (i = 0; i < number_of_variables; ++i) {
    data->z[i] = 0;
    for (j = 0; j < number_of_variables; ++j) {
      data->z[i] += data->rot1[i][j] * tmpvect[j];
    }
  }
  /* Computation core */
  for (i = 0; i < number_of_variables; ++i) {
    sum1 += (data->x_hat[i] - mu0) * (data->x_hat[i] - mu0);
    sum2 += (data->x_hat[i] - mu1) * (data->x_hat[i] - mu1);
    sum3 += cos(2 * coco_pi * data->z[i]);
  }
  result = coco_double_min(sum1, d * (double) number_of_variables + s * sum2)
      + 10. * ((double) number_of_variables - sum3) + 1e4 * penalty;
  coco_free_memory(tmpvect);

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_lunacek_bi_rastrigin_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_lunacek_bi_rastrigin_raw(x, problem->number_of_variables, (f_lunacek_bi_rastrigin_data_t *) problem->data);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the Lunacek bi-Rastrigin data object.
 */
static void f_lunacek_bi_rastrigin_free(coco_problem_t *problem) {
  f_lunacek_bi_rastrigin_data_t *data;
  data = (f_lunacek_bi_rastrigin_data_t *) problem->data;
  coco_free_memory(data->x_hat);
  coco_free_memory(data->z);
  coco_free_memory(data->xopt);
  bbob2009_free_matrix(data->rot1, problem->number_of_variables);
  bbob2009_free_matrix(data->rot2, problem->number_of_variables);

  /* Let the generic free problem code deal with all of the
   * coco_problem_t fields.
   */
  problem->problem_free_function = NULL;
  coco_problem_free(problem);
}

/**
 * @brief Creates the BBOB Lunacek bi-Rastrigin problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_t *f_lunacek_bi_rastrigin_bbob_problem_allocate(const size_t function,
                                                                    const size_t dimension,
                                                                    const size_t instance,
                                                                    const long rseed,
                                                                    const char *problem_id_template,
                                                                    const char *problem_name_template) {

  f_lunacek_bi_rastrigin_data_t *data;
  coco_problem_t *problem = coco_problem_allocate_from_scalars("Lunacek's bi-Rastrigin function",
      f_lunacek_bi_rastrigin_evaluate, f_lunacek_bi_rastrigin_free, dimension, -5.0, 5.0, 0.0);

  const double mu0 = 2.5;

  double fopt, *tmpvect;
  size_t i;

  data = (f_lunacek_bi_rastrigin_data_t *) coco_allocate_memory(sizeof(*data));
  /* Allocate temporary storage and space for the rotation matrices */
  data->x_hat = coco_allocate_vector(dimension);
  data->z = coco_allocate_vector(dimension);
  data->xopt = coco_allocate_vector(dimension);
  data->rot1 = bbob2009_allocate_matrix(dimension, dimension);
  data->rot2 = bbob2009_allocate_matrix(dimension, dimension);
  data->rseed = rseed;

  data->fopt = bbob2009_compute_fopt(24, instance);
  bbob2009_compute_xopt(data->xopt, rseed, dimension);
  bbob2009_compute_rotation(data->rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(data->rot2, rseed, dimension);

  problem->data = data;

  /* Compute best solution */
  tmpvect = coco_allocate_vector(dimension);
  bbob2009_gauss(tmpvect, dimension, rseed);
  for (i = 0; i < dimension; ++i) {
    data->xopt[i] = 0.5 * mu0;
    if (tmpvect[i] < 0.0) {
      data->xopt[i] *= -1.0;
    }
    problem->best_parameter[i] = data->xopt[i];
  }
  coco_free_memory(tmpvect);
  f_lunacek_bi_rastrigin_evaluate(problem, problem->best_parameter, problem->best_value);

  fopt = bbob2009_compute_fopt(function, instance);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "5-weakly-structured");

  return problem;
}
#line 21 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_rastrigin.c"
/**
 * @file f_rastrigin.c
 * @brief Implementation of the Rastrigin function and problem.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#line 11 "code-experiments/src/f_rastrigin.c"
#line 12 "code-experiments/src/f_rastrigin.c"
#line 13 "code-experiments/src/f_rastrigin.c"
#line 1 "code-experiments/src/transform_vars_conditioning.c"
/**
 * @file transform_vars_conditioning.c
 * @brief Implementation of conditioning decision values.
 */

#include <math.h>
#include <assert.h>

#line 10 "code-experiments/src/transform_vars_conditioning.c"
#line 11 "code-experiments/src/transform_vars_conditioning.c"

/**
 * @brief Data type for transform_vars_conditioning.
 */
typedef struct {
  double *x;
  double alpha;
} transform_vars_conditioning_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_conditioning_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_conditioning_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_conditioning_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  for (i = 0; i < problem->number_of_variables; ++i) {
    /* OME: We could precalculate the scaling coefficients if we
     * really wanted to.
     */
    data->x[i] = pow(data->alpha, 0.5 * (double) (long) i / ((double) (long) problem->number_of_variables - 1.0))
        * x[i];
  }
  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

static void transform_vars_conditioning_free(void *thing) {
  transform_vars_conditioning_data_t *data = (transform_vars_conditioning_data_t *) thing;
  coco_free_memory(data->x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_conditioning(coco_problem_t *inner_problem, const double alpha) {
  transform_vars_conditioning_data_t *data;
  coco_problem_t *problem;

  data = (transform_vars_conditioning_data_t *) coco_allocate_memory(sizeof(*data));
  data->x = coco_allocate_vector(inner_problem->number_of_variables);
  data->alpha = alpha;
  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_conditioning_free, "transform_vars_conditioning");
  problem->evaluate_function = transform_vars_conditioning_evaluate;

  if (coco_problem_best_parameter_not_zero(inner_problem)) {
    coco_warning("transform_vars_conditioning(): 'best_parameter' not updated, set to NAN");
    coco_vector_set_to_nan(inner_problem->best_parameter, inner_problem->number_of_variables);
  }  return problem;
}
#line 14 "code-experiments/src/f_rastrigin.c"
#line 15 "code-experiments/src/f_rastrigin.c"
#line 16 "code-experiments/src/f_rastrigin.c"
#line 17 "code-experiments/src/f_rastrigin.c"
#line 18 "code-experiments/src/f_rastrigin.c"
#line 19 "code-experiments/src/f_rastrigin.c"

/**
 * @brief Implements the Rastrigin function without connections to any COCO structures.
 */
static double f_rastrigin_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;
  double result;
  double sum1 = 0.0, sum2 = 0.0;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    sum1 += cos(coco_two_pi * x[i]);
    sum2 += x[i] * x[i];
  }
  result = 10.0 * ((double) (long) number_of_variables - sum1) + sum2;

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_rastrigin_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_rastrigin_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Rastrigin problem.
 */
static coco_problem_t *f_rastrigin_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Rastrigin function",
      f_rastrigin_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "rastrigin", number_of_variables);

  /* Compute best solution */
  f_rastrigin_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Rastrigin problem.
 */
static coco_problem_t *f_rastrigin_bbob_problem_allocate(const size_t function,
                                                         const size_t dimension,
                                                         const size_t instance,
                                                         const long rseed,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_conditioning(problem, 10.0);
  problem = transform_vars_asymmetric(problem, 0.2);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "1-separable");

  coco_free_memory(xopt);
  return problem;
}

/**
 * @brief Creates the BBOB rotated Rastrigin problem.
 */
static coco_problem_t *f_rastrigin_rotated_bbob_problem_allocate(const size_t function,
                                                                 const size_t dimension,
                                                                 const size_t instance,
                                                                 const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j, k;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, **rot2;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  rot2 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2, rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      for (k = 0; k < dimension; ++k) {
        double exponent = 1.0 * (int) k / ((double) (long) dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(sqrt(10), exponent) * rot2[k][j];
      }
    }
  }

  problem = f_rastrigin_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_asymmetric(problem, 0.2);
  problem = transform_vars_oscillate(problem);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);

  bbob2009_free_matrix(rot1, dimension);
  bbob2009_free_matrix(rot2, dimension);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}

#line 22 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_rosenbrock.c"
/**
 * @file f_rosenbrock.c
 * @brief Implementation of the Rosenbrock function and problem.
 */

#include <assert.h>

#line 9 "code-experiments/src/f_rosenbrock.c"
#line 10 "code-experiments/src/f_rosenbrock.c"
#line 11 "code-experiments/src/f_rosenbrock.c"
#line 12 "code-experiments/src/f_rosenbrock.c"
#line 1 "code-experiments/src/transform_vars_scale.c"
/**
 * @file transform_vars_scale.c
 * @brief Implementation of scaling decision values by a given factor.
 */

#include <assert.h>

#line 9 "code-experiments/src/transform_vars_scale.c"
#line 10 "code-experiments/src/transform_vars_scale.c"

/**
 * @brief Data type for transform_vars_scale.
 */
typedef struct {
  double factor;
  double *x;
} transform_vars_scale_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_scale_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_scale_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_scale_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);
  do {
    const double factor = data->factor;

    for (i = 0; i < problem->number_of_variables; ++i) {
      data->x[i] = factor * x[i];
    }
    coco_evaluate_function(inner_problem, data->x, y);
    assert(y[0] + 1e-13 >= problem->best_value[0]);
  } while (0);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_scale_free(void *thing) {
  transform_vars_scale_data_t *data = (transform_vars_scale_data_t *) thing;
  coco_free_memory(data->x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_scale(coco_problem_t *inner_problem, const double factor) {
  transform_vars_scale_data_t *data;
  coco_problem_t *problem;
  size_t i;
  data = (transform_vars_scale_data_t *) coco_allocate_memory(sizeof(*data));
  data->factor = factor;
  data->x = coco_allocate_vector(inner_problem->number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_scale_free, "transform_vars_scale");
  problem->evaluate_function = transform_vars_scale_evaluate;
  /* Compute best parameter */
  if (data->factor != 0.) {
      for (i = 0; i < problem->number_of_variables; i++) {
          problem->best_parameter[i] /= data->factor;
      }
  } /* else error? */
  return problem;
}
#line 13 "code-experiments/src/f_rosenbrock.c"
#line 14 "code-experiments/src/f_rosenbrock.c"
#line 15 "code-experiments/src/f_rosenbrock.c"

/**
 * @brief Implements the Rosenbrock function without connections to any COCO structures.
 */
static double f_rosenbrock_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;
  double result;
  double s1 = 0.0, s2 = 0.0, tmp;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  for (i = 0; i < number_of_variables - 1; ++i) {
    tmp = (x[i] * x[i] - x[i + 1]);
    s1 += tmp * tmp;
    tmp = (x[i] - 1.0);
    s2 += tmp * tmp;
  }
  result = 100.0 * s1 + s2;

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_rosenbrock_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_rosenbrock_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Rosenbrock problem.
 */
static coco_problem_t *f_rosenbrock_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Rosenbrock function",
      f_rosenbrock_evaluate, NULL, number_of_variables, -5.0, 5.0, 1.0);
  coco_problem_set_id(problem, "%s_d%02lu", "rosenbrock", number_of_variables);

  /* Compute best solution */
  f_rosenbrock_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Rosenbrock problem.
 */
static coco_problem_t *f_rosenbrock_bbob_problem_allocate(const size_t function,
                                                          const size_t dimension,
                                                          const size_t instance,
                                                          const long rseed,
                                                          const char *problem_id_template,
                                                          const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i;
  double *minus_one, factor;

  minus_one = coco_allocate_vector(dimension);
  xopt = coco_allocate_vector(dimension);
  bbob2009_compute_xopt(xopt, rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    minus_one[i] = -1.0;
    xopt[i] *= 0.75;
  }
  fopt = bbob2009_compute_fopt(function, instance);
  factor = coco_double_max(1.0, sqrt((double) dimension) / 8.0);

  problem = f_rosenbrock_allocate(dimension);
  problem = transform_vars_shift(problem, minus_one, 0);
  problem = transform_vars_scale(problem, factor);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "2-moderate");

  coco_free_memory(minus_one);
  coco_free_memory(xopt);
  return problem;
}

/**
 * @brief Creates the BBOB rotated Rosenbrock problem.
 */
static coco_problem_t *f_rosenbrock_rotated_bbob_problem_allocate(const size_t function,
                                                                  const size_t dimension,
                                                                  const size_t instance,
                                                                  const long rseed,
                                                                  const char *problem_id_template,
                                                                  const char *problem_name_template) {

  double fopt;
  coco_problem_t *problem = NULL;
  size_t row, column;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, factor;

  fopt = bbob2009_compute_fopt(function, instance);
  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed, dimension);

  factor = coco_double_max(1.0, sqrt((double) dimension) / 8.0);
  /* Compute affine transformation */
  for (row = 0; row < dimension; ++row) {
    current_row = M + row * dimension;
    for (column = 0; column < dimension; ++column) {
      current_row[column] = factor * rot1[row][column];
    }
    b[row] = 0.5;
  }
  bbob2009_free_matrix(rot1, dimension);

  problem = f_rosenbrock_allocate(dimension);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "2-moderate");

  coco_free_memory(M);
  coco_free_memory(b);
  return problem;
}
#line 23 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_schaffers.c"
/**
 * @file f_schaffers.c
 * @brief Implementation of the Schaffer's F7 function and problem, transformations not implemented for the
 * moment.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#line 12 "code-experiments/src/f_schaffers.c"
#line 13 "code-experiments/src/f_schaffers.c"
#line 14 "code-experiments/src/f_schaffers.c"
#line 15 "code-experiments/src/f_schaffers.c"
#line 16 "code-experiments/src/f_schaffers.c"
#line 17 "code-experiments/src/f_schaffers.c"
#line 18 "code-experiments/src/f_schaffers.c"
#line 19 "code-experiments/src/f_schaffers.c"

/**
 * @brief Implements the Schaffer's F7 function without connections to any COCO structures.
 */
static double f_schaffers_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;
  double result;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  /* Computation core */
  result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double tmp = x[i] * x[i] + x[i + 1] * x[i + 1];
    result += pow(tmp, 0.25) * (1.0 + pow(sin(50.0 * pow(tmp, 0.1)), 2.0));
  }
  result = pow(result / ((double) (long) number_of_variables - 1.0), 2.0);

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_schaffers_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_schaffers_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Schaffer's F7 problem.
 */
static coco_problem_t *f_schaffers_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Schaffer's function",
      f_schaffers_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "schaffers", number_of_variables);

  /* Compute best solution */
  f_schaffers_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Schaffer's F7 problem.
 */
static coco_problem_t *f_schaffers_bbob_problem_allocate(const size_t function,
                                                         const size_t dimension,
                                                         const size_t instance,
                                                         const long rseed,
                                                         const double conditioning,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, **rot2;

  const double penalty_factor = 10.0;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  rot2 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2, rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      double exponent = 1.0 * (int) i / ((double) (long) dimension - 1.0);
      current_row[j] = rot2[i][j] * pow(sqrt(conditioning), exponent);
    }
  }

  problem = f_schaffers_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_asymmetric(problem, 0.5);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_penalize(problem, penalty_factor);

  bbob2009_free_matrix(rot1, dimension);
  bbob2009_free_matrix(rot2, dimension);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 24 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_schwefel.c"
/**
 * @file f_schwefel.c
 * @brief Implementation of the Schwefel function and problem.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#line 11 "code-experiments/src/f_schwefel.c"
#line 12 "code-experiments/src/f_schwefel.c"
#line 13 "code-experiments/src/f_schwefel.c"
#line 14 "code-experiments/src/f_schwefel.c"
#line 15 "code-experiments/src/f_schwefel.c"
#line 16 "code-experiments/src/f_schwefel.c"
#line 17 "code-experiments/src/f_schwefel.c"
#line 1 "code-experiments/src/transform_vars_z_hat.c"
/**
 * @file transform_vars_z_hat.c
 * @brief Implementation of the z^hat transformation of decision values for the BBOB Schwefel problem.
 */

#include <assert.h>

#line 9 "code-experiments/src/transform_vars_z_hat.c"
#line 10 "code-experiments/src/transform_vars_z_hat.c"

/**
 * @brief Data type for transform_vars_z_hat.
 */
typedef struct {
  double *xopt;
  double *z;
  coco_problem_free_function_t old_free_problem;
} transform_vars_z_hat_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_z_hat_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_z_hat_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

  data = (transform_vars_z_hat_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  data->z[0] = x[0];

  for (i = 1; i < problem->number_of_variables; ++i) {
    data->z[i] = x[i] + 0.25 * (x[i - 1] - 2.0 * fabs(data->xopt[i - 1]));
  }
  coco_evaluate_function(inner_problem, data->z, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_z_hat_free(void *thing) {
  transform_vars_z_hat_data_t *data = (transform_vars_z_hat_data_t *) thing;
  coco_free_memory(data->xopt);
  coco_free_memory(data->z);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_z_hat(coco_problem_t *inner_problem, const double *xopt) {
  transform_vars_z_hat_data_t *data;
  coco_problem_t *problem;
  data = (transform_vars_z_hat_data_t *) coco_allocate_memory(sizeof(*data));
  data->xopt = coco_duplicate_vector(xopt, inner_problem->number_of_variables);
  data->z = coco_allocate_vector(inner_problem->number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_z_hat_free, "transform_vars_z_hat");
  problem->evaluate_function = transform_vars_z_hat_evaluate;
  /* TODO: When should this warning be output?
   coco_warning("transform_vars_z_hat(): 'best_parameter' not updated"); */
  return problem;
}
#line 18 "code-experiments/src/f_schwefel.c"
#line 1 "code-experiments/src/transform_vars_x_hat.c"
/**
 * @file transform_vars_x_hat.c
 * @brief Implementation of multiplying the decision values by the vector 1+-.
 */

#include <assert.h>

#line 9 "code-experiments/src/transform_vars_x_hat.c"
#line 10 "code-experiments/src/transform_vars_x_hat.c"
#line 11 "code-experiments/src/transform_vars_x_hat.c"

/**
 * @brief Data type for transform_vars_x_hat.
 */
typedef struct {
  long seed;
  double *x;
  coco_problem_free_function_t old_free_problem;
} transform_vars_x_hat_data_t;

/**
 * @brief Evaluates the transformation.
 */
static void transform_vars_x_hat_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_x_hat_data_t *data;
  coco_problem_t *inner_problem;

  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
  	coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
  	return;
  }

 data = (transform_vars_x_hat_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);
  do {
    bbob2009_unif(data->x, problem->number_of_variables, data->seed);

    for (i = 0; i < problem->number_of_variables; ++i) {
      if (data->x[i] - 0.5 < 0.0) {
        data->x[i] = -x[i];
      } else {
        data->x[i] = x[i];
      }
    }
    coco_evaluate_function(inner_problem, data->x, y);
    assert(y[0] + 1e-13 >= problem->best_value[0]);
  } while (0);
}

/**
 * @brief Frees the data object.
 */
static void transform_vars_x_hat_free(void *thing) {
  transform_vars_x_hat_data_t *data = (transform_vars_x_hat_data_t *) thing;
  coco_free_memory(data->x);
}

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_x_hat(coco_problem_t *inner_problem, const long seed) {
  transform_vars_x_hat_data_t *data;
  coco_problem_t *problem;
  const char *result;
  size_t i;

  data = (transform_vars_x_hat_data_t *) coco_allocate_memory(sizeof(*data));
  data->seed = seed;
  data->x = coco_allocate_vector(inner_problem->number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_x_hat_free, "transform_vars_x_hat");
  problem->evaluate_function = transform_vars_x_hat_evaluate;
  /* Dirty way of setting the best parameter of the transformed f_schwefel... */
  bbob2009_unif(data->x, problem->number_of_variables, data->seed);
  result = strstr(coco_problem_get_id(inner_problem), "schwefel");
	if (result != NULL) {
		for (i = 0; i < problem->number_of_variables; ++i) {
			if (data->x[i] - 0.5 < 0.0) {
				problem->best_parameter[i] = -0.5 * 4.2096874633;
			} else {
				problem->best_parameter[i] = 0.5 * 4.2096874633;
			}
		}
	} else if (coco_problem_best_parameter_not_zero(inner_problem)) {
		coco_warning("transform_vars_x_hat(): 'best_parameter' not updated, set to NAN");
		coco_vector_set_to_nan(inner_problem->best_parameter, inner_problem->number_of_variables);
	}
	return problem;
}
#line 19 "code-experiments/src/f_schwefel.c"

/**
 * @brief Implements the Schwefel function without connections to any COCO structures.
 */
static double f_schwefel_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;
  double result;
  double penalty, sum;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  /* Boundary handling*/
  penalty = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    const double tmp = fabs(x[i]) - 500.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  /* Computation core */
  sum = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    sum += x[i] * sin(sqrt(fabs(x[i])));
  }
  result = 0.01 * (penalty + 418.9828872724339 - sum / (double) number_of_variables);

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_schwefel_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_schwefel_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Schwefel problem.
 */
static coco_problem_t *f_schwefel_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("Schwefel function",
      f_schwefel_evaluate, NULL, number_of_variables, -5.0, 5.0, 420.96874633);
  coco_problem_set_id(problem, "%s_d%02lu", "schwefel", number_of_variables);

  /* Compute best solution: best_parameter[i] = 200 * fabs(xopt[i]) */
  f_schwefel_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Schwefel problem.
 */
static coco_problem_t *f_schwefel_bbob_problem_allocate(const size_t function,
                                                        const size_t dimension,
                                                        const size_t instance,
                                                        const long rseed,
                                                        const char *problem_id_template,
                                                        const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j;

  const double condition = 10.;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row;

  double *tmp1 = coco_allocate_vector(dimension);
  double *tmp2 = coco_allocate_vector(dimension);

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_unif(tmp1, dimension, rseed);
  for (i = 0; i < dimension; ++i) {
    xopt[i] = 0.5 * 4.2096874637;
    if (tmp1[i] - 0.5 < 0) {
      xopt[i] *= -1;
    }
  }

  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      if (i == j) {
        double exponent = 1.0 * (int) i / ((double) (long) dimension - 1);
        current_row[j] = pow(sqrt(condition), exponent);
      }
    }
  }

  for (i = 0; i < dimension; ++i) {
    tmp1[i] = -2 * fabs(xopt[i]);
    tmp2[i] = 2 * fabs(xopt[i]);
  }

  problem = f_schwefel_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_scale(problem, 100);
  problem = transform_vars_shift(problem, tmp1, 0);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, tmp2, 0);
  problem = transform_vars_z_hat(problem, xopt);
  problem = transform_vars_scale(problem, 2);
  problem = transform_vars_x_hat(problem, rseed);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "5-weakly-structured");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(tmp1);
  coco_free_memory(tmp2);
  coco_free_memory(xopt);
  return problem;
}
#line 25 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_sharp_ridge.c"
/**
 * @file f_sharp_ridge.c
 * @brief Implementation of the sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#line 10 "code-experiments/src/f_sharp_ridge.c"
#line 11 "code-experiments/src/f_sharp_ridge.c"
#line 12 "code-experiments/src/f_sharp_ridge.c"
#line 13 "code-experiments/src/f_sharp_ridge.c"
#line 14 "code-experiments/src/f_sharp_ridge.c"
#line 15 "code-experiments/src/f_sharp_ridge.c"

/**
 * @brief Implements the sharp ridge function without connections to any COCO structures.
 */
static double f_sharp_ridge_raw(const double *x, const size_t number_of_variables) {

  static const double alpha = 100.0;
  const double vars_40 = number_of_variables <= 40 ? 1 : number_of_variables / 40.0;
  size_t i = 0;
  double result;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = 0.0;
  for (i = ceil(vars_40); i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result / vars_40);
  for (i = 0; i < ceil(vars_40); ++i)
    result += x[i] * x[i] / vars_40;

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_sharp_ridge_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_sharp_ridge_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic sharp ridge problem.
 */
static coco_problem_t *f_sharp_ridge_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("sharp ridge function",
      f_sharp_ridge_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "sharp_ridge", number_of_variables);

  /* Compute best solution */
  f_sharp_ridge_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB sharp ridge problem.
 */
static coco_problem_t *f_sharp_ridge_bbob_problem_allocate(const size_t function,
                                                           const size_t dimension,
                                                           const size_t instance,
                                                           const long rseed,
                                                           const char *problem_id_template,
                                                           const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j, k;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, **rot2;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  rot2 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2, rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      for (k = 0; k < dimension; ++k) {
        double exponent = 1.0 * (int) k / ((double) (long) dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(sqrt(10), exponent) * rot2[k][j];
      }
    }
  }
  bbob2009_free_matrix(rot1, dimension);
  bbob2009_free_matrix(rot2, dimension);
  problem = f_sharp_ridge_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
#line 26 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_sphere.c"
/**
 * @file f_sphere.c
 * @brief Implementation of the sphere function and problem.
 */

#include <stdio.h>
#include <assert.h>

#line 10 "code-experiments/src/f_sphere.c"
#line 11 "code-experiments/src/f_sphere.c"
#line 12 "code-experiments/src/f_sphere.c"
#line 13 "code-experiments/src/f_sphere.c"
#line 14 "code-experiments/src/f_sphere.c"

/**
 * @brief Implements the sphere function without connections to any COCO structures.
 */
static double f_sphere_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_sphere_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_sphere_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic sphere problem.
 */
static coco_problem_t *f_sphere_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("sphere function",
      f_sphere_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "sphere", number_of_variables);

  /* Compute best solution */
  f_sphere_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}

/**
 * @brief Creates the BBOB sphere problem.
 */
static coco_problem_t *f_sphere_bbob_problem_allocate(const size_t function,
                                                      const size_t dimension,
                                                      const size_t instance,
                                                      const long rseed,
                                                      const char *problem_id_template,
                                                      const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  xopt = coco_allocate_vector(dimension);
  bbob2009_compute_xopt(xopt, rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);

  problem = f_sphere_allocate(dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "1-separable");

  coco_free_memory(xopt);
  return problem;
}

#line 27 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_step_ellipsoid.c"
/**
 * @file f_step_ellipsoid.c
 * @brief Implementation of the step ellipsoid function and problem.
 *
 * The BBOB step ellipsoid function intertwines the variable and objective transformations in such a way
 * that it is hard to devise a composition of generic transformations to implement it. In the end one would
 * have to implement several custom transformations which would be used solely by this problem. Therefore
 * we opt to implement it as a monolithic function instead.
 *
 * TODO: It would be nice to have a generic step ellipsoid function to complement this one.
 */
#include <assert.h>

#line 15 "code-experiments/src/f_step_ellipsoid.c"
#line 16 "code-experiments/src/f_step_ellipsoid.c"
#line 17 "code-experiments/src/f_step_ellipsoid.c"
#line 18 "code-experiments/src/f_step_ellipsoid.c"

/**
 * @brief Data type for the step ellipsoid problem.
 */
typedef struct {
  double *x, *xx;
  double *xopt, fopt;
  double **rot1, **rot2;
} f_step_ellipsoid_data_t;

/**
 * @brief Implements the step ellipsoid function without connections to any COCO structures.
 */
static double f_step_ellipsoid_raw(const double *x, const size_t number_of_variables, f_step_ellipsoid_data_t *data) {

  static const double condition = 100;
  static const double alpha = 10.0;
  size_t i, j;
  double penalty = 0.0, x1;
  double result;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double tmp;
    tmp = fabs(x[i]) - 5.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  for (i = 0; i < number_of_variables; ++i) {
    double c1;
    data->x[i] = 0.0;
    c1 = sqrt(pow(condition / 10., (double) i / (double) (number_of_variables - 1)));
    for (j = 0; j < number_of_variables; ++j) {
      data->x[i] += c1 * data->rot2[i][j] * (x[j] - data->xopt[j]);
    }
  }
  x1 = data->x[0];

  for (i = 0; i < number_of_variables; ++i) {
    if (fabs(data->x[i]) > 0.5)
      data->x[i] = coco_double_round(data->x[i]);
    else
      data->x[i] = coco_double_round(alpha * data->x[i]) / alpha;
  }

  for (i = 0; i < number_of_variables; ++i) {
    data->xx[i] = 0.0;
    for (j = 0; j < number_of_variables; ++j) {
      data->xx[i] += data->rot1[i][j] * data->x[j];
    }
  }

  /* Computation core */
  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    double exponent;
    exponent = (double) (long) i / ((double) (long) number_of_variables - 1.0);
    result += pow(condition, exponent) * data->xx[i] * data->xx[i];
    ;
  }
  result = 0.1 * coco_double_max(fabs(x1) * 1.0e-4, result) + penalty + data->fopt;

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_step_ellipsoid_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_step_ellipsoid_raw(x, problem->number_of_variables, (f_step_ellipsoid_data_t *) problem->data);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Frees the step ellipsoid data object.
 */
static void f_step_ellipsoid_free(coco_problem_t *problem) {
  f_step_ellipsoid_data_t *data;
  data = (f_step_ellipsoid_data_t *) problem->data;
  coco_free_memory(data->x);
  coco_free_memory(data->xx);
  coco_free_memory(data->xopt);
  bbob2009_free_matrix(data->rot1, problem->number_of_variables);
  bbob2009_free_matrix(data->rot2, problem->number_of_variables);
  /* Let the generic free problem code deal with all of the coco_problem_t fields */
  problem->problem_free_function = NULL;
  coco_problem_free(problem);
}

/**
 * @brief Creates the BBOB step ellipsoid problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_t *f_step_ellipsoid_bbob_problem_allocate(const size_t function,
                                                              const size_t dimension,
                                                              const size_t instance,
                                                              const long rseed,
                                                              const char *problem_id_template,
                                                              const char *problem_name_template) {

  f_step_ellipsoid_data_t *data;
  size_t i;
  coco_problem_t *problem = coco_problem_allocate_from_scalars("step ellipsoid function",
      f_step_ellipsoid_evaluate, f_step_ellipsoid_free, dimension, -5.0, 5.0, 0);

  data = (f_step_ellipsoid_data_t *) coco_allocate_memory(sizeof(*data));
  /* Allocate temporary storage and space for the rotation matrices */
  data->x = coco_allocate_vector(dimension);
  data->xx = coco_allocate_vector(dimension);
  data->xopt = coco_allocate_vector(dimension);
  data->rot1 = bbob2009_allocate_matrix(dimension, dimension);
  data->rot2 = bbob2009_allocate_matrix(dimension, dimension);

  data->fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(data->xopt, rseed, dimension);
  bbob2009_compute_rotation(data->rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(data->rot2, rseed, dimension);

  problem->data = data;
  
  /* Compute best solution
   *
   * OME: Dirty hack for now because I did not want to invert the
   * transformations to find the best_parameter :/
   */
  for (i = 0; i < problem->number_of_variables; i++) {
      problem->best_parameter[i] = data->xopt[i];
  }
  problem->best_value[0] = data->fopt;

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "2-moderate");

  return problem;
}
#line 28 "code-experiments/src/suite_bbob.c"
#line 1 "code-experiments/src/f_weierstrass.c"
/**
 * @file f_weierstrass.c
 * @brief Implementation of the Weierstrass function and problem.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#line 11 "code-experiments/src/f_weierstrass.c"
#line 12 "code-experiments/src/f_weierstrass.c"
#line 13 "code-experiments/src/f_weierstrass.c"
#line 14 "code-experiments/src/f_weierstrass.c"
#line 15 "code-experiments/src/f_weierstrass.c"
#line 16 "code-experiments/src/f_weierstrass.c"
#line 17 "code-experiments/src/f_weierstrass.c"
#line 18 "code-experiments/src/f_weierstrass.c"

/** @brief Number of summands in the Weierstrass problem. */
#define F_WEIERSTRASS_SUMMANDS 12

/**
 * @brief Data type for the Weierstrass problem.
 */
typedef struct {
  double f0;
  double ak[F_WEIERSTRASS_SUMMANDS];
  double bk[F_WEIERSTRASS_SUMMANDS];
} f_weierstrass_data_t;

/**
 * @brief Implements the Weierstrass function without connections to any COCO structures.
 */
static double f_weierstrass_raw(const double *x, const size_t number_of_variables, f_weierstrass_data_t *data) {

  size_t i, j;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
  	return NAN;

  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    for (j = 0; j < F_WEIERSTRASS_SUMMANDS; ++j) {
      result += cos(2 * coco_pi * (x[i] + 0.5) * data->bk[j]) * data->ak[j];
    }
  }
  result = 10.0 * pow(result / (double) (long) number_of_variables - data->f0, 3.0);

  return result;
}

/**
 * @brief Uses the raw function to evaluate the COCO problem.
 */
static void f_weierstrass_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_weierstrass_raw(x, problem->number_of_variables, (f_weierstrass_data_t *) problem->data);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic Weierstrass problem.
 */
static coco_problem_t *f_weierstrass_allocate(const size_t number_of_variables) {

  f_weierstrass_data_t *data;
  size_t i;
  double *non_unique_best_value;
  coco_problem_t *problem = coco_problem_allocate_from_scalars("Weierstrass function",
      f_weierstrass_evaluate, NULL, number_of_variables, -5.0, 5.0, NAN);
  coco_problem_set_id(problem, "%s_d%02lu", "weierstrass", number_of_variables);

  data = (f_weierstrass_data_t *) coco_allocate_memory(sizeof(*data));
  data->f0 = 0.0;
  for (i = 0; i < F_WEIERSTRASS_SUMMANDS; ++i) {
    data->ak[i] = pow(0.5, (double) i);
    data->bk[i] = pow(3., (double) i);
    data->f0 += data->ak[i] * cos(2 * coco_pi * data->bk[i] * 0.5);
  }
  problem->data = data;

  /* Compute best solution */
  non_unique_best_value = coco_allocate_vector(number_of_variables);
  for (i = 0; i < number_of_variables; i++)
    non_unique_best_value[i] = 0.0;
  f_weierstrass_evaluate(problem, non_unique_best_value, problem->best_value);
  coco_free_memory(non_unique_best_value);
  return problem;
}

/**
 * @brief Creates the BBOB Weierstrass problem.
 */
static coco_problem_t *f_weierstrass_bbob_problem_allocate(const size_t function,
                                                           const size_t dimension,
                                                           const size_t instance,
                                                           const long rseed,
                                                           const char *problem_id_template,
                                                           const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  size_t i, j, k;
  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double *current_row, **rot1, **rot2;

  const double condition = 100.0;
  const double penalty_factor = 10.0 / (double) dimension;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_compute_xopt(xopt, rseed, dimension);

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  rot2 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2, rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      for (k = 0; k < dimension; ++k) {
        const double base = 1.0 / sqrt(condition);
        const double exponent = 1.0 * (int) k / ((double) (long) dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(base, exponent) * rot2[k][j];
      }
    }
  }

  problem = f_weierstrass_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_oscillate(problem);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_penalize(problem, penalty_factor);

  bbob2009_free_matrix(rot1, dimension);
  bbob2009_free_matrix(rot2, dimension);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}

#undef F_WEIERSTRASS_SUMMANDS
#line 29 "code-experiments/src/suite_bbob.c"

static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                         const size_t number_of_functions,
                                         const size_t number_of_dimensions,
                                         const size_t *dimensions,
                                         const char *default_instances);

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static coco_suite_t *suite_bbob_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = { 2, 3, 5, 10, 20, 40 };

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = coco_suite_allocate("bbob", 24, 6, dimensions, "year: 2016");

  return suite;
}

/**
 * @brief Sets the instances associated with years for the bbob suite.
 */
static const char *suite_bbob_get_instances_by_year(const int year) {

  if (year == 2009) {
    return "1-5,1-5,1-5";
  }
  else if (year == 2010) {
    return "1-15";
  }
  else if (year == 2012) {
    return "1-5,21-30";
  }
  else if (year == 2013) {
    return "1-5,31-40";
  }
  else if (year == 2015) {
    return "1-5,41-50";
  }
  else if (year == 2016) {
    return "1-5,51-60";
  }
  else {
    coco_error("suite_bbob_get_instances_by_year(): year %d not defined for suite_bbob", year);
    return NULL;
  }
}

/**
 * @brief Creates and returns a BBOB problem without needing the actual bbob suite.
 *
 * Useful for other suites as well (see for example suite_biobj.c).
 */
static coco_problem_t *coco_get_bbob_problem(const size_t function,
                                             const size_t dimension,
                                             const size_t instance) {
  coco_problem_t *problem = NULL;

  const char *problem_id_template = "bbob_f%03lu_i%02lu_d%02lu";
  const char *problem_name_template = "BBOB suite problem f%lu instance %lu in %luD";

  const long rseed = (long) (function + 10000 * instance);
  const long rseed_3 = (long) (3 + 10000 * instance);
  const long rseed_17 = (long) (17 + 10000 * instance);

  if (function == 1) {
    problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 2) {
    problem = f_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 3) {
    problem = f_rastrigin_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 4) {
    problem = f_bueche_rastrigin_bbob_problem_allocate(function, dimension, instance, rseed_3,
        problem_id_template, problem_name_template);
  } else if (function == 5) {
    problem = f_linear_slope_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 6) {
    problem = f_attractive_sector_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 7) {
    problem = f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 8) {
    problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 9) {
    problem = f_rosenbrock_rotated_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 10) {
    problem = f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 11) {
    problem = f_discus_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 12) {
    problem = f_bent_cigar_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 13) {
    problem = f_sharp_ridge_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 14) {
    problem = f_different_powers_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 15) {
    problem = f_rastrigin_rotated_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 16) {
    problem = f_weierstrass_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 17) {
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, 10,
        problem_id_template, problem_name_template);
  } else if (function == 18) {
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed_17, 1000,
        problem_id_template, problem_name_template);
  } else if (function == 19) {
    problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 20) {
    problem = f_schwefel_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 21) {
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, 101,
        problem_id_template, problem_name_template);
  } else if (function == 22) {
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, 21,
        problem_id_template, problem_name_template);
  } else if (function == 23) {
    problem = f_katsuura_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else if (function == 24) {
    problem = f_lunacek_bi_rastrigin_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else {
    coco_error("get_bbob_problem(): cannot retrieve problem f%lu instance %lu in %luD",
    		(unsigned long) function, (unsigned long) instance, (unsigned long) dimension);
    return NULL; /* Never reached */
  }

  return problem;
}

/**
 * @brief Returns the problem from the bbob suite that corresponds to the given parameters.
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_t *suite_bbob_get_problem(coco_suite_t *suite,
                                              const size_t function_idx,
                                              const size_t dimension_idx,
                                              const size_t instance_idx) {

  coco_problem_t *problem = NULL;

  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];

  problem = coco_get_bbob_problem(function, dimension, instance);

  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(suite, function_idx, dimension_idx, instance_idx);

  return problem;
}
#line 19 "code-experiments/src/coco_suite.c"
#line 1 "code-experiments/src/suite_biobj.c"
/**
 * @file suite_biobj.c
 * @brief Implementation of the bbob-biobj suite containing 55 functions and 6 dimensions.
 *
 * The bi-objective suite was created by combining two single-objective problems from the bbob suite.
 *
 * @note Because some bi-objective problems constructed from two single-objective ones have a single optimal
 * value, some care must be taken when selecting the instances. The already verified instances are stored in
 * suite_biobj_instances. If a new instance of the problem is called, a check ensures that the two underlying
 * single-objective instances create a true bi-objective problem. However, these new instances need to be
 * manually added to suite_biobj_instances, otherwise they will be computed each time the suite constructor
 * is invoked with these instances.
 */

#line 16 "code-experiments/src/suite_biobj.c"
#line 1 "code-experiments/src/mo_utilities.c"
/**
 * @file mo_utilities.c
 * @brief Definitions of miscellaneous functions used for multi-objective problems.
 */

#include <stdlib.h>
#include <stdio.h>
#line 9 "code-experiments/src/mo_utilities.c"

/**
 * @brief Precision used when comparing multi-objective solutions.
 *
 * Two solutions are considered equal in objective space when their normalized difference is smaller than
 * mo_precision.
 *
 * @note mo_precision needs to be smaller than mo_discretization
 */
static const double mo_precision = 1e-13;

/**
 * @brief Discretization interval used for rounding normalized multi-objective solutions.
 *
 * @note mo_discretization needs to be larger than mo_precision
 */
static const double mo_discretization = 5 * 1e-13;

/**
 * @brief Computes and returns the Euclidean norm of two dim-dimensional points first and second.
 */
static double mo_get_norm(const double *first, const double *second, const size_t dim) {

  size_t i;
  double norm = 0;

  for (i = 0; i < dim; i++) {
    norm += pow(first[i] - second[i], 2);
  }

  return sqrt(norm);
}

/**
 * @brief Creates a rounded normalized version of the given solution w.r.t. the given ROI.
 *
 * If the solution seems to be better than the extremes it is corrected (2 objectives are assumed).
 * The caller is responsible for freeing the allocated memory using coco_free_memory().
 */
static double *mo_normalize(const double *y, const double *ideal, const double *nadir, const size_t num_obj) {

  size_t i;
  double *normalized_y = coco_allocate_vector(num_obj);

  for (i = 0; i < num_obj; i++) {
    assert((nadir[i] - ideal[i]) > mo_discretization);
    normalized_y[i] = (y[i] - ideal[i]) / (nadir[i] - ideal[i]);
    normalized_y[i] = coco_double_round(normalized_y[i] / mo_discretization) * mo_discretization;
    if (normalized_y[i] < 0) {
      coco_debug("Adjusting %.15e to %.15e", y[i], ideal[i]);
      normalized_y[i] = 0;
    }
  }

  for (i = 0; i < num_obj; i++) {
    assert(num_obj == 2);
    if (coco_double_almost_equal(normalized_y[i], 0, mo_precision) && (normalized_y[1-i] < 1)) {
      coco_debug("Adjusting %.15e to %.15e", y[1-i], nadir[1-i]);
      normalized_y[1-i] = 1;
    }
  }

  return normalized_y;
}

/**
 * @brief Checks the dominance relation in the unconstrained minimization case between two normalized
 * solutions in the objective space.
 *
 * If two values are closer together than mo_precision, they are treated as equal.
 *
 * @return
 *  1 if normalized_y1 dominates normalized_y2 <br>
 *  0 if normalized_y1 and normalized_y2 are non-dominated <br>
 * -1 if normalized_y2 dominates normalized_y1 <br>
 * -2 if normalized_y1 is identical to normalized_y2
 */
static int mo_get_dominance(const double *normalized_y1, const double *normalized_y2, const size_t num_obj) {

  size_t i;
  int flag1 = 0;
  int flag2 = 0;

  for (i = 0; i < num_obj; i++) {
    if (coco_double_almost_equal(normalized_y1[i], normalized_y2[i], mo_precision)) {
      continue;
    } else if (normalized_y1[i] < normalized_y2[i]) {
      flag1 = 1;
    } else if (normalized_y1[i] > normalized_y2[i]) {
      flag2 = 1;
    }
  }

  if (flag1 && !flag2) {
    return 1;
  } else if (!flag1 && flag2) {
    return -1;
  } else if (flag1 && flag2) {
    return 0;
  } else { /* (!flag1 && !flag2) */
    return -2;
  }
}

/**
 * @brief Checks whether the normalized solution is within [0, 1]^num_obj.
 */
static int mo_is_within_ROI(const double *normalized_y, const size_t num_obj) {

  size_t i;
  int within = 1;

  for (i = 0; i < num_obj; i++) {
    if (coco_double_almost_equal(normalized_y[i], 0, mo_precision) ||
        coco_double_almost_equal(normalized_y[i], 1, mo_precision) ||
        (normalized_y[i] > 0 && normalized_y[i] < 1))
      continue;
    else
      within = 0;
  }
  return within;
}

/**
 * @brief Computes and returns the minimal normalized distance of the point normalized_y from the ROI
 * (equals 0 if within the ROI).
 *
 *  @note Assumes num_obj = 2 and normalized_y >= 0
 */
static double mo_get_distance_to_ROI(const double *normalized_y, const size_t num_obj) {

  double diff_0, diff_1;

  if (mo_is_within_ROI(normalized_y, num_obj))
    return 0;

  assert(num_obj == 2);
  assert(normalized_y[0] >= 0);
  assert(normalized_y[1] >= 0);

  diff_0 = normalized_y[0] - 1;
  diff_1 = normalized_y[1] - 1;
  if ((diff_0 > 0) && (diff_1 > 0)) {
    return sqrt(pow(diff_0, 2) + pow(diff_1, 2));
  }
  else if (diff_0 > 0)
    return diff_0;
  else
    return diff_1;
}
#line 17 "code-experiments/src/suite_biobj.c"
#line 18 "code-experiments/src/suite_biobj.c"
#line 1 "code-experiments/src/suite_biobj_best_values_hyp.c"
/**
 * @file suite_biobj_best_values_hyp.c
 * @brief Contains the best known hypervolume values for the bbob-biobj suite problems.
 */

/**
 * @brief The best known hypervolume values for the bbob-biobj suite problems.
 *
 * @note Because this file is used for automatically retrieving the existing best hypervolume values for
 * pre-processing purposes, its formatting should not be altered. This means that there must be exactly one
 * string per line, the first string appearing on the next line after "static const char..." (no comments 
 * allowed in between). Nothing should be placed on the last line (line with };).
 */
static const char *suite_biobj_best_values_hyp[] = { /* Best values on 29.03.2016 16:29:20 */
  "bbob-biobj_f01_i01_d02 0.833332421874055",
  "bbob-biobj_f01_i01_d03 0.833330861151064",
  "bbob-biobj_f01_i01_d05 0.833326252335351",
  "bbob-biobj_f01_i01_d10 0.833329281931691",
  "bbob-biobj_f01_i01_d20 0.833327121321942",
  "bbob-biobj_f01_i01_d40 0.833285390254100",
  "bbob-biobj_f01_i02_d02 0.833332720751734",
  "bbob-biobj_f01_i02_d03 0.833331088535920",
  "bbob-biobj_f01_i02_d05 0.833326851066133",
  "bbob-biobj_f01_i02_d10 0.833328914056953",
  "bbob-biobj_f01_i02_d20 0.833327494791052",
  "bbob-biobj_f01_i02_d40 0.833285385028953",
  "bbob-biobj_f01_i03_d02 0.833332349386930",
  "bbob-biobj_f01_i03_d03 0.833330739812853",
  "bbob-biobj_f01_i03_d05 0.833326555969321",
  "bbob-biobj_f01_i03_d10 0.833328819605771",
  "bbob-biobj_f01_i03_d20 0.833327457761651",
  "bbob-biobj_f01_i03_d40 0.833284614608700",
  "bbob-biobj_f01_i04_d02 0.833332499407236",
  "bbob-biobj_f01_i04_d03 0.833330724870132",
  "bbob-biobj_f01_i04_d05 0.833326640960598",
  "bbob-biobj_f01_i04_d10 0.833329375153094",
  "bbob-biobj_f01_i04_d20 0.833327270681467",
  "bbob-biobj_f01_i04_d40 0.833285857552923",
  "bbob-biobj_f01_i05_d02 0.833332512519638",
  "bbob-biobj_f01_i05_d03 0.833330830192043",
  "bbob-biobj_f01_i05_d05 0.833327050557970",
  "bbob-biobj_f01_i05_d10 0.833329037419504",
  "bbob-biobj_f01_i05_d20 0.833327476018437",
  "bbob-biobj_f01_i05_d40 0.833285593696607",
  "bbob-biobj_f01_i06_d02 0.833310937590996",
  "bbob-biobj_f01_i06_d03 0.833311215439312",
  "bbob-biobj_f01_i06_d05 0.833311611756160",
  "bbob-biobj_f01_i06_d10 0.833322275092647",
  "bbob-biobj_f01_i06_d20 0.833323776843884",
  "bbob-biobj_f01_i06_d40 0.832871381871906",
  "bbob-biobj_f01_i07_d02 0.833312631705348",
  "bbob-biobj_f01_i07_d03 0.833311900839859",
  "bbob-biobj_f01_i07_d05 0.833311287507580",
  "bbob-biobj_f01_i07_d10 0.833322245618162",
  "bbob-biobj_f01_i07_d20 0.833323789560072",
  "bbob-biobj_f01_i07_d40 0.832843499826846",
  "bbob-biobj_f01_i08_d02 0.833311659421989",
  "bbob-biobj_f01_i08_d03 0.833317899895486",
  "bbob-biobj_f01_i08_d05 0.833309999210594",
  "bbob-biobj_f01_i08_d10 0.833311245739727",
  "bbob-biobj_f01_i08_d20 0.833321907248139",
  "bbob-biobj_f01_i08_d40 0.832840913556568",
  "bbob-biobj_f01_i09_d02 0.833310991347934",
  "bbob-biobj_f01_i09_d03 0.833317530366279",
  "bbob-biobj_f01_i09_d05 0.833311307842818",
  "bbob-biobj_f01_i09_d10 0.833321913943093",
  "bbob-biobj_f01_i09_d20 0.833323790069312",
  "bbob-biobj_f01_i09_d40 0.832837360764175",
  "bbob-biobj_f01_i10_d02 0.833312368748217",
  "bbob-biobj_f01_i10_d03 0.833312231262573",
  "bbob-biobj_f01_i10_d05 0.833311675561901",
  "bbob-biobj_f01_i10_d10 0.833320870526659",
  "bbob-biobj_f01_i10_d20 0.833323597040930",
  "bbob-biobj_f01_i10_d40 0.832841129178451",
  "bbob-biobj_f02_i01_d02 0.995822556060843",
  "bbob-biobj_f02_i01_d03 0.879310058709001",
  "bbob-biobj_f02_i01_d05 0.953382600215256",
  "bbob-biobj_f02_i01_d10 0.978187425331875",
  "bbob-biobj_f02_i01_d20 0.951897811219661",
  "bbob-biobj_f02_i01_d40 0.949759312151031",
  "bbob-biobj_f02_i02_d02 0.917892766375743",
  "bbob-biobj_f02_i02_d03 0.981135803628767",
  "bbob-biobj_f02_i02_d05 0.966472109431022",
  "bbob-biobj_f02_i02_d10 0.954018561983206",
  "bbob-biobj_f02_i02_d20 0.980910929306251",
  "bbob-biobj_f02_i02_d40 0.967748482711833",
  "bbob-biobj_f02_i03_d02 0.990979154816691",
  "bbob-biobj_f02_i03_d03 0.952600931288877",
  "bbob-biobj_f02_i03_d05 0.950362077519103",
  "bbob-biobj_f02_i03_d10 0.890661671505104",
  "bbob-biobj_f02_i03_d20 0.971896136710556",
  "bbob-biobj_f02_i03_d40 0.962163282727804",
  "bbob-biobj_f02_i04_d02 0.956280220644638",
  "bbob-biobj_f02_i04_d03 0.889456798390216",
  "bbob-biobj_f02_i04_d05 0.881361353479823",
  "bbob-biobj_f02_i04_d10 0.972587804484677",
  "bbob-biobj_f02_i04_d20 0.977350713649177",
  "bbob-biobj_f02_i04_d40 0.976434229826955",
  "bbob-biobj_f02_i05_d02 0.960749286074231",
  "bbob-biobj_f02_i05_d03 0.924022338253155",
  "bbob-biobj_f02_i05_d05 0.861589171038206",
  "bbob-biobj_f02_i05_d10 0.933909238971449",
  "bbob-biobj_f02_i05_d20 0.962916791231407",
  "bbob-biobj_f02_i05_d40 0.960214010040557",
  "bbob-biobj_f02_i06_d02 0.875339986379539",
  "bbob-biobj_f02_i06_d03 0.948354623940449",
  "bbob-biobj_f02_i06_d05 0.910398345181518",
  "bbob-biobj_f02_i06_d10 0.977159747327357",
  "bbob-biobj_f02_i06_d20 0.971855499769274",
  "bbob-biobj_f02_i06_d40 0.972239569031113",
  "bbob-biobj_f02_i07_d02 0.832353335558318",
  "bbob-biobj_f02_i07_d03 0.875157618659083",
  "bbob-biobj_f02_i07_d05 0.912510812335122",
  "bbob-biobj_f02_i07_d10 0.958795661268572",
  "bbob-biobj_f02_i07_d20 0.977371875490846",
  "bbob-biobj_f02_i07_d40 0.967907676344765",
  "bbob-biobj_f02_i08_d02 0.829416511270900",
  "bbob-biobj_f02_i08_d03 0.974507626263228",
  "bbob-biobj_f02_i08_d05 0.935980638569705",
  "bbob-biobj_f02_i08_d10 0.922547833535203",
  "bbob-biobj_f02_i08_d20 0.967053592877665",
  "bbob-biobj_f02_i08_d40 0.962870102717631",
  "bbob-biobj_f02_i09_d02 0.991876232377447",
  "bbob-biobj_f02_i09_d03 0.978460183352107",
  "bbob-biobj_f02_i09_d05 0.942483669233472",
  "bbob-biobj_f02_i09_d10 0.953981656721534",
  "bbob-biobj_f02_i09_d20 0.961021472337196",
  "bbob-biobj_f02_i09_d40 0.973754082039648",
  "bbob-biobj_f02_i10_d02 0.947412366488691",
  "bbob-biobj_f02_i10_d03 0.991594440578253",
  "bbob-biobj_f02_i10_d05 0.948401779285958",
  "bbob-biobj_f02_i10_d10 0.965506887203965",
  "bbob-biobj_f02_i10_d20 0.970927452537725",
  "bbob-biobj_f02_i10_d40 0.963522134737799",
  "bbob-biobj_f03_i01_d02 0.811763857049145",
  "bbob-biobj_f03_i01_d03 0.974987112638223",
  "bbob-biobj_f03_i01_d05 0.846749947485651",
  "bbob-biobj_f03_i01_d10 0.916890141903379",
  "bbob-biobj_f03_i01_d20 0.887081921593395",
  "bbob-biobj_f03_i01_d40 0.876756825730484",
  "bbob-biobj_f03_i02_d02 0.870718437701562",
  "bbob-biobj_f03_i02_d03 0.845124188803593",
  "bbob-biobj_f03_i02_d05 0.961673728477631",
  "bbob-biobj_f03_i02_d10 0.980161858850881",
  "bbob-biobj_f03_i02_d20 0.950048245936200",
  "bbob-biobj_f03_i02_d40 0.940791241893415",
  "bbob-biobj_f03_i03_d02 0.843026945062331",
  "bbob-biobj_f03_i03_d03 0.860125236870242",
  "bbob-biobj_f03_i03_d05 0.836864473044929",
  "bbob-biobj_f03_i03_d10 0.985498403006374",
  "bbob-biobj_f03_i03_d20 0.867649235640338",
  "bbob-biobj_f03_i03_d40 0.885306102324499",
  "bbob-biobj_f03_i04_d02 0.816336921376108",
  "bbob-biobj_f03_i04_d03 0.965016698720033",
  "bbob-biobj_f03_i04_d05 0.832971705407539",
  "bbob-biobj_f03_i04_d10 0.908651290167361",
  "bbob-biobj_f03_i04_d20 0.932774984670678",
  "bbob-biobj_f03_i04_d40 0.910804664993381",
  "bbob-biobj_f03_i05_d02 0.854018820793507",
  "bbob-biobj_f03_i05_d03 0.879234581540502",
  "bbob-biobj_f03_i05_d05 0.959274894132135",
  "bbob-biobj_f03_i05_d10 0.881578736863151",
  "bbob-biobj_f03_i05_d20 0.875606096077933",
  "bbob-biobj_f03_i05_d40 0.897334402104937",
  "bbob-biobj_f03_i06_d02 0.830299197227874",
  "bbob-biobj_f03_i06_d03 0.944928349862207",
  "bbob-biobj_f03_i06_d05 0.970494904084209",
  "bbob-biobj_f03_i06_d10 0.844312649003636",
  "bbob-biobj_f03_i06_d20 0.910346449561183",
  "bbob-biobj_f03_i06_d40 0.912768614220291",
  "bbob-biobj_f03_i07_d02 0.868088162374295",
  "bbob-biobj_f03_i07_d03 0.869292897359135",
  "bbob-biobj_f03_i07_d05 0.903824108643888",
  "bbob-biobj_f03_i07_d10 0.845677498117250",
  "bbob-biobj_f03_i07_d20 0.923282074457029",
  "bbob-biobj_f03_i07_d40 0.907237469842650",
  "bbob-biobj_f03_i08_d02 0.990906422555066",
  "bbob-biobj_f03_i08_d03 0.835063452780516",
  "bbob-biobj_f03_i08_d05 0.922451009670776",
  "bbob-biobj_f03_i08_d10 0.886849156132951",
  "bbob-biobj_f03_i08_d20 0.924580691210482",
  "bbob-biobj_f03_i08_d40 0.901577541904608",
  "bbob-biobj_f03_i09_d02 0.813018786454735",
  "bbob-biobj_f03_i09_d03 0.929413275221192",
  "bbob-biobj_f03_i09_d05 0.852617246373087",
  "bbob-biobj_f03_i09_d10 0.988145373016528",
  "bbob-biobj_f03_i09_d20 0.891244651866651",
  "bbob-biobj_f03_i09_d40 0.960184519519258",
  "bbob-biobj_f03_i10_d02 0.806775641141177",
  "bbob-biobj_f03_i10_d03 0.889790706676965",
  "bbob-biobj_f03_i10_d05 0.872221017258980",
  "bbob-biobj_f03_i10_d10 0.838762983764005",
  "bbob-biobj_f03_i10_d20 0.932819414252023",
  "bbob-biobj_f03_i10_d40 0.919022006355222",
  "bbob-biobj_f04_i01_d02 0.965338549772770",
  "bbob-biobj_f04_i01_d03 0.968687234742656",
  "bbob-biobj_f04_i01_d05 0.943987114220437",
  "bbob-biobj_f04_i01_d10 0.944834478127607",
  "bbob-biobj_f04_i01_d20 0.935871868899059",
  "bbob-biobj_f04_i01_d40 0.938323088236502",
  "bbob-biobj_f04_i02_d02 0.970390625425527",
  "bbob-biobj_f04_i02_d03 0.955014013214034",
  "bbob-biobj_f04_i02_d05 0.963228111629313",
  "bbob-biobj_f04_i02_d10 0.954672682142632",
  "bbob-biobj_f04_i02_d20 0.941757018716809",
  "bbob-biobj_f04_i02_d40 0.942529855761297",
  "bbob-biobj_f04_i03_d02 0.971131579953310",
  "bbob-biobj_f04_i03_d03 0.910479816491977",
  "bbob-biobj_f04_i03_d05 0.937876509647381",
  "bbob-biobj_f04_i03_d10 0.951452121813764",
  "bbob-biobj_f04_i03_d20 0.931485725184838",
  "bbob-biobj_f04_i03_d40 0.935898816310602",
  "bbob-biobj_f04_i04_d02 0.977012414280338",
  "bbob-biobj_f04_i04_d03 0.994699593709363",
  "bbob-biobj_f04_i04_d05 0.944467965817465",
  "bbob-biobj_f04_i04_d10 0.936535080068134",
  "bbob-biobj_f04_i04_d20 0.942530799149335",
  "bbob-biobj_f04_i04_d40 0.936521563276061",
  "bbob-biobj_f04_i05_d02 0.924874195613136",
  "bbob-biobj_f04_i05_d03 0.923160754609184",
  "bbob-biobj_f04_i05_d05 0.942088788049264",
  "bbob-biobj_f04_i05_d10 0.941629230113325",
  "bbob-biobj_f04_i05_d20 0.948792969456780",
  "bbob-biobj_f04_i05_d40 0.941318449117935",
  "bbob-biobj_f04_i06_d02 0.972691840870603",
  "bbob-biobj_f04_i06_d03 0.943232435272800",
  "bbob-biobj_f04_i06_d05 0.950574098979926",
  "bbob-biobj_f04_i06_d10 0.952949532886138",
  "bbob-biobj_f04_i06_d20 0.950717154295874",
  "bbob-biobj_f04_i06_d40 0.941893144828156",
  "bbob-biobj_f04_i07_d02 0.955654394311334",
  "bbob-biobj_f04_i07_d03 0.967689784486838",
  "bbob-biobj_f04_i07_d05 0.954731380416701",
  "bbob-biobj_f04_i07_d10 0.963145382230048",
  "bbob-biobj_f04_i07_d20 0.945560736705184",
  "bbob-biobj_f04_i07_d40 0.946037040689202",
  "bbob-biobj_f04_i08_d02 0.907722989346049",
  "bbob-biobj_f04_i08_d03 0.921824431251827",
  "bbob-biobj_f04_i08_d05 0.959550859031763",
  "bbob-biobj_f04_i08_d10 0.948087750040469",
  "bbob-biobj_f04_i08_d20 0.941557710236498",
  "bbob-biobj_f04_i08_d40 0.933853737661553",
  "bbob-biobj_f04_i09_d02 0.810183898248207",
  "bbob-biobj_f04_i09_d03 0.940436461325365",
  "bbob-biobj_f04_i09_d05 0.933551325277965",
  "bbob-biobj_f04_i09_d10 0.942886500254279",
  "bbob-biobj_f04_i09_d20 0.932995494849010",
  "bbob-biobj_f04_i09_d40 0.935741488687228",
  "bbob-biobj_f04_i10_d02 0.954293699220526",
  "bbob-biobj_f04_i10_d03 0.918850718457585",
  "bbob-biobj_f04_i10_d05 0.935579264509031",
  "bbob-biobj_f04_i10_d10 0.926094048836512",
  "bbob-biobj_f04_i10_d20 0.933732740883168",
  "bbob-biobj_f04_i10_d40 0.935305333133984",
  "bbob-biobj_f05_i01_d02 0.754127865991100",
  "bbob-biobj_f05_i01_d03 0.728471324736124",
  "bbob-biobj_f05_i01_d05 0.732282352556790",
  "bbob-biobj_f05_i01_d10 0.714066803669735",
  "bbob-biobj_f05_i01_d20 0.694475751705236",
  "bbob-biobj_f05_i01_d40 0.709846837542285",
  "bbob-biobj_f05_i02_d02 0.954989359408118",
  "bbob-biobj_f05_i02_d03 0.688040297931016",
  "bbob-biobj_f05_i02_d05 0.714728659711816",
  "bbob-biobj_f05_i02_d10 0.730732920139007",
  "bbob-biobj_f05_i02_d20 0.689198535536362",
  "bbob-biobj_f05_i02_d40 0.698095234085181",
  "bbob-biobj_f05_i03_d02 0.684482309631958",
  "bbob-biobj_f05_i03_d03 0.802884315600055",
  "bbob-biobj_f05_i03_d05 0.699887927831400",
  "bbob-biobj_f05_i03_d10 0.683926423161017",
  "bbob-biobj_f05_i03_d20 0.697132994973134",
  "bbob-biobj_f05_i03_d40 0.694169722678116",
  "bbob-biobj_f05_i04_d02 0.878629177779160",
  "bbob-biobj_f05_i04_d03 0.744993525663924",
  "bbob-biobj_f05_i04_d05 0.776082723541421",
  "bbob-biobj_f05_i04_d10 0.716302292985648",
  "bbob-biobj_f05_i04_d20 0.705263541171486",
  "bbob-biobj_f05_i04_d40 0.699934703162273",
  "bbob-biobj_f05_i05_d02 0.926274375941586",
  "bbob-biobj_f05_i05_d03 0.701514838673870",
  "bbob-biobj_f05_i05_d05 0.737160017743642",
  "bbob-biobj_f05_i05_d10 0.749489719851661",
  "bbob-biobj_f05_i05_d20 0.695487080627162",
  "bbob-biobj_f05_i05_d40 0.698154791518541",
  "bbob-biobj_f05_i06_d02 0.885902522164927",
  "bbob-biobj_f05_i06_d03 0.752872122704611",
  "bbob-biobj_f05_i06_d05 0.777039555624318",
  "bbob-biobj_f05_i06_d10 0.760017111363263",
  "bbob-biobj_f05_i06_d20 0.701404536132544",
  "bbob-biobj_f05_i06_d40 0.703103765834964",
  "bbob-biobj_f05_i07_d02 0.733263314113369",
  "bbob-biobj_f05_i07_d03 0.837078460675662",
  "bbob-biobj_f05_i07_d05 0.732314245319365",
  "bbob-biobj_f05_i07_d10 0.704066489270897",
  "bbob-biobj_f05_i07_d20 0.714267651335895",
  "bbob-biobj_f05_i07_d40 0.701650339436692",
  "bbob-biobj_f05_i08_d02 0.720804608854457",
  "bbob-biobj_f05_i08_d03 0.718897770204007",
  "bbob-biobj_f05_i08_d05 0.719008134327259",
  "bbob-biobj_f05_i08_d10 0.695227127767466",
  "bbob-biobj_f05_i08_d20 0.721370358468108",
  "bbob-biobj_f05_i08_d40 0.698384816370664",
  "bbob-biobj_f05_i09_d02 0.843335036506577",
  "bbob-biobj_f05_i09_d03 0.678630261985751",
  "bbob-biobj_f05_i09_d05 0.774319592228313",
  "bbob-biobj_f05_i09_d10 0.702659596290559",
  "bbob-biobj_f05_i09_d20 0.700085059778838",
  "bbob-biobj_f05_i09_d40 0.695297011513272",
  "bbob-biobj_f05_i10_d02 0.781255377607552",
  "bbob-biobj_f05_i10_d03 0.935404825799813",
  "bbob-biobj_f05_i10_d05 0.765104988680107",
  "bbob-biobj_f05_i10_d10 0.703626540120432",
  "bbob-biobj_f05_i10_d20 0.695325429591975",
  "bbob-biobj_f05_i10_d40 0.698878500876222",
  "bbob-biobj_f06_i01_d02 0.667254107678143",
  "bbob-biobj_f06_i01_d03 0.954291256239000",
  "bbob-biobj_f06_i01_d05 0.846010493810613",
  "bbob-biobj_f06_i01_d10 0.937006022241783",
  "bbob-biobj_f06_i01_d20 0.931057894098642",
  "bbob-biobj_f06_i01_d40 0.923065329882169",
  "bbob-biobj_f06_i02_d02 0.901470247301924",
  "bbob-biobj_f06_i02_d03 0.863891197188673",
  "bbob-biobj_f06_i02_d05 0.867718960905413",
  "bbob-biobj_f06_i02_d10 0.875236429261449",
  "bbob-biobj_f06_i02_d20 0.910529843280032",
  "bbob-biobj_f06_i02_d40 0.942992942865892",
  "bbob-biobj_f06_i03_d02 0.884299612080317",
  "bbob-biobj_f06_i03_d03 0.833636623657147",
  "bbob-biobj_f06_i03_d05 0.848588634199072",
  "bbob-biobj_f06_i03_d10 0.895936158446357",
  "bbob-biobj_f06_i03_d20 0.932801430955019",
  "bbob-biobj_f06_i03_d40 0.902577785053624",
  "bbob-biobj_f06_i04_d02 0.945804104911204",
  "bbob-biobj_f06_i04_d03 0.921370202901199",
  "bbob-biobj_f06_i04_d05 0.935941754702140",
  "bbob-biobj_f06_i04_d10 0.930057004748218",
  "bbob-biobj_f06_i04_d20 0.901076360810312",
  "bbob-biobj_f06_i04_d40 0.876237901084423",
  "bbob-biobj_f06_i05_d02 0.942603152334283",
  "bbob-biobj_f06_i05_d03 0.899135548418015",
  "bbob-biobj_f06_i05_d05 0.930564101030568",
  "bbob-biobj_f06_i05_d10 0.743358909511397",
  "bbob-biobj_f06_i05_d20 0.918102018528550",
  "bbob-biobj_f06_i05_d40 0.858518894302918",
  "bbob-biobj_f06_i06_d02 0.899055246706590",
  "bbob-biobj_f06_i06_d03 0.836219832639499",
  "bbob-biobj_f06_i06_d05 0.811598424960342",
  "bbob-biobj_f06_i06_d10 0.928885711526853",
  "bbob-biobj_f06_i06_d20 0.868865550589963",
  "bbob-biobj_f06_i06_d40 0.913780766426345",
  "bbob-biobj_f06_i07_d02 0.813313852808380",
  "bbob-biobj_f06_i07_d03 0.899712678629379",
  "bbob-biobj_f06_i07_d05 0.876974920873286",
  "bbob-biobj_f06_i07_d10 0.815085734775849",
  "bbob-biobj_f06_i07_d20 0.935537984112798",
  "bbob-biobj_f06_i07_d40 0.917999648217449",
  "bbob-biobj_f06_i08_d02 0.910550218635067",
  "bbob-biobj_f06_i08_d03 0.667359764307318",
  "bbob-biobj_f06_i08_d05 0.937284426078644",
  "bbob-biobj_f06_i08_d10 0.930145807981647",
  "bbob-biobj_f06_i08_d20 0.910941583830995",
  "bbob-biobj_f06_i08_d40 0.925749531689669",
  "bbob-biobj_f06_i09_d02 0.675896080189426",
  "bbob-biobj_f06_i09_d03 0.867811302877810",
  "bbob-biobj_f06_i09_d05 0.897281153808301",
  "bbob-biobj_f06_i09_d10 0.845108748422141",
  "bbob-biobj_f06_i09_d20 0.949712352781937",
  "bbob-biobj_f06_i09_d40 0.945237374202699",
  "bbob-biobj_f06_i10_d02 0.882420980315420",
  "bbob-biobj_f06_i10_d03 0.907395391261030",
  "bbob-biobj_f06_i10_d05 0.905171717500438",
  "bbob-biobj_f06_i10_d10 0.906485259223307",
  "bbob-biobj_f06_i10_d20 0.901299605350245",
  "bbob-biobj_f06_i10_d40 0.921746795175045",
  "bbob-biobj_f07_i01_d02 0.936972399930132",
  "bbob-biobj_f07_i01_d03 0.937571173307641",
  "bbob-biobj_f07_i01_d05 0.860201953501511",
  "bbob-biobj_f07_i01_d10 0.896577394043478",
  "bbob-biobj_f07_i01_d20 0.942485386120226",
  "bbob-biobj_f07_i01_d40 0.909059154290133",
  "bbob-biobj_f07_i02_d02 0.906340757914227",
  "bbob-biobj_f07_i02_d03 0.923759786548036",
  "bbob-biobj_f07_i02_d05 0.893273257747443",
  "bbob-biobj_f07_i02_d10 0.895698212813317",
  "bbob-biobj_f07_i02_d20 0.898580560160935",
  "bbob-biobj_f07_i02_d40 0.861964006417624",
  "bbob-biobj_f07_i03_d02 0.886134179761151",
  "bbob-biobj_f07_i03_d03 0.921397888403983",
  "bbob-biobj_f07_i03_d05 0.868176173139134",
  "bbob-biobj_f07_i03_d10 0.893625883670611",
  "bbob-biobj_f07_i03_d20 0.896944660506832",
  "bbob-biobj_f07_i03_d40 0.887841572402023",
  "bbob-biobj_f07_i04_d02 0.870759578365989",
  "bbob-biobj_f07_i04_d03 0.933979945799374",
  "bbob-biobj_f07_i04_d05 0.870194778631520",
  "bbob-biobj_f07_i04_d10 0.884013788375882",
  "bbob-biobj_f07_i04_d20 0.893221823670617",
  "bbob-biobj_f07_i04_d40 0.903232224011270",
  "bbob-biobj_f07_i05_d02 0.911523142200378",
  "bbob-biobj_f07_i05_d03 0.887626374448897",
  "bbob-biobj_f07_i05_d05 0.911637907560092",
  "bbob-biobj_f07_i05_d10 0.868092909928177",
  "bbob-biobj_f07_i05_d20 0.886708365225772",
  "bbob-biobj_f07_i05_d40 0.912451735893751",
  "bbob-biobj_f07_i06_d02 0.937825227004125",
  "bbob-biobj_f07_i06_d03 0.945843590726311",
  "bbob-biobj_f07_i06_d05 0.906530098104496",
  "bbob-biobj_f07_i06_d10 0.889579864412070",
  "bbob-biobj_f07_i06_d20 0.881646714085581",
  "bbob-biobj_f07_i06_d40 0.872446728348711",
  "bbob-biobj_f07_i07_d02 0.871120237220164",
  "bbob-biobj_f07_i07_d03 0.910072602845794",
  "bbob-biobj_f07_i07_d05 0.881805758304167",
  "bbob-biobj_f07_i07_d10 0.894012109079692",
  "bbob-biobj_f07_i07_d20 0.893700664927148",
  "bbob-biobj_f07_i07_d40 0.824135589500768",
  "bbob-biobj_f07_i08_d02 0.849304845123815",
  "bbob-biobj_f07_i08_d03 0.909048667993127",
  "bbob-biobj_f07_i08_d05 0.838503958734642",
  "bbob-biobj_f07_i08_d10 0.916245433495188",
  "bbob-biobj_f07_i08_d20 0.886236839114243",
  "bbob-biobj_f07_i08_d40 0.893434654917281",
  "bbob-biobj_f07_i09_d02 0.877375163908561",
  "bbob-biobj_f07_i09_d03 0.928121321538075",
  "bbob-biobj_f07_i09_d05 0.886433217520399",
  "bbob-biobj_f07_i09_d10 0.909555880081429",
  "bbob-biobj_f07_i09_d20 0.897381981901614",
  "bbob-biobj_f07_i09_d40 0.839855674570544",
  "bbob-biobj_f07_i10_d02 0.907365720229927",
  "bbob-biobj_f07_i10_d03 0.918960996777705",
  "bbob-biobj_f07_i10_d05 0.891774404805287",
  "bbob-biobj_f07_i10_d10 0.897479586502742",
  "bbob-biobj_f07_i10_d20 0.932496622583228",
  "bbob-biobj_f07_i10_d40 0.859515671273887",
  "bbob-biobj_f08_i01_d02 0.903849269051818",
  "bbob-biobj_f08_i01_d03 0.911798619642005",
  "bbob-biobj_f08_i01_d05 0.942655000145460",
  "bbob-biobj_f08_i01_d10 0.982213909300857",
  "bbob-biobj_f08_i01_d20 0.967978562063254",
  "bbob-biobj_f08_i01_d40 0.941343715004475",
  "bbob-biobj_f08_i02_d02 0.784765078034902",
  "bbob-biobj_f08_i02_d03 0.882018389110311",
  "bbob-biobj_f08_i02_d05 0.909387465399333",
  "bbob-biobj_f08_i02_d10 0.915933949542192",
  "bbob-biobj_f08_i02_d20 0.902120755033762",
  "bbob-biobj_f08_i02_d40 0.942550840384823",
  "bbob-biobj_f08_i03_d02 0.748604742133252",
  "bbob-biobj_f08_i03_d03 0.850945639689013",
  "bbob-biobj_f08_i03_d05 0.805258279611532",
  "bbob-biobj_f08_i03_d10 0.930897297074430",
  "bbob-biobj_f08_i03_d20 0.948815915384450",
  "bbob-biobj_f08_i03_d40 0.882347452604998",
  "bbob-biobj_f08_i04_d02 0.743267154821628",
  "bbob-biobj_f08_i04_d03 0.667005341586434",
  "bbob-biobj_f08_i04_d05 0.927686052419946",
  "bbob-biobj_f08_i04_d10 0.950780081828940",
  "bbob-biobj_f08_i04_d20 0.955613618200345",
  "bbob-biobj_f08_i04_d40 0.923529281274225",
  "bbob-biobj_f08_i05_d02 0.865136515479480",
  "bbob-biobj_f08_i05_d03 0.893995583789640",
  "bbob-biobj_f08_i05_d05 0.917809142365534",
  "bbob-biobj_f08_i05_d10 0.930294235467496",
  "bbob-biobj_f08_i05_d20 0.906624299145002",
  "bbob-biobj_f08_i05_d40 0.944010808321168",
  "bbob-biobj_f08_i06_d02 0.829435542844206",
  "bbob-biobj_f08_i06_d03 0.895110159139378",
  "bbob-biobj_f08_i06_d05 0.879081523914896",
  "bbob-biobj_f08_i06_d10 0.906270697142924",
  "bbob-biobj_f08_i06_d20 0.923422132859498",
  "bbob-biobj_f08_i06_d40 0.898758559470351",
  "bbob-biobj_f08_i07_d02 0.933678913852236",
  "bbob-biobj_f08_i07_d03 0.909557034401524",
  "bbob-biobj_f08_i07_d05 0.895697751653331",
  "bbob-biobj_f08_i07_d10 0.869943414311157",
  "bbob-biobj_f08_i07_d20 0.947297927998108",
  "bbob-biobj_f08_i07_d40 0.867734452605647",
  "bbob-biobj_f08_i08_d02 0.901598955987436",
  "bbob-biobj_f08_i08_d03 0.931901261513024",
  "bbob-biobj_f08_i08_d05 0.903388612365953",
  "bbob-biobj_f08_i08_d10 0.921322705536751",
  "bbob-biobj_f08_i08_d20 0.930910689290388",
  "bbob-biobj_f08_i08_d40 0.818778993301297",
  "bbob-biobj_f08_i09_d02 0.904686562833679",
  "bbob-biobj_f08_i09_d03 0.938036286897513",
  "bbob-biobj_f08_i09_d05 0.826898944668083",
  "bbob-biobj_f08_i09_d10 0.916335831477329",
  "bbob-biobj_f08_i09_d20 0.903325713468571",
  "bbob-biobj_f08_i09_d40 0.926009566367219",
  "bbob-biobj_f08_i10_d02 0.904852352345612",
  "bbob-biobj_f08_i10_d03 0.931047258452528",
  "bbob-biobj_f08_i10_d05 0.946338993945063",
  "bbob-biobj_f08_i10_d10 0.915755506148993",
  "bbob-biobj_f08_i10_d20 0.912376835181859",
  "bbob-biobj_f08_i10_d40 0.959408270528453",
  "bbob-biobj_f09_i01_d02 0.925657494081605",
  "bbob-biobj_f09_i01_d03 0.904195153094803",
  "bbob-biobj_f09_i01_d05 0.932178589563142",
  "bbob-biobj_f09_i01_d10 0.940798324900223",
  "bbob-biobj_f09_i01_d20 0.960314050378320",
  "bbob-biobj_f09_i01_d40 0.966647252965589",
  "bbob-biobj_f09_i02_d02 0.977793679012675",
  "bbob-biobj_f09_i02_d03 0.992207203222584",
  "bbob-biobj_f09_i02_d05 0.961852919392127",
  "bbob-biobj_f09_i02_d10 0.975989706087452",
  "bbob-biobj_f09_i02_d20 0.962606945889951",
  "bbob-biobj_f09_i02_d40 0.963751992610271",
  "bbob-biobj_f09_i03_d02 0.968705721628269",
  "bbob-biobj_f09_i03_d03 0.986085301970526",
  "bbob-biobj_f09_i03_d05 0.930450094656895",
  "bbob-biobj_f09_i03_d10 0.955041683216728",
  "bbob-biobj_f09_i03_d20 0.970330750194347",
  "bbob-biobj_f09_i03_d40 0.969381844759192",
  "bbob-biobj_f09_i04_d02 0.948342390905787",
  "bbob-biobj_f09_i04_d03 0.940843190553506",
  "bbob-biobj_f09_i04_d05 0.950694330245476",
  "bbob-biobj_f09_i04_d10 0.944319962876144",
  "bbob-biobj_f09_i04_d20 0.961726548135361",
  "bbob-biobj_f09_i04_d40 0.963091309867986",
  "bbob-biobj_f09_i05_d02 0.860780287881932",
  "bbob-biobj_f09_i05_d03 0.939788217960798",
  "bbob-biobj_f09_i05_d05 0.968560017482823",
  "bbob-biobj_f09_i05_d10 0.940730585183576",
  "bbob-biobj_f09_i05_d20 0.963733002210023",
  "bbob-biobj_f09_i05_d40 0.967130700974806",
  "bbob-biobj_f09_i06_d02 0.957368168761571",
  "bbob-biobj_f09_i06_d03 0.987784971852043",
  "bbob-biobj_f09_i06_d05 0.910205061988753",
  "bbob-biobj_f09_i06_d10 0.958607574186893",
  "bbob-biobj_f09_i06_d20 0.969286622421876",
  "bbob-biobj_f09_i06_d40 0.960757674665563",
  "bbob-biobj_f09_i07_d02 0.990545336056762",
  "bbob-biobj_f09_i07_d03 0.962688656913570",
  "bbob-biobj_f09_i07_d05 0.974588312387505",
  "bbob-biobj_f09_i07_d10 0.974596764128842",
  "bbob-biobj_f09_i07_d20 0.975332759310894",
  "bbob-biobj_f09_i07_d40 0.972038954692242",
  "bbob-biobj_f09_i08_d02 0.889467229726116",
  "bbob-biobj_f09_i08_d03 0.964097402483471",
  "bbob-biobj_f09_i08_d05 0.959671642037049",
  "bbob-biobj_f09_i08_d10 0.954108400973619",
  "bbob-biobj_f09_i08_d20 0.950910459097299",
  "bbob-biobj_f09_i08_d40 0.959804989720909",
  "bbob-biobj_f09_i09_d02 0.950445732912183",
  "bbob-biobj_f09_i09_d03 0.983613374221624",
  "bbob-biobj_f09_i09_d05 0.964952044615142",
  "bbob-biobj_f09_i09_d10 0.979819305252771",
  "bbob-biobj_f09_i09_d20 0.978341291145915",
  "bbob-biobj_f09_i09_d40 0.967545119559943",
  "bbob-biobj_f09_i10_d02 0.869847511977600",
  "bbob-biobj_f09_i10_d03 0.943417198143025",
  "bbob-biobj_f09_i10_d05 0.920128104933707",
  "bbob-biobj_f09_i10_d10 0.944139414485540",
  "bbob-biobj_f09_i10_d20 0.957501235248264",
  "bbob-biobj_f09_i10_d40 0.964012416299468",
  "bbob-biobj_f10_i01_d02 0.922987549727637",
  "bbob-biobj_f10_i01_d03 0.927567696810629",
  "bbob-biobj_f10_i01_d05 0.867906861359309",
  "bbob-biobj_f10_i01_d10 0.879649116695308",
  "bbob-biobj_f10_i01_d20 0.840938246208460",
  "bbob-biobj_f10_i01_d40 0.729356054801270",
  "bbob-biobj_f10_i02_d02 0.889244311731490",
  "bbob-biobj_f10_i02_d03 0.883786041156248",
  "bbob-biobj_f10_i02_d05 0.898112903164327",
  "bbob-biobj_f10_i02_d10 0.798467735698137",
  "bbob-biobj_f10_i02_d20 0.812103327654335",
  "bbob-biobj_f10_i02_d40 0.663812162941286",
  "bbob-biobj_f10_i03_d02 0.921417927383868",
  "bbob-biobj_f10_i03_d03 0.901133084432264",
  "bbob-biobj_f10_i03_d05 0.890596927609057",
  "bbob-biobj_f10_i03_d10 0.710806047115800",
  "bbob-biobj_f10_i03_d20 0.798300946578455",
  "bbob-biobj_f10_i03_d40 0.747898316341296",
  "bbob-biobj_f10_i04_d02 0.942089279409840",
  "bbob-biobj_f10_i04_d03 0.932701730874906",
  "bbob-biobj_f10_i04_d05 0.946969102718957",
  "bbob-biobj_f10_i04_d10 0.906381933781305",
  "bbob-biobj_f10_i04_d20 0.800179081775505",
  "bbob-biobj_f10_i04_d40 0.638788079762735",
  "bbob-biobj_f10_i05_d02 0.940003470249227",
  "bbob-biobj_f10_i05_d03 0.934362602797827",
  "bbob-biobj_f10_i05_d05 0.934945006905130",
  "bbob-biobj_f10_i05_d10 0.842819089279240",
  "bbob-biobj_f10_i05_d20 0.778529368085970",
  "bbob-biobj_f10_i05_d40 0.638345357359405",
  "bbob-biobj_f10_i06_d02 0.884469500115369",
  "bbob-biobj_f10_i06_d03 0.929295668097127",
  "bbob-biobj_f10_i06_d05 0.937687370315965",
  "bbob-biobj_f10_i06_d10 0.916289993621501",
  "bbob-biobj_f10_i06_d20 0.803794793811259",
  "bbob-biobj_f10_i06_d40 0.690318695065605",
  "bbob-biobj_f10_i07_d02 0.972906803928232",
  "bbob-biobj_f10_i07_d03 0.971443583283319",
  "bbob-biobj_f10_i07_d05 0.941658545474087",
  "bbob-biobj_f10_i07_d10 0.817874642361040",
  "bbob-biobj_f10_i07_d20 0.763203076667597",
  "bbob-biobj_f10_i07_d40 0.726574213658138",
  "bbob-biobj_f10_i08_d02 0.926499922535728",
  "bbob-biobj_f10_i08_d03 0.950920243869703",
  "bbob-biobj_f10_i08_d05 0.975224077117786",
  "bbob-biobj_f10_i08_d10 0.913936722938377",
  "bbob-biobj_f10_i08_d20 0.857597806783516",
  "bbob-biobj_f10_i08_d40 0.775190013320872",
  "bbob-biobj_f10_i09_d02 0.663675998347694",
  "bbob-biobj_f10_i09_d03 0.878892141679031",
  "bbob-biobj_f10_i09_d05 0.941987305956976",
  "bbob-biobj_f10_i09_d10 0.904811868533760",
  "bbob-biobj_f10_i09_d20 0.792532110086720",
  "bbob-biobj_f10_i09_d40 0.713839757359807",
  "bbob-biobj_f10_i10_d02 0.909582510068885",
  "bbob-biobj_f10_i10_d03 0.950967867587539",
  "bbob-biobj_f10_i10_d05 0.936923102991904",
  "bbob-biobj_f10_i10_d10 0.885441698175989",
  "bbob-biobj_f10_i10_d20 0.748638975223095",
  "bbob-biobj_f10_i10_d40 0.532515062577884",
  "bbob-biobj_f11_i01_d02 0.823972802033229",
  "bbob-biobj_f11_i01_d03 0.878620383750947",
  "bbob-biobj_f11_i01_d05 0.812583126216868",
  "bbob-biobj_f11_i01_d10 0.836574058545642",
  "bbob-biobj_f11_i01_d20 0.836376248355817",
  "bbob-biobj_f11_i01_d40 0.824905961028570",
  "bbob-biobj_f11_i02_d02 0.834474630879233",
  "bbob-biobj_f11_i02_d03 0.833333076244992",
  "bbob-biobj_f11_i02_d05 0.813661496792840",
  "bbob-biobj_f11_i02_d10 0.829487777342274",
  "bbob-biobj_f11_i02_d20 0.835464005007358",
  "bbob-biobj_f11_i02_d40 0.836265387756231",
  "bbob-biobj_f11_i03_d02 0.817436856949172",
  "bbob-biobj_f11_i03_d03 0.827369830625197",
  "bbob-biobj_f11_i03_d05 0.841063442990087",
  "bbob-biobj_f11_i03_d10 0.821924762364347",
  "bbob-biobj_f11_i03_d20 0.835205098414205",
  "bbob-biobj_f11_i03_d40 0.841616125488769",
  "bbob-biobj_f11_i04_d02 0.883087518329401",
  "bbob-biobj_f11_i04_d03 0.841523521058958",
  "bbob-biobj_f11_i04_d05 0.886487909384851",
  "bbob-biobj_f11_i04_d10 0.834674298423295",
  "bbob-biobj_f11_i04_d20 0.838290601518139",
  "bbob-biobj_f11_i04_d40 0.836113969942348",
  "bbob-biobj_f11_i05_d02 0.849343563564765",
  "bbob-biobj_f11_i05_d03 0.775580568622602",
  "bbob-biobj_f11_i05_d05 0.834554877740455",
  "bbob-biobj_f11_i05_d10 0.840894436332211",
  "bbob-biobj_f11_i05_d20 0.841626618613620",
  "bbob-biobj_f11_i05_d40 0.833031491068498",
  "bbob-biobj_f11_i06_d02 0.826920112167869",
  "bbob-biobj_f11_i06_d03 0.829287041102090",
  "bbob-biobj_f11_i06_d05 0.835720395191891",
  "bbob-biobj_f11_i06_d10 0.834251559345367",
  "bbob-biobj_f11_i06_d20 0.835707658355499",
  "bbob-biobj_f11_i06_d40 0.833786750883675",
  "bbob-biobj_f11_i07_d02 0.827116149813326",
  "bbob-biobj_f11_i07_d03 0.834074532916097",
  "bbob-biobj_f11_i07_d05 0.822460831442904",
  "bbob-biobj_f11_i07_d10 0.828482154972000",
  "bbob-biobj_f11_i07_d20 0.847672481820259",
  "bbob-biobj_f11_i07_d40 0.839251236466114",
  "bbob-biobj_f11_i08_d02 0.816775818075309",
  "bbob-biobj_f11_i08_d03 0.828518557019658",
  "bbob-biobj_f11_i08_d05 0.821643047395599",
  "bbob-biobj_f11_i08_d10 0.837582904190785",
  "bbob-biobj_f11_i08_d20 0.842724563815014",
  "bbob-biobj_f11_i08_d40 0.840480833441533",
  "bbob-biobj_f11_i09_d02 0.832754119980936",
  "bbob-biobj_f11_i09_d03 0.824016448830023",
  "bbob-biobj_f11_i09_d05 0.792480182887908",
  "bbob-biobj_f11_i09_d10 0.824396112083681",
  "bbob-biobj_f11_i09_d20 0.829209747434282",
  "bbob-biobj_f11_i09_d40 0.826889632605898",
  "bbob-biobj_f11_i10_d02 0.826860960440728",
  "bbob-biobj_f11_i10_d03 0.815591906718715",
  "bbob-biobj_f11_i10_d05 0.813067282203090",
  "bbob-biobj_f11_i10_d10 0.845229160749538",
  "bbob-biobj_f11_i10_d20 0.827367969570966",
  "bbob-biobj_f11_i10_d40 0.834461054280065",
  "bbob-biobj_f12_i01_d02 0.981018923952081",
  "bbob-biobj_f12_i01_d03 0.997338224021507",
  "bbob-biobj_f12_i01_d05 0.999992978794068",
  "bbob-biobj_f12_i01_d10 0.997429994637158",
  "bbob-biobj_f12_i01_d20 0.999547772931986",
  "bbob-biobj_f12_i01_d40 0.998374235929858",
  "bbob-biobj_f12_i02_d02 0.999911393192580",
  "bbob-biobj_f12_i02_d03 0.999935181229594",
  "bbob-biobj_f12_i02_d05 0.999910205093077",
  "bbob-biobj_f12_i02_d10 0.999751285425627",
  "bbob-biobj_f12_i02_d20 0.999572608820679",
  "bbob-biobj_f12_i02_d40 0.999213065148323",
  "bbob-biobj_f12_i03_d02 0.999915738926235",
  "bbob-biobj_f12_i03_d03 0.999886884015880",
  "bbob-biobj_f12_i03_d05 0.947357122629376",
  "bbob-biobj_f12_i03_d10 0.999904744702820",
  "bbob-biobj_f12_i03_d20 0.999772131274245",
  "bbob-biobj_f12_i03_d40 0.998237726952336",
  "bbob-biobj_f12_i04_d02 0.972296052339516",
  "bbob-biobj_f12_i04_d03 0.934335952177352",
  "bbob-biobj_f12_i04_d05 0.973685836073498",
  "bbob-biobj_f12_i04_d10 0.999925890363863",
  "bbob-biobj_f12_i04_d20 0.999900094345513",
  "bbob-biobj_f12_i04_d40 0.998362712511718",
  "bbob-biobj_f12_i05_d02 0.908759829060693",
  "bbob-biobj_f12_i05_d03 0.934982160085602",
  "bbob-biobj_f12_i05_d05 0.999480706194336",
  "bbob-biobj_f12_i05_d10 0.999713655577313",
  "bbob-biobj_f12_i05_d20 0.999223761178916",
  "bbob-biobj_f12_i05_d40 0.993030763487402",
  "bbob-biobj_f12_i06_d02 0.950885811476826",
  "bbob-biobj_f12_i06_d03 0.994616678835083",
  "bbob-biobj_f12_i06_d05 0.999940426090094",
  "bbob-biobj_f12_i06_d10 0.998777097328148",
  "bbob-biobj_f12_i06_d20 0.999850413404899",
  "bbob-biobj_f12_i06_d40 0.998533844858078",
  "bbob-biobj_f12_i07_d02 0.999976180189943",
  "bbob-biobj_f12_i07_d03 0.999953419858868",
  "bbob-biobj_f12_i07_d05 0.969366146413401",
  "bbob-biobj_f12_i07_d10 0.999954529278590",
  "bbob-biobj_f12_i07_d20 0.999764951492469",
  "bbob-biobj_f12_i07_d40 0.999563478959354",
  "bbob-biobj_f12_i08_d02 0.995571869066140",
  "bbob-biobj_f12_i08_d03 0.848124048856883",
  "bbob-biobj_f12_i08_d05 0.999927201099838",
  "bbob-biobj_f12_i08_d10 0.999740805374047",
  "bbob-biobj_f12_i08_d20 0.999406896621783",
  "bbob-biobj_f12_i08_d40 0.993531902131157",
  "bbob-biobj_f12_i09_d02 0.927605782925467",
  "bbob-biobj_f12_i09_d03 0.999273159571939",
  "bbob-biobj_f12_i09_d05 0.999973053677131",
  "bbob-biobj_f12_i09_d10 0.999911920855163",
  "bbob-biobj_f12_i09_d20 0.999705990242430",
  "bbob-biobj_f12_i09_d40 0.997636634472918",
  "bbob-biobj_f12_i10_d02 0.865737626755513",
  "bbob-biobj_f12_i10_d03 0.999978980513732",
  "bbob-biobj_f12_i10_d05 0.999675831298970",
  "bbob-biobj_f12_i10_d10 0.998438307553193",
  "bbob-biobj_f12_i10_d20 0.998088338358802",
  "bbob-biobj_f12_i10_d40 0.999237524283657",
  "bbob-biobj_f13_i01_d02 0.944887998977745",
  "bbob-biobj_f13_i01_d03 0.999495772338908",
  "bbob-biobj_f13_i01_d05 0.999802988484274",
  "bbob-biobj_f13_i01_d10 0.999291169614413",
  "bbob-biobj_f13_i01_d20 0.949989015858468",
  "bbob-biobj_f13_i01_d40 0.988688562193543",
  "bbob-biobj_f13_i02_d02 0.999759633726896",
  "bbob-biobj_f13_i02_d03 0.998700053055222",
  "bbob-biobj_f13_i02_d05 0.998548418781503",
  "bbob-biobj_f13_i02_d10 0.999205061000801",
  "bbob-biobj_f13_i02_d20 0.970026298418277",
  "bbob-biobj_f13_i02_d40 0.989827293270395",
  "bbob-biobj_f13_i03_d02 0.999951135499806",
  "bbob-biobj_f13_i03_d03 0.962145762074075",
  "bbob-biobj_f13_i03_d05 0.999818983161909",
  "bbob-biobj_f13_i03_d10 0.999218956587019",
  "bbob-biobj_f13_i03_d20 0.945577704630341",
  "bbob-biobj_f13_i03_d40 0.990935490401149",
  "bbob-biobj_f13_i04_d02 0.999963191344402",
  "bbob-biobj_f13_i04_d03 0.999978299131414",
  "bbob-biobj_f13_i04_d05 0.999810914634710",
  "bbob-biobj_f13_i04_d10 0.991556775275303",
  "bbob-biobj_f13_i04_d20 0.997755513117971",
  "bbob-biobj_f13_i04_d40 0.986208111607261",
  "bbob-biobj_f13_i05_d02 0.999818303891124",
  "bbob-biobj_f13_i05_d03 0.999920715948849",
  "bbob-biobj_f13_i05_d05 0.999370136583932",
  "bbob-biobj_f13_i05_d10 0.994188642308816",
  "bbob-biobj_f13_i05_d20 0.998297557674693",
  "bbob-biobj_f13_i05_d40 0.986630965803090",
  "bbob-biobj_f13_i06_d02 0.999785884998514",
  "bbob-biobj_f13_i06_d03 0.999441987331643",
  "bbob-biobj_f13_i06_d05 0.996121377467114",
  "bbob-biobj_f13_i06_d10 0.999309176171936",
  "bbob-biobj_f13_i06_d20 0.971211251821416",
  "bbob-biobj_f13_i06_d40 0.975933691925670",
  "bbob-biobj_f13_i07_d02 0.998280886452156",
  "bbob-biobj_f13_i07_d03 0.999946376863723",
  "bbob-biobj_f13_i07_d05 0.957743677469676",
  "bbob-biobj_f13_i07_d10 0.997017991345758",
  "bbob-biobj_f13_i07_d20 0.986193834743328",
  "bbob-biobj_f13_i07_d40 0.992637426638605",
  "bbob-biobj_f13_i08_d02 0.986243529678102",
  "bbob-biobj_f13_i08_d03 0.987229023477311",
  "bbob-biobj_f13_i08_d05 0.998520022008726",
  "bbob-biobj_f13_i08_d10 0.997546714662799",
  "bbob-biobj_f13_i08_d20 0.992790837432136",
  "bbob-biobj_f13_i08_d40 0.993105523847946",
  "bbob-biobj_f13_i09_d02 0.999837701686230",
  "bbob-biobj_f13_i09_d03 0.999564311440944",
  "bbob-biobj_f13_i09_d05 0.993062622405883",
  "bbob-biobj_f13_i09_d10 0.996828245542617",
  "bbob-biobj_f13_i09_d20 0.998226560559077",
  "bbob-biobj_f13_i09_d40 0.984098714771724",
  "bbob-biobj_f13_i10_d02 0.999981625511090",
  "bbob-biobj_f13_i10_d03 0.998116473112330",
  "bbob-biobj_f13_i10_d05 0.993810109440758",
  "bbob-biobj_f13_i10_d10 0.999182553764442",
  "bbob-biobj_f13_i10_d20 0.999266377199578",
  "bbob-biobj_f13_i10_d40 0.988665896583187",
  "bbob-biobj_f14_i01_d02 0.912471830778958",
  "bbob-biobj_f14_i01_d03 0.832549702621248",
  "bbob-biobj_f14_i01_d05 0.935775796288249",
  "bbob-biobj_f14_i01_d10 0.856110147066621",
  "bbob-biobj_f14_i01_d20 0.838321488493349",
  "bbob-biobj_f14_i01_d40 0.858605910259305",
  "bbob-biobj_f14_i02_d02 0.997754487379014",
  "bbob-biobj_f14_i02_d03 0.996012919660629",
  "bbob-biobj_f14_i02_d05 0.842166972572312",
  "bbob-biobj_f14_i02_d10 0.880375618604676",
  "bbob-biobj_f14_i02_d20 0.909292899014422",
  "bbob-biobj_f14_i02_d40 0.872240987947634",
  "bbob-biobj_f14_i03_d02 0.966432433375896",
  "bbob-biobj_f14_i03_d03 0.983419702052425",
  "bbob-biobj_f14_i03_d05 0.792211190962942",
  "bbob-biobj_f14_i03_d10 0.931283688496531",
  "bbob-biobj_f14_i03_d20 0.918540286100306",
  "bbob-biobj_f14_i03_d40 0.866546462827371",
  "bbob-biobj_f14_i04_d02 0.982899450472848",
  "bbob-biobj_f14_i04_d03 0.999118278368303",
  "bbob-biobj_f14_i04_d05 0.854072197986231",
  "bbob-biobj_f14_i04_d10 0.871663522094536",
  "bbob-biobj_f14_i04_d20 0.902461674417968",
  "bbob-biobj_f14_i04_d40 0.905168523836626",
  "bbob-biobj_f14_i05_d02 0.987027892777520",
  "bbob-biobj_f14_i05_d03 0.986701943293328",
  "bbob-biobj_f14_i05_d05 0.966158837894660",
  "bbob-biobj_f14_i05_d10 0.940516556356812",
  "bbob-biobj_f14_i05_d20 0.885647412872639",
  "bbob-biobj_f14_i05_d40 0.880302038613729",
  "bbob-biobj_f14_i06_d02 0.989610875602541",
  "bbob-biobj_f14_i06_d03 0.967419428698016",
  "bbob-biobj_f14_i06_d05 0.947962195487106",
  "bbob-biobj_f14_i06_d10 0.969663885718599",
  "bbob-biobj_f14_i06_d20 0.877562875248050",
  "bbob-biobj_f14_i06_d40 0.882619053834477",
  "bbob-biobj_f14_i07_d02 0.887837769659544",
  "bbob-biobj_f14_i07_d03 0.994839805885289",
  "bbob-biobj_f14_i07_d05 0.993450701783038",
  "bbob-biobj_f14_i07_d10 0.918295794384147",
  "bbob-biobj_f14_i07_d20 0.845976409577271",
  "bbob-biobj_f14_i07_d40 0.886017282342959",
  "bbob-biobj_f14_i08_d02 0.785251939709812",
  "bbob-biobj_f14_i08_d03 0.811943440895144",
  "bbob-biobj_f14_i08_d05 0.962729480511351",
  "bbob-biobj_f14_i08_d10 0.837638468859587",
  "bbob-biobj_f14_i08_d20 0.864554912097360",
  "bbob-biobj_f14_i08_d40 0.886397089215907",
  "bbob-biobj_f14_i09_d02 0.999594648530017",
  "bbob-biobj_f14_i09_d03 0.835363745633698",
  "bbob-biobj_f14_i09_d05 0.827847936226102",
  "bbob-biobj_f14_i09_d10 0.903105742694003",
  "bbob-biobj_f14_i09_d20 0.899934539970372",
  "bbob-biobj_f14_i09_d40 0.861840932681692",
  "bbob-biobj_f14_i10_d02 0.980317958903754",
  "bbob-biobj_f14_i10_d03 0.870552037123357",
  "bbob-biobj_f14_i10_d05 0.987681201109085",
  "bbob-biobj_f14_i10_d10 0.804330255533714",
  "bbob-biobj_f14_i10_d20 0.948226545458335",
  "bbob-biobj_f14_i10_d40 0.877439564739814",
  "bbob-biobj_f15_i01_d02 0.978329744139448",
  "bbob-biobj_f15_i01_d03 0.928761871647275",
  "bbob-biobj_f15_i01_d05 0.948213253768694",
  "bbob-biobj_f15_i01_d10 0.965007748158915",
  "bbob-biobj_f15_i01_d20 0.988167933478117",
  "bbob-biobj_f15_i01_d40 0.979525774276022",
  "bbob-biobj_f15_i02_d02 0.941436758811989",
  "bbob-biobj_f15_i02_d03 0.954421257030267",
  "bbob-biobj_f15_i02_d05 0.979712604400776",
  "bbob-biobj_f15_i02_d10 0.992133953121318",
  "bbob-biobj_f15_i02_d20 0.957729642814764",
  "bbob-biobj_f15_i02_d40 0.979293415431107",
  "bbob-biobj_f15_i03_d02 0.998672053297861",
  "bbob-biobj_f15_i03_d03 0.994036318703754",
  "bbob-biobj_f15_i03_d05 0.980683737423164",
  "bbob-biobj_f15_i03_d10 0.991361399690812",
  "bbob-biobj_f15_i03_d20 0.982818016023940",
  "bbob-biobj_f15_i03_d40 0.980570309277482",
  "bbob-biobj_f15_i04_d02 0.825078978774112",
  "bbob-biobj_f15_i04_d03 0.975014462586473",
  "bbob-biobj_f15_i04_d05 0.944759757043127",
  "bbob-biobj_f15_i04_d10 0.994604317225863",
  "bbob-biobj_f15_i04_d20 0.973598426097798",
  "bbob-biobj_f15_i04_d40 0.993777255236331",
  "bbob-biobj_f15_i05_d02 0.988337532739167",
  "bbob-biobj_f15_i05_d03 0.932552115048332",
  "bbob-biobj_f15_i05_d05 0.867415698494500",
  "bbob-biobj_f15_i05_d10 0.985324596574554",
  "bbob-biobj_f15_i05_d20 0.988079063451025",
  "bbob-biobj_f15_i05_d40 0.992416053597612",
  "bbob-biobj_f15_i06_d02 0.969085878571844",
  "bbob-biobj_f15_i06_d03 0.955382906442424",
  "bbob-biobj_f15_i06_d05 0.975599518082786",
  "bbob-biobj_f15_i06_d10 0.985307869919625",
  "bbob-biobj_f15_i06_d20 0.985003970094727",
  "bbob-biobj_f15_i06_d40 0.977089099269068",
  "bbob-biobj_f15_i07_d02 0.982625119376436",
  "bbob-biobj_f15_i07_d03 0.904495724916453",
  "bbob-biobj_f15_i07_d05 0.993506742600270",
  "bbob-biobj_f15_i07_d10 0.956841407597722",
  "bbob-biobj_f15_i07_d20 0.977280331164342",
  "bbob-biobj_f15_i07_d40 0.982248710216346",
  "bbob-biobj_f15_i08_d02 0.925488457367322",
  "bbob-biobj_f15_i08_d03 0.886626169500458",
  "bbob-biobj_f15_i08_d05 0.977518440597604",
  "bbob-biobj_f15_i08_d10 0.987522911600190",
  "bbob-biobj_f15_i08_d20 0.972961886598506",
  "bbob-biobj_f15_i08_d40 0.987668544194751",
  "bbob-biobj_f15_i09_d02 0.985671370725154",
  "bbob-biobj_f15_i09_d03 0.988738880947552",
  "bbob-biobj_f15_i09_d05 0.932479987870587",
  "bbob-biobj_f15_i09_d10 0.974229585800356",
  "bbob-biobj_f15_i09_d20 0.951562817077932",
  "bbob-biobj_f15_i09_d40 0.980888112549418",
  "bbob-biobj_f15_i10_d02 0.932445823916236",
  "bbob-biobj_f15_i10_d03 0.760223252492131",
  "bbob-biobj_f15_i10_d05 0.962365594477710",
  "bbob-biobj_f15_i10_d10 0.968501543202303",
  "bbob-biobj_f15_i10_d20 0.983957139266508",
  "bbob-biobj_f15_i10_d40 0.989438827925627",
  "bbob-biobj_f16_i01_d02 0.981136700654277",
  "bbob-biobj_f16_i01_d03 0.999116510727310",
  "bbob-biobj_f16_i01_d05 0.991904924593850",
  "bbob-biobj_f16_i01_d10 0.985347609776800",
  "bbob-biobj_f16_i01_d20 0.974971231310992",
  "bbob-biobj_f16_i01_d40 0.956921741797102",
  "bbob-biobj_f16_i02_d02 0.965606135498190",
  "bbob-biobj_f16_i02_d03 0.934170201491467",
  "bbob-biobj_f16_i02_d05 0.954677025450060",
  "bbob-biobj_f16_i02_d10 0.960199586818040",
  "bbob-biobj_f16_i02_d20 0.977149756950747",
  "bbob-biobj_f16_i02_d40 0.957646353413762",
  "bbob-biobj_f16_i03_d02 0.959401149990272",
  "bbob-biobj_f16_i03_d03 0.943744234862334",
  "bbob-biobj_f16_i03_d05 0.990029096692041",
  "bbob-biobj_f16_i03_d10 0.968399176349063",
  "bbob-biobj_f16_i03_d20 0.987412412519172",
  "bbob-biobj_f16_i03_d40 0.970592596401184",
  "bbob-biobj_f16_i04_d02 0.995345563313723",
  "bbob-biobj_f16_i04_d03 0.997825549300744",
  "bbob-biobj_f16_i04_d05 0.984584481171368",
  "bbob-biobj_f16_i04_d10 0.984721586998123",
  "bbob-biobj_f16_i04_d20 0.979094578910330",
  "bbob-biobj_f16_i04_d40 0.966636363162328",
  "bbob-biobj_f16_i05_d02 0.999043228368928",
  "bbob-biobj_f16_i05_d03 0.953342518024141",
  "bbob-biobj_f16_i05_d05 0.995251101113935",
  "bbob-biobj_f16_i05_d10 0.974136855929559",
  "bbob-biobj_f16_i05_d20 0.977421249413489",
  "bbob-biobj_f16_i05_d40 0.981290213354151",
  "bbob-biobj_f16_i06_d02 0.139712205427015",
  "bbob-biobj_f16_i06_d03 0.931037964250559",
  "bbob-biobj_f16_i06_d05 0.948669437105361",
  "bbob-biobj_f16_i06_d10 0.977958207843702",
  "bbob-biobj_f16_i06_d20 0.982494249192594",
  "bbob-biobj_f16_i06_d40 0.936635770451486",
  "bbob-biobj_f16_i07_d02 0.997707544883722",
  "bbob-biobj_f16_i07_d03 0.953474102846906",
  "bbob-biobj_f16_i07_d05 0.991371726429119",
  "bbob-biobj_f16_i07_d10 0.951308997934433",
  "bbob-biobj_f16_i07_d20 0.970433007701376",
  "bbob-biobj_f16_i07_d40 0.946969469111869",
  "bbob-biobj_f16_i08_d02 0.968053065624156",
  "bbob-biobj_f16_i08_d03 0.984027399418469",
  "bbob-biobj_f16_i08_d05 0.970406708176693",
  "bbob-biobj_f16_i08_d10 0.981887161503345",
  "bbob-biobj_f16_i08_d20 0.991067110914981",
  "bbob-biobj_f16_i08_d40 0.906801599239074",
  "bbob-biobj_f16_i09_d02 0.971189887986641",
  "bbob-biobj_f16_i09_d03 0.970827719065357",
  "bbob-biobj_f16_i09_d05 0.969883088508387",
  "bbob-biobj_f16_i09_d10 0.992543440275702",
  "bbob-biobj_f16_i09_d20 0.979320895550577",
  "bbob-biobj_f16_i09_d40 0.964134468147611",
  "bbob-biobj_f16_i10_d02 0.928205956144776",
  "bbob-biobj_f16_i10_d03 0.991573103524712",
  "bbob-biobj_f16_i10_d05 0.947646759618908",
  "bbob-biobj_f16_i10_d10 0.993479223551865",
  "bbob-biobj_f16_i10_d20 0.985344436021621",
  "bbob-biobj_f16_i10_d40 0.949572615434425",
  "bbob-biobj_f17_i01_d02 0.979958696582261",
  "bbob-biobj_f17_i01_d03 0.899961483197216",
  "bbob-biobj_f17_i01_d05 0.980212408948254",
  "bbob-biobj_f17_i01_d10 0.991499141850806",
  "bbob-biobj_f17_i01_d20 0.973803574766593",
  "bbob-biobj_f17_i01_d40 0.977329568943416",
  "bbob-biobj_f17_i02_d02 0.942990097100050",
  "bbob-biobj_f17_i02_d03 0.931420423215006",
  "bbob-biobj_f17_i02_d05 0.971547274550899",
  "bbob-biobj_f17_i02_d10 0.981775672319523",
  "bbob-biobj_f17_i02_d20 0.968195623409150",
  "bbob-biobj_f17_i02_d40 0.977915877795807",
  "bbob-biobj_f17_i03_d02 0.941726842813370",
  "bbob-biobj_f17_i03_d03 0.956408723394532",
  "bbob-biobj_f17_i03_d05 0.967807670314565",
  "bbob-biobj_f17_i03_d10 0.990814719626047",
  "bbob-biobj_f17_i03_d20 0.985631075660359",
  "bbob-biobj_f17_i03_d40 0.955265611433480",
  "bbob-biobj_f17_i04_d02 0.766767384446464",
  "bbob-biobj_f17_i04_d03 0.861029503563583",
  "bbob-biobj_f17_i04_d05 0.995625749291626",
  "bbob-biobj_f17_i04_d10 0.993472114071952",
  "bbob-biobj_f17_i04_d20 0.993745049323583",
  "bbob-biobj_f17_i04_d40 0.985588911394250",
  "bbob-biobj_f17_i05_d02 0.940686728497900",
  "bbob-biobj_f17_i05_d03 0.855772681754819",
  "bbob-biobj_f17_i05_d05 0.989638126261116",
  "bbob-biobj_f17_i05_d10 0.988021171884147",
  "bbob-biobj_f17_i05_d20 0.998300078076616",
  "bbob-biobj_f17_i05_d40 0.970726773312885",
  "bbob-biobj_f17_i06_d02 0.987713985194397",
  "bbob-biobj_f17_i06_d03 0.933596473754802",
  "bbob-biobj_f17_i06_d05 0.993060463789679",
  "bbob-biobj_f17_i06_d10 0.986631085585597",
  "bbob-biobj_f17_i06_d20 0.984883683291481",
  "bbob-biobj_f17_i06_d40 0.929787924268597",
  "bbob-biobj_f17_i07_d02 0.996487340351512",
  "bbob-biobj_f17_i07_d03 0.953930870414184",
  "bbob-biobj_f17_i07_d05 0.995823634879168",
  "bbob-biobj_f17_i07_d10 0.993081757419179",
  "bbob-biobj_f17_i07_d20 0.998770776346561",
  "bbob-biobj_f17_i07_d40 0.960616338166501",
  "bbob-biobj_f17_i08_d02 0.889841257972338",
  "bbob-biobj_f17_i08_d03 0.823820493238124",
  "bbob-biobj_f17_i08_d05 0.906608013582778",
  "bbob-biobj_f17_i08_d10 0.991080354922109",
  "bbob-biobj_f17_i08_d20 0.980730037656149",
  "bbob-biobj_f17_i08_d40 0.946053183664007",
  "bbob-biobj_f17_i09_d02 0.939017725322369",
  "bbob-biobj_f17_i09_d03 0.984693732150837",
  "bbob-biobj_f17_i09_d05 0.939456807354725",
  "bbob-biobj_f17_i09_d10 0.990885474758620",
  "bbob-biobj_f17_i09_d20 0.993408423512887",
  "bbob-biobj_f17_i09_d40 0.918142489996231",
  "bbob-biobj_f17_i10_d02 0.998455044767241",
  "bbob-biobj_f17_i10_d03 0.952850070908516",
  "bbob-biobj_f17_i10_d05 0.976079006404020",
  "bbob-biobj_f17_i10_d10 0.994011795242328",
  "bbob-biobj_f17_i10_d20 0.985724108347979",
  "bbob-biobj_f17_i10_d40 0.945974036468278",
  "bbob-biobj_f18_i01_d02 0.969204867140222",
  "bbob-biobj_f18_i01_d03 0.998688622271819",
  "bbob-biobj_f18_i01_d05 0.998492819648938",
  "bbob-biobj_f18_i01_d10 0.992629645389604",
  "bbob-biobj_f18_i01_d20 0.947912375007042",
  "bbob-biobj_f18_i01_d40 0.946337627383827",
  "bbob-biobj_f18_i02_d02 0.953409533896478",
  "bbob-biobj_f18_i02_d03 0.993777473977030",
  "bbob-biobj_f18_i02_d05 0.991742402096738",
  "bbob-biobj_f18_i02_d10 0.962762745769322",
  "bbob-biobj_f18_i02_d20 0.950134125189844",
  "bbob-biobj_f18_i02_d40 0.985121020300274",
  "bbob-biobj_f18_i03_d02 0.999403466523437",
  "bbob-biobj_f18_i03_d03 0.999304353432177",
  "bbob-biobj_f18_i03_d05 0.946581786017729",
  "bbob-biobj_f18_i03_d10 0.992474105619029",
  "bbob-biobj_f18_i03_d20 0.960658596619811",
  "bbob-biobj_f18_i03_d40 0.977564674019477",
  "bbob-biobj_f18_i04_d02 0.990014933357924",
  "bbob-biobj_f18_i04_d03 0.998331130733417",
  "bbob-biobj_f18_i04_d05 0.973768285471165",
  "bbob-biobj_f18_i04_d10 0.950023687137867",
  "bbob-biobj_f18_i04_d20 0.933587248934312",
  "bbob-biobj_f18_i04_d40 0.982988682493896",
  "bbob-biobj_f18_i05_d02 0.999935347951619",
  "bbob-biobj_f18_i05_d03 0.975931617728692",
  "bbob-biobj_f18_i05_d05 0.941220553067101",
  "bbob-biobj_f18_i05_d10 0.954459542216504",
  "bbob-biobj_f18_i05_d20 0.975441006994456",
  "bbob-biobj_f18_i05_d40 0.972422258828605",
  "bbob-biobj_f18_i06_d02 0.999890514137715",
  "bbob-biobj_f18_i06_d03 0.951926054626427",
  "bbob-biobj_f18_i06_d05 0.972128197216178",
  "bbob-biobj_f18_i06_d10 0.983782771464663",
  "bbob-biobj_f18_i06_d20 0.988454060894957",
  "bbob-biobj_f18_i06_d40 0.943709125442017",
  "bbob-biobj_f18_i07_d02 0.999722776570118",
  "bbob-biobj_f18_i07_d03 0.939207716362455",
  "bbob-biobj_f18_i07_d05 0.973384396295115",
  "bbob-biobj_f18_i07_d10 0.956602236799538",
  "bbob-biobj_f18_i07_d20 0.975080974637268",
  "bbob-biobj_f18_i07_d40 0.963040285384617",
  "bbob-biobj_f18_i08_d02 0.979510876330516",
  "bbob-biobj_f18_i08_d03 0.983564503887299",
  "bbob-biobj_f18_i08_d05 0.972260953789183",
  "bbob-biobj_f18_i08_d10 0.943588491564128",
  "bbob-biobj_f18_i08_d20 0.984980809657981",
  "bbob-biobj_f18_i08_d40 0.967628799027729",
  "bbob-biobj_f18_i09_d02 0.826509241045423",
  "bbob-biobj_f18_i09_d03 0.955224693894568",
  "bbob-biobj_f18_i09_d05 0.990329407682473",
  "bbob-biobj_f18_i09_d10 0.973673240837011",
  "bbob-biobj_f18_i09_d20 0.927149821301475",
  "bbob-biobj_f18_i09_d40 0.940546921139129",
  "bbob-biobj_f18_i10_d02 0.981156284083054",
  "bbob-biobj_f18_i10_d03 0.977461510507809",
  "bbob-biobj_f18_i10_d05 0.973048715805323",
  "bbob-biobj_f18_i10_d10 0.930460976380332",
  "bbob-biobj_f18_i10_d20 0.958165972527928",
  "bbob-biobj_f18_i10_d40 0.936488898975624",
  "bbob-biobj_f19_i01_d02 0.865811588625181",
  "bbob-biobj_f19_i01_d03 0.973016869560970",
  "bbob-biobj_f19_i01_d05 0.992521219434583",
  "bbob-biobj_f19_i01_d10 0.992012342550646",
  "bbob-biobj_f19_i01_d20 0.985705751927878",
  "bbob-biobj_f19_i01_d40 0.968412445641202",
  "bbob-biobj_f19_i02_d02 0.920899376364515",
  "bbob-biobj_f19_i02_d03 0.986110860778247",
  "bbob-biobj_f19_i02_d05 0.989653718047710",
  "bbob-biobj_f19_i02_d10 0.998068968304825",
  "bbob-biobj_f19_i02_d20 0.971805630202519",
  "bbob-biobj_f19_i02_d40 0.947182109899931",
  "bbob-biobj_f19_i03_d02 0.904545240983088",
  "bbob-biobj_f19_i03_d03 0.957586890418128",
  "bbob-biobj_f19_i03_d05 0.982448581261593",
  "bbob-biobj_f19_i03_d10 0.991728876050119",
  "bbob-biobj_f19_i03_d20 0.991085360226958",
  "bbob-biobj_f19_i03_d40 0.962863247619966",
  "bbob-biobj_f19_i04_d02 0.999884820914529",
  "bbob-biobj_f19_i04_d03 0.996380904635380",
  "bbob-biobj_f19_i04_d05 0.992031336531309",
  "bbob-biobj_f19_i04_d10 0.992957066422246",
  "bbob-biobj_f19_i04_d20 0.976691229613891",
  "bbob-biobj_f19_i04_d40 0.954311507194530",
  "bbob-biobj_f19_i05_d02 0.997258704039593",
  "bbob-biobj_f19_i05_d03 0.959361638755915",
  "bbob-biobj_f19_i05_d05 0.993758195885880",
  "bbob-biobj_f19_i05_d10 0.992125940103881",
  "bbob-biobj_f19_i05_d20 0.984452338080095",
  "bbob-biobj_f19_i05_d40 0.938647177698134",
  "bbob-biobj_f19_i06_d02 0.955179571170891",
  "bbob-biobj_f19_i06_d03 0.991784618477664",
  "bbob-biobj_f19_i06_d05 0.995285000486707",
  "bbob-biobj_f19_i06_d10 0.993449923825602",
  "bbob-biobj_f19_i06_d20 0.988866349458714",
  "bbob-biobj_f19_i06_d40 0.926458477834549",
  "bbob-biobj_f19_i07_d02 0.916731595171689",
  "bbob-biobj_f19_i07_d03 0.976691255967938",
  "bbob-biobj_f19_i07_d05 0.998883278152017",
  "bbob-biobj_f19_i07_d10 0.984710607481539",
  "bbob-biobj_f19_i07_d20 0.983632123694558",
  "bbob-biobj_f19_i07_d40 0.977504251577860",
  "bbob-biobj_f19_i08_d02 0.954018164994102",
  "bbob-biobj_f19_i08_d03 0.958826344007053",
  "bbob-biobj_f19_i08_d05 0.980883278818685",
  "bbob-biobj_f19_i08_d10 0.988434270645912",
  "bbob-biobj_f19_i08_d20 0.975025754015807",
  "bbob-biobj_f19_i08_d40 0.931399074334377",
  "bbob-biobj_f19_i09_d02 0.949449196092764",
  "bbob-biobj_f19_i09_d03 0.968277545349826",
  "bbob-biobj_f19_i09_d05 0.955094713161426",
  "bbob-biobj_f19_i09_d10 0.990769409987048",
  "bbob-biobj_f19_i09_d20 0.980825420068253",
  "bbob-biobj_f19_i09_d40 0.961450631623366",
  "bbob-biobj_f19_i10_d02 0.946053790931798",
  "bbob-biobj_f19_i10_d03 0.980549229995464",
  "bbob-biobj_f19_i10_d05 0.999353208304047",
  "bbob-biobj_f19_i10_d10 0.987379497096524",
  "bbob-biobj_f19_i10_d20 0.988854714616387",
  "bbob-biobj_f19_i10_d40 0.896565493720926",
  "bbob-biobj_f20_i01_d02 0.995717988201693",
  "bbob-biobj_f20_i01_d03 0.813226741480784",
  "bbob-biobj_f20_i01_d05 0.967868469380036",
  "bbob-biobj_f20_i01_d10 0.999902663760163",
  "bbob-biobj_f20_i01_d20 0.999771502228653",
  "bbob-biobj_f20_i01_d40 0.999489958018389",
  "bbob-biobj_f20_i02_d02 0.910959916528987",
  "bbob-biobj_f20_i02_d03 0.974016949278183",
  "bbob-biobj_f20_i02_d05 0.813996181562145",
  "bbob-biobj_f20_i02_d10 0.999676294376271",
  "bbob-biobj_f20_i02_d20 0.999651451997351",
  "bbob-biobj_f20_i02_d40 0.999733893382635",
  "bbob-biobj_f20_i03_d02 0.985777664983236",
  "bbob-biobj_f20_i03_d03 0.999442770435194",
  "bbob-biobj_f20_i03_d05 0.892401200948616",
  "bbob-biobj_f20_i03_d10 0.972681061470658",
  "bbob-biobj_f20_i03_d20 0.999807557922657",
  "bbob-biobj_f20_i03_d40 0.996958424990987",
  "bbob-biobj_f20_i04_d02 0.890179820278766",
  "bbob-biobj_f20_i04_d03 0.925391276015625",
  "bbob-biobj_f20_i04_d05 0.999111453517483",
  "bbob-biobj_f20_i04_d10 0.999930922553685",
  "bbob-biobj_f20_i04_d20 0.999628162387754",
  "bbob-biobj_f20_i04_d40 0.999868192453144",
  "bbob-biobj_f20_i05_d02 0.980549187818064",
  "bbob-biobj_f20_i05_d03 0.963674542281100",
  "bbob-biobj_f20_i05_d05 0.946993788339496",
  "bbob-biobj_f20_i05_d10 0.999846387265367",
  "bbob-biobj_f20_i05_d20 0.958581965640984",
  "bbob-biobj_f20_i05_d40 0.999726220269249",
  "bbob-biobj_f20_i06_d02 0.868960770544449",
  "bbob-biobj_f20_i06_d03 0.999652403872231",
  "bbob-biobj_f20_i06_d05 0.888518523925641",
  "bbob-biobj_f20_i06_d10 0.921048822325153",
  "bbob-biobj_f20_i06_d20 0.999864117303474",
  "bbob-biobj_f20_i06_d40 0.996969194081855",
  "bbob-biobj_f20_i07_d02 0.950863247036010",
  "bbob-biobj_f20_i07_d03 0.864246023092120",
  "bbob-biobj_f20_i07_d05 0.999629856094957",
  "bbob-biobj_f20_i07_d10 0.972922950949932",
  "bbob-biobj_f20_i07_d20 0.982938280870815",
  "bbob-biobj_f20_i07_d40 0.989945773024964",
  "bbob-biobj_f20_i08_d02 0.840670532011553",
  "bbob-biobj_f20_i08_d03 0.849546600060598",
  "bbob-biobj_f20_i08_d05 0.994125044141308",
  "bbob-biobj_f20_i08_d10 0.999722234140680",
  "bbob-biobj_f20_i08_d20 0.999799187932440",
  "bbob-biobj_f20_i08_d40 0.997154124092010",
  "bbob-biobj_f20_i09_d02 0.998335552575094",
  "bbob-biobj_f20_i09_d03 0.999065130079023",
  "bbob-biobj_f20_i09_d05 0.963412153353338",
  "bbob-biobj_f20_i09_d10 0.999432425503625",
  "bbob-biobj_f20_i09_d20 0.999759178540532",
  "bbob-biobj_f20_i09_d40 0.999794888401685",
  "bbob-biobj_f20_i10_d02 0.994830610836821",
  "bbob-biobj_f20_i10_d03 0.945475716626870",
  "bbob-biobj_f20_i10_d05 0.859699717321543",
  "bbob-biobj_f20_i10_d10 0.999802870884534",
  "bbob-biobj_f20_i10_d20 0.999303800757148",
  "bbob-biobj_f20_i10_d40 0.992311894295089",
  "bbob-biobj_f21_i01_d02 0.999736611084103",
  "bbob-biobj_f21_i01_d03 0.911539551415472",
  "bbob-biobj_f21_i01_d05 0.912540567945498",
  "bbob-biobj_f21_i01_d10 0.993978538692239",
  "bbob-biobj_f21_i01_d20 0.981406877967379",
  "bbob-biobj_f21_i01_d40 0.983499264936933",
  "bbob-biobj_f21_i02_d02 0.985471376910841",
  "bbob-biobj_f21_i02_d03 0.980784048920377",
  "bbob-biobj_f21_i02_d05 0.998160286632846",
  "bbob-biobj_f21_i02_d10 0.996909063057436",
  "bbob-biobj_f21_i02_d20 0.993493283068370",
  "bbob-biobj_f21_i02_d40 0.996169084540127",
  "bbob-biobj_f21_i03_d02 0.973889757463946",
  "bbob-biobj_f21_i03_d03 0.968953532443341",
  "bbob-biobj_f21_i03_d05 0.929903444365651",
  "bbob-biobj_f21_i03_d10 0.998452863451952",
  "bbob-biobj_f21_i03_d20 0.993060764992102",
  "bbob-biobj_f21_i03_d40 0.995895079588410",
  "bbob-biobj_f21_i04_d02 0.999788748975104",
  "bbob-biobj_f21_i04_d03 0.954073986608832",
  "bbob-biobj_f21_i04_d05 0.928301923789326",
  "bbob-biobj_f21_i04_d10 0.904321643215245",
  "bbob-biobj_f21_i04_d20 0.996038359624565",
  "bbob-biobj_f21_i04_d40 0.968247897361166",
  "bbob-biobj_f21_i05_d02 0.890724797691638",
  "bbob-biobj_f21_i05_d03 0.999573344629867",
  "bbob-biobj_f21_i05_d05 0.997869894633878",
  "bbob-biobj_f21_i05_d10 0.958079232794853",
  "bbob-biobj_f21_i05_d20 0.982086183547206",
  "bbob-biobj_f21_i05_d40 0.985493219192335",
  "bbob-biobj_f21_i06_d02 0.998904552168893",
  "bbob-biobj_f21_i06_d03 0.898627121689058",
  "bbob-biobj_f21_i06_d05 0.999919475038921",
  "bbob-biobj_f21_i06_d10 0.993107312405050",
  "bbob-biobj_f21_i06_d20 0.995669286497430",
  "bbob-biobj_f21_i06_d40 0.990788892044220",
  "bbob-biobj_f21_i07_d02 0.862039447736370",
  "bbob-biobj_f21_i07_d03 0.994758030627227",
  "bbob-biobj_f21_i07_d05 0.999834546461715",
  "bbob-biobj_f21_i07_d10 0.997655833365488",
  "bbob-biobj_f21_i07_d20 0.953933077250989",
  "bbob-biobj_f21_i07_d40 0.977666056653005",
  "bbob-biobj_f21_i08_d02 0.986052380752108",
  "bbob-biobj_f21_i08_d03 0.999625266184618",
  "bbob-biobj_f21_i08_d05 0.942302148770571",
  "bbob-biobj_f21_i08_d10 0.996894931837987",
  "bbob-biobj_f21_i08_d20 0.988122140029128",
  "bbob-biobj_f21_i08_d40 0.997775728057085",
  "bbob-biobj_f21_i09_d02 0.972462187334999",
  "bbob-biobj_f21_i09_d03 0.998497439307589",
  "bbob-biobj_f21_i09_d05 0.998950788096672",
  "bbob-biobj_f21_i09_d10 0.981086470119500",
  "bbob-biobj_f21_i09_d20 0.981962233114373",
  "bbob-biobj_f21_i09_d40 0.989294182298165",
  "bbob-biobj_f21_i10_d02 0.940849980475781",
  "bbob-biobj_f21_i10_d03 0.949332925002064",
  "bbob-biobj_f21_i10_d05 0.979555325255604",
  "bbob-biobj_f21_i10_d10 0.977252537271692",
  "bbob-biobj_f21_i10_d20 0.972947454557864",
  "bbob-biobj_f21_i10_d40 0.997031890987017",
  "bbob-biobj_f22_i01_d02 0.700915214438678",
  "bbob-biobj_f22_i01_d03 0.694480395902342",
  "bbob-biobj_f22_i01_d05 0.986971595315615",
  "bbob-biobj_f22_i01_d10 0.837890262731969",
  "bbob-biobj_f22_i01_d20 0.771210398849058",
  "bbob-biobj_f22_i01_d40 0.795334918259521",
  "bbob-biobj_f22_i02_d02 0.999055910320585",
  "bbob-biobj_f22_i02_d03 0.742599610984376",
  "bbob-biobj_f22_i02_d05 0.764268896053655",
  "bbob-biobj_f22_i02_d10 0.728735699616278",
  "bbob-biobj_f22_i02_d20 0.756653207913305",
  "bbob-biobj_f22_i02_d40 0.854446084484735",
  "bbob-biobj_f22_i03_d02 0.678522922384644",
  "bbob-biobj_f22_i03_d03 0.951234163557838",
  "bbob-biobj_f22_i03_d05 0.735416447436036",
  "bbob-biobj_f22_i03_d10 0.862571148987637",
  "bbob-biobj_f22_i03_d20 0.863107380943036",
  "bbob-biobj_f22_i03_d40 0.774848997128731",
  "bbob-biobj_f22_i04_d02 0.846365044335506",
  "bbob-biobj_f22_i04_d03 0.803804267073623",
  "bbob-biobj_f22_i04_d05 0.834636451792339",
  "bbob-biobj_f22_i04_d10 0.842296787779240",
  "bbob-biobj_f22_i04_d20 0.915870982856891",
  "bbob-biobj_f22_i04_d40 0.804689191550325",
  "bbob-biobj_f22_i05_d02 0.856038328100935",
  "bbob-biobj_f22_i05_d03 0.929852773705256",
  "bbob-biobj_f22_i05_d05 0.892867161779887",
  "bbob-biobj_f22_i05_d10 0.819615576069317",
  "bbob-biobj_f22_i05_d20 0.789423786912445",
  "bbob-biobj_f22_i05_d40 0.763549979203388",
  "bbob-biobj_f22_i06_d02 0.977959593880942",
  "bbob-biobj_f22_i06_d03 0.724307485265540",
  "bbob-biobj_f22_i06_d05 0.811223433694945",
  "bbob-biobj_f22_i06_d10 0.861024567753003",
  "bbob-biobj_f22_i06_d20 0.761732191531536",
  "bbob-biobj_f22_i06_d40 0.890315126237409",
  "bbob-biobj_f22_i07_d02 0.909936185125206",
  "bbob-biobj_f22_i07_d03 0.690658544562536",
  "bbob-biobj_f22_i07_d05 0.722187169778075",
  "bbob-biobj_f22_i07_d10 0.786251401254339",
  "bbob-biobj_f22_i07_d20 0.763669328950710",
  "bbob-biobj_f22_i07_d40 0.770601878006756",
  "bbob-biobj_f22_i08_d02 0.906823592358557",
  "bbob-biobj_f22_i08_d03 0.835715891680142",
  "bbob-biobj_f22_i08_d05 0.771934815309923",
  "bbob-biobj_f22_i08_d10 0.940052340139338",
  "bbob-biobj_f22_i08_d20 0.769308735348968",
  "bbob-biobj_f22_i08_d40 0.783637202717775",
  "bbob-biobj_f22_i09_d02 0.968955004271021",
  "bbob-biobj_f22_i09_d03 0.916880432118617",
  "bbob-biobj_f22_i09_d05 0.950105977970102",
  "bbob-biobj_f22_i09_d10 0.899068097141043",
  "bbob-biobj_f22_i09_d20 0.812635898363424",
  "bbob-biobj_f22_i09_d40 0.825404007682545",
  "bbob-biobj_f22_i10_d02 0.928374894313840",
  "bbob-biobj_f22_i10_d03 0.652282247044987",
  "bbob-biobj_f22_i10_d05 0.761309421584228",
  "bbob-biobj_f22_i10_d10 0.743016739860994",
  "bbob-biobj_f22_i10_d20 0.700240282578198",
  "bbob-biobj_f22_i10_d40 0.822666948073165",
  "bbob-biobj_f23_i01_d02 0.992533953395154",
  "bbob-biobj_f23_i01_d03 0.872234688902190",
  "bbob-biobj_f23_i01_d05 0.980181632277668",
  "bbob-biobj_f23_i01_d10 0.991306818909370",
  "bbob-biobj_f23_i01_d20 0.941113208773929",
  "bbob-biobj_f23_i01_d40 0.968950422821488",
  "bbob-biobj_f23_i02_d02 0.996445249589204",
  "bbob-biobj_f23_i02_d03 0.862093973604522",
  "bbob-biobj_f23_i02_d05 0.919716165222663",
  "bbob-biobj_f23_i02_d10 0.957100598517084",
  "bbob-biobj_f23_i02_d20 0.916201049788741",
  "bbob-biobj_f23_i02_d40 0.975995263825478",
  "bbob-biobj_f23_i03_d02 0.928099254530052",
  "bbob-biobj_f23_i03_d03 0.880342180563437",
  "bbob-biobj_f23_i03_d05 0.894398815403811",
  "bbob-biobj_f23_i03_d10 0.985387963274745",
  "bbob-biobj_f23_i03_d20 0.969060786370051",
  "bbob-biobj_f23_i03_d40 0.978431958926778",
  "bbob-biobj_f23_i04_d02 0.965417925195475",
  "bbob-biobj_f23_i04_d03 0.925755922761600",
  "bbob-biobj_f23_i04_d05 0.887971184395093",
  "bbob-biobj_f23_i04_d10 0.980012506333671",
  "bbob-biobj_f23_i04_d20 0.983093267579092",
  "bbob-biobj_f23_i04_d40 0.970220673432683",
  "bbob-biobj_f23_i05_d02 0.738908860896508",
  "bbob-biobj_f23_i05_d03 0.840964704668912",
  "bbob-biobj_f23_i05_d05 0.891337667840918",
  "bbob-biobj_f23_i05_d10 0.954453677708731",
  "bbob-biobj_f23_i05_d20 0.958512795338489",
  "bbob-biobj_f23_i05_d40 0.951750674629426",
  "bbob-biobj_f23_i06_d02 0.897572283126700",
  "bbob-biobj_f23_i06_d03 0.850273341814831",
  "bbob-biobj_f23_i06_d05 0.902791204372538",
  "bbob-biobj_f23_i06_d10 0.953859913756311",
  "bbob-biobj_f23_i06_d20 0.961791809256483",
  "bbob-biobj_f23_i06_d40 0.973433862523446",
  "bbob-biobj_f23_i07_d02 0.911688067688116",
  "bbob-biobj_f23_i07_d03 0.945019575684005",
  "bbob-biobj_f23_i07_d05 0.984510821980045",
  "bbob-biobj_f23_i07_d10 0.892785278602108",
  "bbob-biobj_f23_i07_d20 0.885709523616749",
  "bbob-biobj_f23_i07_d40 0.957374915872667",
  "bbob-biobj_f23_i08_d02 0.980393684736155",
  "bbob-biobj_f23_i08_d03 0.940130600749823",
  "bbob-biobj_f23_i08_d05 0.912641230372312",
  "bbob-biobj_f23_i08_d10 0.957177567359608",
  "bbob-biobj_f23_i08_d20 0.968121334827907",
  "bbob-biobj_f23_i08_d40 0.974012325493684",
  "bbob-biobj_f23_i09_d02 0.965437295935279",
  "bbob-biobj_f23_i09_d03 0.949151881326183",
  "bbob-biobj_f23_i09_d05 0.975209388126653",
  "bbob-biobj_f23_i09_d10 0.957139500910925",
  "bbob-biobj_f23_i09_d20 0.930214045145034",
  "bbob-biobj_f23_i09_d40 0.971823376224681",
  "bbob-biobj_f23_i10_d02 0.965968828591755",
  "bbob-biobj_f23_i10_d03 0.935147679436830",
  "bbob-biobj_f23_i10_d05 0.842741197682145",
  "bbob-biobj_f23_i10_d10 0.926175701991374",
  "bbob-biobj_f23_i10_d20 0.896674707003784",
  "bbob-biobj_f23_i10_d40 0.968907787674549",
  "bbob-biobj_f24_i01_d02 0.886262905923173",
  "bbob-biobj_f24_i01_d03 0.987113955696040",
  "bbob-biobj_f24_i01_d05 0.869771497008233",
  "bbob-biobj_f24_i01_d10 0.954128581545177",
  "bbob-biobj_f24_i01_d20 0.943099103083434",
  "bbob-biobj_f24_i01_d40 0.946967357198072",
  "bbob-biobj_f24_i02_d02 0.988139628830382",
  "bbob-biobj_f24_i02_d03 0.906306287240695",
  "bbob-biobj_f24_i02_d05 0.947908084579683",
  "bbob-biobj_f24_i02_d10 0.896877403629461",
  "bbob-biobj_f24_i02_d20 0.986633685158978",
  "bbob-biobj_f24_i02_d40 0.966496587220493",
  "bbob-biobj_f24_i03_d02 0.952811572384392",
  "bbob-biobj_f24_i03_d03 0.855065248383182",
  "bbob-biobj_f24_i03_d05 0.935400761731271",
  "bbob-biobj_f24_i03_d10 0.909965560486538",
  "bbob-biobj_f24_i03_d20 0.981205469113305",
  "bbob-biobj_f24_i03_d40 0.932568910096235",
  "bbob-biobj_f24_i04_d02 0.953732218651997",
  "bbob-biobj_f24_i04_d03 0.931779822752282",
  "bbob-biobj_f24_i04_d05 0.880902711714886",
  "bbob-biobj_f24_i04_d10 0.968604785778216",
  "bbob-biobj_f24_i04_d20 0.949134746316553",
  "bbob-biobj_f24_i04_d40 0.941456056510291",
  "bbob-biobj_f24_i05_d02 0.823876907014444",
  "bbob-biobj_f24_i05_d03 0.892822575129825",
  "bbob-biobj_f24_i05_d05 0.958073507719365",
  "bbob-biobj_f24_i05_d10 0.939524316051582",
  "bbob-biobj_f24_i05_d20 0.926050594381492",
  "bbob-biobj_f24_i05_d40 0.948938112360789",
  "bbob-biobj_f24_i06_d02 0.961536381747650",
  "bbob-biobj_f24_i06_d03 0.992941283427474",
  "bbob-biobj_f24_i06_d05 0.866140818108140",
  "bbob-biobj_f24_i06_d10 0.954002012140290",
  "bbob-biobj_f24_i06_d20 0.968423660145754",
  "bbob-biobj_f24_i06_d40 0.966375402700913",
  "bbob-biobj_f24_i07_d02 0.880960114309728",
  "bbob-biobj_f24_i07_d03 0.926062674039788",
  "bbob-biobj_f24_i07_d05 0.949558951423851",
  "bbob-biobj_f24_i07_d10 0.927700808926320",
  "bbob-biobj_f24_i07_d20 0.935952364581234",
  "bbob-biobj_f24_i07_d40 0.925822002554914",
  "bbob-biobj_f24_i08_d02 0.746518610595370",
  "bbob-biobj_f24_i08_d03 0.974952781282284",
  "bbob-biobj_f24_i08_d05 0.904378583668934",
  "bbob-biobj_f24_i08_d10 0.972240400082175",
  "bbob-biobj_f24_i08_d20 0.907590679991813",
  "bbob-biobj_f24_i08_d40 0.903923697283298",
  "bbob-biobj_f24_i09_d02 0.941247141981918",
  "bbob-biobj_f24_i09_d03 0.996923371703027",
  "bbob-biobj_f24_i09_d05 0.966952435510239",
  "bbob-biobj_f24_i09_d10 0.908312375022889",
  "bbob-biobj_f24_i09_d20 0.940494648083169",
  "bbob-biobj_f24_i09_d40 0.900797679025189",
  "bbob-biobj_f24_i10_d02 0.929419348168966",
  "bbob-biobj_f24_i10_d03 0.995031721016635",
  "bbob-biobj_f24_i10_d05 0.870847461913636",
  "bbob-biobj_f24_i10_d10 0.947260657098951",
  "bbob-biobj_f24_i10_d20 0.957173792777159",
  "bbob-biobj_f24_i10_d40 0.929454939103428",
  "bbob-biobj_f25_i01_d02 0.890542925273445",
  "bbob-biobj_f25_i01_d03 0.996715655062759",
  "bbob-biobj_f25_i01_d05 0.992998562669765",
  "bbob-biobj_f25_i01_d10 0.986544179648822",
  "bbob-biobj_f25_i01_d20 0.978803085418218",
  "bbob-biobj_f25_i01_d40 0.965224595149623",
  "bbob-biobj_f25_i02_d02 0.944370794547814",
  "bbob-biobj_f25_i02_d03 0.961644138487843",
  "bbob-biobj_f25_i02_d05 0.983787947692708",
  "bbob-biobj_f25_i02_d10 0.972192877055343",
  "bbob-biobj_f25_i02_d20 0.934574214838880",
  "bbob-biobj_f25_i02_d40 0.979905509952750",
  "bbob-biobj_f25_i03_d02 0.978559983179905",
  "bbob-biobj_f25_i03_d03 0.957710702880073",
  "bbob-biobj_f25_i03_d05 0.969867210317412",
  "bbob-biobj_f25_i03_d10 0.987686732157505",
  "bbob-biobj_f25_i03_d20 0.985188585159227",
  "bbob-biobj_f25_i03_d40 0.966476888605183",
  "bbob-biobj_f25_i04_d02 0.827380317430740",
  "bbob-biobj_f25_i04_d03 0.985238959938757",
  "bbob-biobj_f25_i04_d05 0.992540046939388",
  "bbob-biobj_f25_i04_d10 0.972361968381773",
  "bbob-biobj_f25_i04_d20 0.991764320720276",
  "bbob-biobj_f25_i04_d40 0.954749790722463",
  "bbob-biobj_f25_i05_d02 0.903625072167231",
  "bbob-biobj_f25_i05_d03 0.987241523779227",
  "bbob-biobj_f25_i05_d05 0.938536624590743",
  "bbob-biobj_f25_i05_d10 0.957341098388940",
  "bbob-biobj_f25_i05_d20 0.958151869758615",
  "bbob-biobj_f25_i05_d40 0.936206497195498",
  "bbob-biobj_f25_i06_d02 0.953579508708741",
  "bbob-biobj_f25_i06_d03 0.973045616809661",
  "bbob-biobj_f25_i06_d05 0.780703874069183",
  "bbob-biobj_f25_i06_d10 0.954670927125577",
  "bbob-biobj_f25_i06_d20 0.996404383701917",
  "bbob-biobj_f25_i06_d40 0.889458931138201",
  "bbob-biobj_f25_i07_d02 0.723354491695644",
  "bbob-biobj_f25_i07_d03 0.961849853998311",
  "bbob-biobj_f25_i07_d05 0.899543698196000",
  "bbob-biobj_f25_i07_d10 0.991683054868503",
  "bbob-biobj_f25_i07_d20 0.995760710737638",
  "bbob-biobj_f25_i07_d40 0.958131259637428",
  "bbob-biobj_f25_i08_d02 0.877809309379157",
  "bbob-biobj_f25_i08_d03 0.872974526629367",
  "bbob-biobj_f25_i08_d05 0.987756096156677",
  "bbob-biobj_f25_i08_d10 0.960433783997125",
  "bbob-biobj_f25_i08_d20 0.995917374508696",
  "bbob-biobj_f25_i08_d40 0.953915373638093",
  "bbob-biobj_f25_i09_d02 0.865107177538313",
  "bbob-biobj_f25_i09_d03 0.995564468344703",
  "bbob-biobj_f25_i09_d05 0.958818899556594",
  "bbob-biobj_f25_i09_d10 0.971309227070695",
  "bbob-biobj_f25_i09_d20 0.945093366817353",
  "bbob-biobj_f25_i09_d40 0.917200882624970",
  "bbob-biobj_f25_i10_d02 0.885965221349101",
  "bbob-biobj_f25_i10_d03 0.980332318101590",
  "bbob-biobj_f25_i10_d05 0.931988275341605",
  "bbob-biobj_f25_i10_d10 0.997948127986080",
  "bbob-biobj_f25_i10_d20 0.967929115276903",
  "bbob-biobj_f25_i10_d40 0.951960079846788",
  "bbob-biobj_f26_i01_d02 0.978921935718709",
  "bbob-biobj_f26_i01_d03 0.999860002687421",
  "bbob-biobj_f26_i01_d05 0.949032783266067",
  "bbob-biobj_f26_i01_d10 0.999623700084350",
  "bbob-biobj_f26_i01_d20 0.999902242763329",
  "bbob-biobj_f26_i01_d40 0.999573599516853",
  "bbob-biobj_f26_i02_d02 0.994494596966963",
  "bbob-biobj_f26_i02_d03 0.988829166290261",
  "bbob-biobj_f26_i02_d05 0.979221374256521",
  "bbob-biobj_f26_i02_d10 0.999635083366585",
  "bbob-biobj_f26_i02_d20 0.996343737621384",
  "bbob-biobj_f26_i02_d40 0.997422279459322",
  "bbob-biobj_f26_i03_d02 0.999888275904127",
  "bbob-biobj_f26_i03_d03 0.996485238965582",
  "bbob-biobj_f26_i03_d05 0.984009970321662",
  "bbob-biobj_f26_i03_d10 0.997282897320552",
  "bbob-biobj_f26_i03_d20 0.999942519852529",
  "bbob-biobj_f26_i03_d40 0.999630688440478",
  "bbob-biobj_f26_i04_d02 0.929918978865514",
  "bbob-biobj_f26_i04_d03 0.996886144160252",
  "bbob-biobj_f26_i04_d05 0.965392151557247",
  "bbob-biobj_f26_i04_d10 0.999901970016775",
  "bbob-biobj_f26_i04_d20 0.999910304957104",
  "bbob-biobj_f26_i04_d40 0.992513325930629",
  "bbob-biobj_f26_i05_d02 0.732358950040875",
  "bbob-biobj_f26_i05_d03 0.919257062435087",
  "bbob-biobj_f26_i05_d05 0.999565003133621",
  "bbob-biobj_f26_i05_d10 0.998146402896562",
  "bbob-biobj_f26_i05_d20 0.994308174891972",
  "bbob-biobj_f26_i05_d40 0.995474176948838",
  "bbob-biobj_f26_i06_d02 0.998258046263559",
  "bbob-biobj_f26_i06_d03 0.999729023660396",
  "bbob-biobj_f26_i06_d05 0.999909507299600",
  "bbob-biobj_f26_i06_d10 0.999959216338841",
  "bbob-biobj_f26_i06_d20 0.996023601678104",
  "bbob-biobj_f26_i06_d40 0.999914106635681",
  "bbob-biobj_f26_i07_d02 0.958124544848224",
  "bbob-biobj_f26_i07_d03 0.982235104842870",
  "bbob-biobj_f26_i07_d05 0.999567163460420",
  "bbob-biobj_f26_i07_d10 0.997383627913456",
  "bbob-biobj_f26_i07_d20 0.970124626322293",
  "bbob-biobj_f26_i07_d40 0.987753037743407",
  "bbob-biobj_f26_i08_d02 0.893953976171903",
  "bbob-biobj_f26_i08_d03 0.904124381246117",
  "bbob-biobj_f26_i08_d05 0.993760975734663",
  "bbob-biobj_f26_i08_d10 0.995838547105752",
  "bbob-biobj_f26_i08_d20 0.999611932519577",
  "bbob-biobj_f26_i08_d40 0.999508639196612",
  "bbob-biobj_f26_i09_d02 0.928157845400416",
  "bbob-biobj_f26_i09_d03 0.948194726494350",
  "bbob-biobj_f26_i09_d05 0.999667008802383",
  "bbob-biobj_f26_i09_d10 0.976172358180182",
  "bbob-biobj_f26_i09_d20 0.991835909147738",
  "bbob-biobj_f26_i09_d40 0.999481548895557",
  "bbob-biobj_f26_i10_d02 0.975775907908519",
  "bbob-biobj_f26_i10_d03 0.792240379786697",
  "bbob-biobj_f26_i10_d05 0.992640176232114",
  "bbob-biobj_f26_i10_d10 0.999411169788605",
  "bbob-biobj_f26_i10_d20 0.986066226299899",
  "bbob-biobj_f26_i10_d40 0.999836115712070",
  "bbob-biobj_f27_i01_d02 0.903502206401833",
  "bbob-biobj_f27_i01_d03 0.987126443717575",
  "bbob-biobj_f27_i01_d05 0.992203655045538",
  "bbob-biobj_f27_i01_d10 0.981631550497566",
  "bbob-biobj_f27_i01_d20 0.977451506783242",
  "bbob-biobj_f27_i01_d40 0.929804876788547",
  "bbob-biobj_f27_i02_d02 0.951958966128448",
  "bbob-biobj_f27_i02_d03 0.953825810996211",
  "bbob-biobj_f27_i02_d05 0.932656649627805",
  "bbob-biobj_f27_i02_d10 0.903587534644867",
  "bbob-biobj_f27_i02_d20 0.966380883152270",
  "bbob-biobj_f27_i02_d40 0.890734610987759",
  "bbob-biobj_f27_i03_d02 0.957759312399377",
  "bbob-biobj_f27_i03_d03 0.964059250114641",
  "bbob-biobj_f27_i03_d05 0.978490126422125",
  "bbob-biobj_f27_i03_d10 0.991361907611424",
  "bbob-biobj_f27_i03_d20 0.966345304052504",
  "bbob-biobj_f27_i03_d40 0.920565628885986",
  "bbob-biobj_f27_i04_d02 0.960857690938292",
  "bbob-biobj_f27_i04_d03 0.952867272760359",
  "bbob-biobj_f27_i04_d05 0.981448840032206",
  "bbob-biobj_f27_i04_d10 0.987121596080223",
  "bbob-biobj_f27_i04_d20 0.978983096941641",
  "bbob-biobj_f27_i04_d40 0.876799658083650",
  "bbob-biobj_f27_i05_d02 0.930274710740560",
  "bbob-biobj_f27_i05_d03 0.960952531243144",
  "bbob-biobj_f27_i05_d05 0.975425675306981",
  "bbob-biobj_f27_i05_d10 0.981695367312035",
  "bbob-biobj_f27_i05_d20 0.926111704541228",
  "bbob-biobj_f27_i05_d40 0.841173541405940",
  "bbob-biobj_f27_i06_d02 0.962259536027885",
  "bbob-biobj_f27_i06_d03 0.993495745010671",
  "bbob-biobj_f27_i06_d05 0.964350154243847",
  "bbob-biobj_f27_i06_d10 0.952265863863045",
  "bbob-biobj_f27_i06_d20 0.949060318808586",
  "bbob-biobj_f27_i06_d40 0.881168390379038",
  "bbob-biobj_f27_i07_d02 0.950861129749424",
  "bbob-biobj_f27_i07_d03 0.946795991866674",
  "bbob-biobj_f27_i07_d05 0.997084717538608",
  "bbob-biobj_f27_i07_d10 0.964754431608966",
  "bbob-biobj_f27_i07_d20 0.931452740418371",
  "bbob-biobj_f27_i07_d40 0.800872576192903",
  "bbob-biobj_f27_i08_d02 0.989189557992834",
  "bbob-biobj_f27_i08_d03 0.944503236391228",
  "bbob-biobj_f27_i08_d05 0.981904737031804",
  "bbob-biobj_f27_i08_d10 0.995675312704574",
  "bbob-biobj_f27_i08_d20 0.940148066198553",
  "bbob-biobj_f27_i08_d40 0.760209371588989",
  "bbob-biobj_f27_i09_d02 0.931408553538173",
  "bbob-biobj_f27_i09_d03 0.986109779620887",
  "bbob-biobj_f27_i09_d05 0.999682888126238",
  "bbob-biobj_f27_i09_d10 0.925545214667907",
  "bbob-biobj_f27_i09_d20 0.983176894276770",
  "bbob-biobj_f27_i09_d40 0.853713505733515",
  "bbob-biobj_f27_i10_d02 0.976300396957638",
  "bbob-biobj_f27_i10_d03 0.995369096891944",
  "bbob-biobj_f27_i10_d05 0.937129382454360",
  "bbob-biobj_f27_i10_d10 0.955167860062906",
  "bbob-biobj_f27_i10_d20 0.980491887614351",
  "bbob-biobj_f27_i10_d40 0.787140818600140",
  "bbob-biobj_f28_i01_d02 0.977401739849146",
  "bbob-biobj_f28_i01_d03 0.998639599209397",
  "bbob-biobj_f28_i01_d05 0.995557284150280",
  "bbob-biobj_f28_i01_d10 0.994068867829939",
  "bbob-biobj_f28_i01_d20 0.992454480815253",
  "bbob-biobj_f28_i01_d40 0.992341193444900",
  "bbob-biobj_f28_i02_d02 0.998930188968113",
  "bbob-biobj_f28_i02_d03 0.993519352925855",
  "bbob-biobj_f28_i02_d05 0.991805294302727",
  "bbob-biobj_f28_i02_d10 0.992228309135279",
  "bbob-biobj_f28_i02_d20 0.990006492841615",
  "bbob-biobj_f28_i02_d40 0.990087803846637",
  "bbob-biobj_f28_i03_d02 0.999666277505892",
  "bbob-biobj_f28_i03_d03 0.977206371838020",
  "bbob-biobj_f28_i03_d05 0.993681814867215",
  "bbob-biobj_f28_i03_d10 0.994078408234904",
  "bbob-biobj_f28_i03_d20 0.993728811755353",
  "bbob-biobj_f28_i03_d40 0.993502900831588",
  "bbob-biobj_f28_i04_d02 0.984721338670538",
  "bbob-biobj_f28_i04_d03 0.992073902063011",
  "bbob-biobj_f28_i04_d05 0.997392752918237",
  "bbob-biobj_f28_i04_d10 0.992800925135662",
  "bbob-biobj_f28_i04_d20 0.994516641552363",
  "bbob-biobj_f28_i04_d40 0.992158708840145",
  "bbob-biobj_f28_i05_d02 0.999901255201244",
  "bbob-biobj_f28_i05_d03 0.990126928654807",
  "bbob-biobj_f28_i05_d05 0.987356824339541",
  "bbob-biobj_f28_i05_d10 0.993737253443776",
  "bbob-biobj_f28_i05_d20 0.995292057751840",
  "bbob-biobj_f28_i05_d40 0.992463241466847",
  "bbob-biobj_f28_i06_d02 0.999625111303881",
  "bbob-biobj_f28_i06_d03 0.958055936488174",
  "bbob-biobj_f28_i06_d05 0.997861945856941",
  "bbob-biobj_f28_i06_d10 0.993863496568245",
  "bbob-biobj_f28_i06_d20 0.992376930928585",
  "bbob-biobj_f28_i06_d40 0.990946665473032",
  "bbob-biobj_f28_i07_d02 0.997311511998985",
  "bbob-biobj_f28_i07_d03 0.991241887296098",
  "bbob-biobj_f28_i07_d05 0.991900279023225",
  "bbob-biobj_f28_i07_d10 0.995102876975957",
  "bbob-biobj_f28_i07_d20 0.991044125475600",
  "bbob-biobj_f28_i07_d40 0.991171901135509",
  "bbob-biobj_f28_i08_d02 0.983531017082843",
  "bbob-biobj_f28_i08_d03 0.969999371854049",
  "bbob-biobj_f28_i08_d05 0.991730556303096",
  "bbob-biobj_f28_i08_d10 0.994108454851434",
  "bbob-biobj_f28_i08_d20 0.990408448797344",
  "bbob-biobj_f28_i08_d40 0.991668928043641",
  "bbob-biobj_f28_i09_d02 0.992137066959289",
  "bbob-biobj_f28_i09_d03 0.978690300285065",
  "bbob-biobj_f28_i09_d05 0.995078829300506",
  "bbob-biobj_f28_i09_d10 0.995257886909117",
  "bbob-biobj_f28_i09_d20 0.992380716938039",
  "bbob-biobj_f28_i09_d40 0.992595150586751",
  "bbob-biobj_f28_i10_d02 0.999840446818862",
  "bbob-biobj_f28_i10_d03 0.994813201637420",
  "bbob-biobj_f28_i10_d05 0.989124585591138",
  "bbob-biobj_f28_i10_d10 0.996349741278776",
  "bbob-biobj_f28_i10_d20 0.993618129845296",
  "bbob-biobj_f28_i10_d40 0.989833089219968",
  "bbob-biobj_f29_i01_d02 0.972456391332674",
  "bbob-biobj_f29_i01_d03 0.866365733646573",
  "bbob-biobj_f29_i01_d05 0.870418605165827",
  "bbob-biobj_f29_i01_d10 0.894362096446415",
  "bbob-biobj_f29_i01_d20 0.808745063728655",
  "bbob-biobj_f29_i01_d40 0.840413324584372",
  "bbob-biobj_f29_i02_d02 0.999165519812747",
  "bbob-biobj_f29_i02_d03 0.939190829078205",
  "bbob-biobj_f29_i02_d05 0.882799523408890",
  "bbob-biobj_f29_i02_d10 0.804946317347259",
  "bbob-biobj_f29_i02_d20 0.859875275398261",
  "bbob-biobj_f29_i02_d40 0.835680010264772",
  "bbob-biobj_f29_i03_d02 0.993280309403803",
  "bbob-biobj_f29_i03_d03 0.980540524729557",
  "bbob-biobj_f29_i03_d05 0.830291271062953",
  "bbob-biobj_f29_i03_d10 0.851081962925671",
  "bbob-biobj_f29_i03_d20 0.830665167637888",
  "bbob-biobj_f29_i03_d40 0.837194092497397",
  "bbob-biobj_f29_i04_d02 0.966908599055087",
  "bbob-biobj_f29_i04_d03 0.973838150777401",
  "bbob-biobj_f29_i04_d05 0.899714748221521",
  "bbob-biobj_f29_i04_d10 0.894919059851961",
  "bbob-biobj_f29_i04_d20 0.865257050397847",
  "bbob-biobj_f29_i04_d40 0.836624327182609",
  "bbob-biobj_f29_i05_d02 0.988529005467012",
  "bbob-biobj_f29_i05_d03 0.949154431613501",
  "bbob-biobj_f29_i05_d05 0.927435967382696",
  "bbob-biobj_f29_i05_d10 0.900608993514346",
  "bbob-biobj_f29_i05_d20 0.835978832033558",
  "bbob-biobj_f29_i05_d40 0.845940916781953",
  "bbob-biobj_f29_i06_d02 0.967527569931253",
  "bbob-biobj_f29_i06_d03 0.967272857789836",
  "bbob-biobj_f29_i06_d05 0.930681778227153",
  "bbob-biobj_f29_i06_d10 0.847546798196211",
  "bbob-biobj_f29_i06_d20 0.842302730724967",
  "bbob-biobj_f29_i06_d40 0.823597979477189",
  "bbob-biobj_f29_i07_d02 0.981106173671041",
  "bbob-biobj_f29_i07_d03 0.957195556159746",
  "bbob-biobj_f29_i07_d05 0.902382575069356",
  "bbob-biobj_f29_i07_d10 0.946969100732708",
  "bbob-biobj_f29_i07_d20 0.821540236470965",
  "bbob-biobj_f29_i07_d40 0.860947536331094",
  "bbob-biobj_f29_i08_d02 0.984725351402661",
  "bbob-biobj_f29_i08_d03 0.924479789647542",
  "bbob-biobj_f29_i08_d05 0.804202269471248",
  "bbob-biobj_f29_i08_d10 0.918960642476287",
  "bbob-biobj_f29_i08_d20 0.842064684675849",
  "bbob-biobj_f29_i08_d40 0.851779889000568",
  "bbob-biobj_f29_i09_d02 0.992958120825111",
  "bbob-biobj_f29_i09_d03 0.997309840895446",
  "bbob-biobj_f29_i09_d05 0.872251456961736",
  "bbob-biobj_f29_i09_d10 0.868432629364904",
  "bbob-biobj_f29_i09_d20 0.832038333787257",
  "bbob-biobj_f29_i09_d40 0.823950193567713",
  "bbob-biobj_f29_i10_d02 0.998395344316965",
  "bbob-biobj_f29_i10_d03 0.957698564955352",
  "bbob-biobj_f29_i10_d05 0.987513385842123",
  "bbob-biobj_f29_i10_d10 0.855735639108079",
  "bbob-biobj_f29_i10_d20 0.862165546155237",
  "bbob-biobj_f29_i10_d40 0.814473985769004",
  "bbob-biobj_f30_i01_d02 0.976212391884318",
  "bbob-biobj_f30_i01_d03 0.795597399545723",
  "bbob-biobj_f30_i01_d05 0.940209615676908",
  "bbob-biobj_f30_i01_d10 0.926781574796662",
  "bbob-biobj_f30_i01_d20 0.962010692649429",
  "bbob-biobj_f30_i01_d40 0.971315600273112",
  "bbob-biobj_f30_i02_d02 0.996222766326675",
  "bbob-biobj_f30_i02_d03 0.931860253086072",
  "bbob-biobj_f30_i02_d05 0.982391031372733",
  "bbob-biobj_f30_i02_d10 0.979266611427294",
  "bbob-biobj_f30_i02_d20 0.948073614979061",
  "bbob-biobj_f30_i02_d40 0.968892101028860",
  "bbob-biobj_f30_i03_d02 0.941596026692009",
  "bbob-biobj_f30_i03_d03 0.913952880246775",
  "bbob-biobj_f30_i03_d05 0.914681224527996",
  "bbob-biobj_f30_i03_d10 0.986318664324346",
  "bbob-biobj_f30_i03_d20 0.968008843765655",
  "bbob-biobj_f30_i03_d40 0.963738386811532",
  "bbob-biobj_f30_i04_d02 0.983029998728216",
  "bbob-biobj_f30_i04_d03 0.977327690333281",
  "bbob-biobj_f30_i04_d05 0.979631549217623",
  "bbob-biobj_f30_i04_d10 0.980901390335724",
  "bbob-biobj_f30_i04_d20 0.980970593729196",
  "bbob-biobj_f30_i04_d40 0.989217314263775",
  "bbob-biobj_f30_i05_d02 0.995608409377838",
  "bbob-biobj_f30_i05_d03 0.938513468728340",
  "bbob-biobj_f30_i05_d05 0.933170411362319",
  "bbob-biobj_f30_i05_d10 0.978390758883721",
  "bbob-biobj_f30_i05_d20 0.986547939416421",
  "bbob-biobj_f30_i05_d40 0.943152144087444",
  "bbob-biobj_f30_i06_d02 0.868995465253104",
  "bbob-biobj_f30_i06_d03 0.941098047213624",
  "bbob-biobj_f30_i06_d05 0.956637817519768",
  "bbob-biobj_f30_i06_d10 0.943048982247906",
  "bbob-biobj_f30_i06_d20 0.983026375287029",
  "bbob-biobj_f30_i06_d40 0.981410556945632",
  "bbob-biobj_f30_i07_d02 0.850916039694341",
  "bbob-biobj_f30_i07_d03 0.992456563611573",
  "bbob-biobj_f30_i07_d05 0.978264480487193",
  "bbob-biobj_f30_i07_d10 0.961489543939775",
  "bbob-biobj_f30_i07_d20 0.977377727125355",
  "bbob-biobj_f30_i07_d40 0.983967640342075",
  "bbob-biobj_f30_i08_d02 0.994797413889769",
  "bbob-biobj_f30_i08_d03 0.971819487123594",
  "bbob-biobj_f30_i08_d05 0.983761767209514",
  "bbob-biobj_f30_i08_d10 0.979219261565601",
  "bbob-biobj_f30_i08_d20 0.962465917325102",
  "bbob-biobj_f30_i08_d40 0.989627476762782",
  "bbob-biobj_f30_i09_d02 0.940336082398701",
  "bbob-biobj_f30_i09_d03 0.992719488303988",
  "bbob-biobj_f30_i09_d05 0.985260179067806",
  "bbob-biobj_f30_i09_d10 0.951799238882504",
  "bbob-biobj_f30_i09_d20 0.958847139142673",
  "bbob-biobj_f30_i09_d40 0.980832045200368",
  "bbob-biobj_f30_i10_d02 0.832977994430622",
  "bbob-biobj_f30_i10_d03 0.970997818435150",
  "bbob-biobj_f30_i10_d05 0.987861046004904",
  "bbob-biobj_f30_i10_d10 0.964527185265092",
  "bbob-biobj_f30_i10_d20 0.969249792773688",
  "bbob-biobj_f30_i10_d40 0.955879199143995",
  "bbob-biobj_f31_i01_d02 0.955984085871603",
  "bbob-biobj_f31_i01_d03 0.988038414994945",
  "bbob-biobj_f31_i01_d05 0.973301456584575",
  "bbob-biobj_f31_i01_d10 0.960797584654700",
  "bbob-biobj_f31_i01_d20 0.964632857214670",
  "bbob-biobj_f31_i01_d40 0.964644929471550",
  "bbob-biobj_f31_i02_d02 0.963697231805374",
  "bbob-biobj_f31_i02_d03 0.923773287001374",
  "bbob-biobj_f31_i02_d05 0.977069222171345",
  "bbob-biobj_f31_i02_d10 0.943110573823707",
  "bbob-biobj_f31_i02_d20 0.966096830199769",
  "bbob-biobj_f31_i02_d40 0.928921378338131",
  "bbob-biobj_f31_i03_d02 0.981722020953158",
  "bbob-biobj_f31_i03_d03 0.984775969014314",
  "bbob-biobj_f31_i03_d05 0.975988405888085",
  "bbob-biobj_f31_i03_d10 0.954640660032854",
  "bbob-biobj_f31_i03_d20 0.970283232049837",
  "bbob-biobj_f31_i03_d40 0.951982517581485",
  "bbob-biobj_f31_i04_d02 0.952958511212198",
  "bbob-biobj_f31_i04_d03 0.994016367179768",
  "bbob-biobj_f31_i04_d05 0.977972687898182",
  "bbob-biobj_f31_i04_d10 0.969513406975955",
  "bbob-biobj_f31_i04_d20 0.974626124154705",
  "bbob-biobj_f31_i04_d40 0.952265968520882",
  "bbob-biobj_f31_i05_d02 0.989988682239425",
  "bbob-biobj_f31_i05_d03 0.975602053920554",
  "bbob-biobj_f31_i05_d05 0.985500067314789",
  "bbob-biobj_f31_i05_d10 0.973373101369873",
  "bbob-biobj_f31_i05_d20 0.960646486040987",
  "bbob-biobj_f31_i05_d40 0.954452695883527",
  "bbob-biobj_f31_i06_d02 0.993336361299947",
  "bbob-biobj_f31_i06_d03 0.997761535471247",
  "bbob-biobj_f31_i06_d05 0.967003360813607",
  "bbob-biobj_f31_i06_d10 0.964797907545882",
  "bbob-biobj_f31_i06_d20 0.965798879303843",
  "bbob-biobj_f31_i06_d40 0.977483934244301",
  "bbob-biobj_f31_i07_d02 0.974719524344369",
  "bbob-biobj_f31_i07_d03 0.980459056469975",
  "bbob-biobj_f31_i07_d05 0.962328735477573",
  "bbob-biobj_f31_i07_d10 0.966870060432519",
  "bbob-biobj_f31_i07_d20 0.964033021159753",
  "bbob-biobj_f31_i07_d40 0.883371649386440",
  "bbob-biobj_f31_i08_d02 0.979768767126209",
  "bbob-biobj_f31_i08_d03 0.941759522629115",
  "bbob-biobj_f31_i08_d05 0.957011027742155",
  "bbob-biobj_f31_i08_d10 0.957146360706000",
  "bbob-biobj_f31_i08_d20 0.965744801674200",
  "bbob-biobj_f31_i08_d40 0.924221235535711",
  "bbob-biobj_f31_i09_d02 0.975491331974885",
  "bbob-biobj_f31_i09_d03 0.945701056779519",
  "bbob-biobj_f31_i09_d05 0.948252164823377",
  "bbob-biobj_f31_i09_d10 0.985685779893663",
  "bbob-biobj_f31_i09_d20 0.961014373658260",
  "bbob-biobj_f31_i09_d40 0.930629055574162",
  "bbob-biobj_f31_i10_d02 0.933792954618640",
  "bbob-biobj_f31_i10_d03 0.961379021599416",
  "bbob-biobj_f31_i10_d05 0.941053258377657",
  "bbob-biobj_f31_i10_d10 0.972487214421872",
  "bbob-biobj_f31_i10_d20 0.970440969764005",
  "bbob-biobj_f31_i10_d40 0.932460885786384",
  "bbob-biobj_f32_i01_d02 0.920330761005554",
  "bbob-biobj_f32_i01_d03 0.915330782748840",
  "bbob-biobj_f32_i01_d05 0.972599849011612",
  "bbob-biobj_f32_i01_d10 0.949204788073746",
  "bbob-biobj_f32_i01_d20 0.981515702413721",
  "bbob-biobj_f32_i01_d40 0.961292288128107",
  "bbob-biobj_f32_i02_d02 0.675234316067362",
  "bbob-biobj_f32_i02_d03 0.922150714718287",
  "bbob-biobj_f32_i02_d05 0.938528527832380",
  "bbob-biobj_f32_i02_d10 0.964298605845870",
  "bbob-biobj_f32_i02_d20 0.960802272582804",
  "bbob-biobj_f32_i02_d40 0.957508035210367",
  "bbob-biobj_f32_i03_d02 0.921767715599282",
  "bbob-biobj_f32_i03_d03 0.968198174673630",
  "bbob-biobj_f32_i03_d05 0.983105781597808",
  "bbob-biobj_f32_i03_d10 0.962342312174357",
  "bbob-biobj_f32_i03_d20 0.979625645472811",
  "bbob-biobj_f32_i03_d40 0.965004188074738",
  "bbob-biobj_f32_i04_d02 0.944003925260494",
  "bbob-biobj_f32_i04_d03 0.906090395552451",
  "bbob-biobj_f32_i04_d05 0.987970425641090",
  "bbob-biobj_f32_i04_d10 0.987762782583046",
  "bbob-biobj_f32_i04_d20 0.989911105052668",
  "bbob-biobj_f32_i04_d40 0.962535091946121",
  "bbob-biobj_f32_i05_d02 0.844964277490244",
  "bbob-biobj_f32_i05_d03 0.983676044504197",
  "bbob-biobj_f32_i05_d05 0.982593769034902",
  "bbob-biobj_f32_i05_d10 0.993085374531154",
  "bbob-biobj_f32_i05_d20 0.989167495743449",
  "bbob-biobj_f32_i05_d40 0.951371297338482",
  "bbob-biobj_f32_i06_d02 0.936185223968449",
  "bbob-biobj_f32_i06_d03 0.955270985849230",
  "bbob-biobj_f32_i06_d05 0.940138737948394",
  "bbob-biobj_f32_i06_d10 0.940548806079393",
  "bbob-biobj_f32_i06_d20 0.988426823934137",
  "bbob-biobj_f32_i06_d40 0.911350104998297",
  "bbob-biobj_f32_i07_d02 0.957426888523200",
  "bbob-biobj_f32_i07_d03 0.941139628762055",
  "bbob-biobj_f32_i07_d05 0.975765546822080",
  "bbob-biobj_f32_i07_d10 0.977883563985928",
  "bbob-biobj_f32_i07_d20 0.977722436173549",
  "bbob-biobj_f32_i07_d40 0.889732441115030",
  "bbob-biobj_f32_i08_d02 0.969116873902084",
  "bbob-biobj_f32_i08_d03 0.995564101875544",
  "bbob-biobj_f32_i08_d05 0.985946793275347",
  "bbob-biobj_f32_i08_d10 0.967235711724125",
  "bbob-biobj_f32_i08_d20 0.963585562806338",
  "bbob-biobj_f32_i08_d40 0.947796848503868",
  "bbob-biobj_f32_i09_d02 0.721545258466030",
  "bbob-biobj_f32_i09_d03 0.986354192081451",
  "bbob-biobj_f32_i09_d05 0.961851852651161",
  "bbob-biobj_f32_i09_d10 0.964787664322737",
  "bbob-biobj_f32_i09_d20 0.944657929223530",
  "bbob-biobj_f32_i09_d40 0.940779793758805",
  "bbob-biobj_f32_i10_d02 0.650792309240297",
  "bbob-biobj_f32_i10_d03 0.975518606982846",
  "bbob-biobj_f32_i10_d05 0.955570275426794",
  "bbob-biobj_f32_i10_d10 0.987721955466459",
  "bbob-biobj_f32_i10_d20 0.983584157948028",
  "bbob-biobj_f32_i10_d40 0.918763000553540",
  "bbob-biobj_f33_i01_d02 0.997888833151880",
  "bbob-biobj_f33_i01_d03 0.997624564773393",
  "bbob-biobj_f33_i01_d05 0.998341905642183",
  "bbob-biobj_f33_i01_d10 0.995514031710703",
  "bbob-biobj_f33_i01_d20 0.990225223031912",
  "bbob-biobj_f33_i01_d40 0.996609751740017",
  "bbob-biobj_f33_i02_d02 0.997492967539081",
  "bbob-biobj_f33_i02_d03 0.998822646695541",
  "bbob-biobj_f33_i02_d05 0.989647422804803",
  "bbob-biobj_f33_i02_d10 0.995930212809998",
  "bbob-biobj_f33_i02_d20 0.995149035345042",
  "bbob-biobj_f33_i02_d40 0.993996885864406",
  "bbob-biobj_f33_i03_d02 0.999823325123336",
  "bbob-biobj_f33_i03_d03 0.935833285400742",
  "bbob-biobj_f33_i03_d05 0.999018147217449",
  "bbob-biobj_f33_i03_d10 0.998455572034838",
  "bbob-biobj_f33_i03_d20 0.998360083847699",
  "bbob-biobj_f33_i03_d40 0.994133731269471",
  "bbob-biobj_f33_i04_d02 0.959892839230054",
  "bbob-biobj_f33_i04_d03 0.989064148620451",
  "bbob-biobj_f33_i04_d05 0.999673770334207",
  "bbob-biobj_f33_i04_d10 0.996076805375278",
  "bbob-biobj_f33_i04_d20 0.996910256375666",
  "bbob-biobj_f33_i04_d40 0.994774774892129",
  "bbob-biobj_f33_i05_d02 0.966110377602840",
  "bbob-biobj_f33_i05_d03 0.998861219424993",
  "bbob-biobj_f33_i05_d05 0.999856944246115",
  "bbob-biobj_f33_i05_d10 0.999704223934855",
  "bbob-biobj_f33_i05_d20 0.994309244425620",
  "bbob-biobj_f33_i05_d40 0.994350807712746",
  "bbob-biobj_f33_i06_d02 0.988007874945376",
  "bbob-biobj_f33_i06_d03 0.998475157629897",
  "bbob-biobj_f33_i06_d05 0.994225205472486",
  "bbob-biobj_f33_i06_d10 0.997405379840233",
  "bbob-biobj_f33_i06_d20 0.995656376599321",
  "bbob-biobj_f33_i06_d40 0.991034860399352",
  "bbob-biobj_f33_i07_d02 0.908411632033230",
  "bbob-biobj_f33_i07_d03 0.893632373750011",
  "bbob-biobj_f33_i07_d05 0.975756013845648",
  "bbob-biobj_f33_i07_d10 0.999358629830707",
  "bbob-biobj_f33_i07_d20 0.986490502693860",
  "bbob-biobj_f33_i07_d40 0.994517293149754",
  "bbob-biobj_f33_i08_d02 0.999898764570977",
  "bbob-biobj_f33_i08_d03 0.995361006858616",
  "bbob-biobj_f33_i08_d05 0.999234525330710",
  "bbob-biobj_f33_i08_d10 0.997546804884822",
  "bbob-biobj_f33_i08_d20 0.998045834082956",
  "bbob-biobj_f33_i08_d40 0.995593060825240",
  "bbob-biobj_f33_i09_d02 0.797548074305266",
  "bbob-biobj_f33_i09_d03 0.817352004149393",
  "bbob-biobj_f33_i09_d05 0.869142706422598",
  "bbob-biobj_f33_i09_d10 0.997371335975271",
  "bbob-biobj_f33_i09_d20 0.995961945805410",
  "bbob-biobj_f33_i09_d40 0.991305093363462",
  "bbob-biobj_f33_i10_d02 0.999797872898550",
  "bbob-biobj_f33_i10_d03 0.999774146966997",
  "bbob-biobj_f33_i10_d05 0.998634770875355",
  "bbob-biobj_f33_i10_d10 0.999968440732627",
  "bbob-biobj_f33_i10_d20 0.994426545919763",
  "bbob-biobj_f33_i10_d40 0.990815922897162",
  "bbob-biobj_f34_i01_d02 0.929758805624314",
  "bbob-biobj_f34_i01_d03 0.914923801674542",
  "bbob-biobj_f34_i01_d05 0.981336494435819",
  "bbob-biobj_f34_i01_d10 0.966127872335118",
  "bbob-biobj_f34_i01_d20 0.950007871264013",
  "bbob-biobj_f34_i01_d40 0.936082386113844",
  "bbob-biobj_f34_i02_d02 0.996825921926922",
  "bbob-biobj_f34_i02_d03 0.914729610596124",
  "bbob-biobj_f34_i02_d05 0.991099691689551",
  "bbob-biobj_f34_i02_d10 0.959851052333613",
  "bbob-biobj_f34_i02_d20 0.945759794539885",
  "bbob-biobj_f34_i02_d40 0.927588070673910",
  "bbob-biobj_f34_i03_d02 0.990964371830805",
  "bbob-biobj_f34_i03_d03 0.962893751153050",
  "bbob-biobj_f34_i03_d05 0.968962646178077",
  "bbob-biobj_f34_i03_d10 0.983103032306111",
  "bbob-biobj_f34_i03_d20 0.940493218694200",
  "bbob-biobj_f34_i03_d40 0.906851732579638",
  "bbob-biobj_f34_i04_d02 0.967251356822567",
  "bbob-biobj_f34_i04_d03 0.985996586463376",
  "bbob-biobj_f34_i04_d05 0.990433096396692",
  "bbob-biobj_f34_i04_d10 0.985660370601553",
  "bbob-biobj_f34_i04_d20 0.942740668741104",
  "bbob-biobj_f34_i04_d40 0.946949316285049",
  "bbob-biobj_f34_i05_d02 0.956107364487255",
  "bbob-biobj_f34_i05_d03 0.990692065957800",
  "bbob-biobj_f34_i05_d05 0.991915439495093",
  "bbob-biobj_f34_i05_d10 0.986747682162292",
  "bbob-biobj_f34_i05_d20 0.946925944116485",
  "bbob-biobj_f34_i05_d40 0.918420250676018",
  "bbob-biobj_f34_i06_d02 0.992423068217570",
  "bbob-biobj_f34_i06_d03 0.928977413475456",
  "bbob-biobj_f34_i06_d05 0.977722500230699",
  "bbob-biobj_f34_i06_d10 0.958552499087304",
  "bbob-biobj_f34_i06_d20 0.893572201285847",
  "bbob-biobj_f34_i06_d40 0.931391786707158",
  "bbob-biobj_f34_i07_d02 0.997952817875845",
  "bbob-biobj_f34_i07_d03 0.968799951050218",
  "bbob-biobj_f34_i07_d05 0.968314701581942",
  "bbob-biobj_f34_i07_d10 0.987711893847965",
  "bbob-biobj_f34_i07_d20 0.967844420222837",
  "bbob-biobj_f34_i07_d40 0.917119551354191",
  "bbob-biobj_f34_i08_d02 0.986421616836380",
  "bbob-biobj_f34_i08_d03 0.987627767463273",
  "bbob-biobj_f34_i08_d05 0.973271016818824",
  "bbob-biobj_f34_i08_d10 0.978710309708377",
  "bbob-biobj_f34_i08_d20 0.953784614896843",
  "bbob-biobj_f34_i08_d40 0.924545814131870",
  "bbob-biobj_f34_i09_d02 0.979369454189214",
  "bbob-biobj_f34_i09_d03 0.985042366556295",
  "bbob-biobj_f34_i09_d05 0.971630663710999",
  "bbob-biobj_f34_i09_d10 0.974264518135895",
  "bbob-biobj_f34_i09_d20 0.953888014994909",
  "bbob-biobj_f34_i09_d40 0.898550158585044",
  "bbob-biobj_f34_i10_d02 0.986274175042359",
  "bbob-biobj_f34_i10_d03 0.995494022043420",
  "bbob-biobj_f34_i10_d05 0.979280689879259",
  "bbob-biobj_f34_i10_d10 0.966869492505104",
  "bbob-biobj_f34_i10_d20 0.977561648845395",
  "bbob-biobj_f34_i10_d40 0.952637742173473",
  "bbob-biobj_f35_i01_d02 0.928163128954411",
  "bbob-biobj_f35_i01_d03 0.534050929347080",
  "bbob-biobj_f35_i01_d05 0.635757565911625",
  "bbob-biobj_f35_i01_d10 0.580861243059274",
  "bbob-biobj_f35_i01_d20 0.554977949908380",
  "bbob-biobj_f35_i01_d40 0.559204648515045",
  "bbob-biobj_f35_i02_d02 0.502441647898970",
  "bbob-biobj_f35_i02_d03 0.788703796796430",
  "bbob-biobj_f35_i02_d05 0.581810861312113",
  "bbob-biobj_f35_i02_d10 0.592285160032464",
  "bbob-biobj_f35_i02_d20 0.562915958081632",
  "bbob-biobj_f35_i02_d40 0.592165291851369",
  "bbob-biobj_f35_i03_d02 0.875495049923549",
  "bbob-biobj_f35_i03_d03 0.573978025502487",
  "bbob-biobj_f35_i03_d05 0.768235426152865",
  "bbob-biobj_f35_i03_d10 0.606821886661125",
  "bbob-biobj_f35_i03_d20 0.608130174377868",
  "bbob-biobj_f35_i03_d40 0.570819141303950",
  "bbob-biobj_f35_i04_d02 0.818249550895629",
  "bbob-biobj_f35_i04_d03 0.945462936313717",
  "bbob-biobj_f35_i04_d05 0.735427134771649",
  "bbob-biobj_f35_i04_d10 0.582583984930894",
  "bbob-biobj_f35_i04_d20 0.561790409553075",
  "bbob-biobj_f35_i04_d40 0.578303629613508",
  "bbob-biobj_f35_i05_d02 0.572775304213331",
  "bbob-biobj_f35_i05_d03 0.773374133410163",
  "bbob-biobj_f35_i05_d05 0.761506440836008",
  "bbob-biobj_f35_i05_d10 0.597825128770010",
  "bbob-biobj_f35_i05_d20 0.540833243230470",
  "bbob-biobj_f35_i05_d40 0.564885532467370",
  "bbob-biobj_f35_i06_d02 0.980526794424884",
  "bbob-biobj_f35_i06_d03 0.835093857831682",
  "bbob-biobj_f35_i06_d05 0.562422717969856",
  "bbob-biobj_f35_i06_d10 0.587956869935200",
  "bbob-biobj_f35_i06_d20 0.610119471720179",
  "bbob-biobj_f35_i06_d40 0.557365525031513",
  "bbob-biobj_f35_i07_d02 0.990094491164167",
  "bbob-biobj_f35_i07_d03 0.582738983892805",
  "bbob-biobj_f35_i07_d05 0.622466269684858",
  "bbob-biobj_f35_i07_d10 0.558632303567902",
  "bbob-biobj_f35_i07_d20 0.581873550670568",
  "bbob-biobj_f35_i07_d40 0.575436782021753",
  "bbob-biobj_f35_i08_d02 0.585023413596274",
  "bbob-biobj_f35_i08_d03 0.650042559787474",
  "bbob-biobj_f35_i08_d05 0.526989002478125",
  "bbob-biobj_f35_i08_d10 0.589577017551821",
  "bbob-biobj_f35_i08_d20 0.560917276677095",
  "bbob-biobj_f35_i08_d40 0.571605429690267",
  "bbob-biobj_f35_i09_d02 0.811980149928068",
  "bbob-biobj_f35_i09_d03 0.853207710551811",
  "bbob-biobj_f35_i09_d05 0.536439303798306",
  "bbob-biobj_f35_i09_d10 0.540687322788878",
  "bbob-biobj_f35_i09_d20 0.596178199415022",
  "bbob-biobj_f35_i09_d40 0.557266787704120",
  "bbob-biobj_f35_i10_d02 0.997530732472314",
  "bbob-biobj_f35_i10_d03 0.587399626683934",
  "bbob-biobj_f35_i10_d05 0.714025877194927",
  "bbob-biobj_f35_i10_d10 0.588731628553982",
  "bbob-biobj_f35_i10_d20 0.602156913481219",
  "bbob-biobj_f35_i10_d40 0.578694053343176",
  "bbob-biobj_f36_i01_d02 0.945913658069173",
  "bbob-biobj_f36_i01_d03 0.700105030960641",
  "bbob-biobj_f36_i01_d05 0.885933163539674",
  "bbob-biobj_f36_i01_d10 0.781245213715888",
  "bbob-biobj_f36_i01_d20 0.804649167961649",
  "bbob-biobj_f36_i01_d40 0.843698667633701",
  "bbob-biobj_f36_i02_d02 0.982338880101769",
  "bbob-biobj_f36_i02_d03 0.737150249400376",
  "bbob-biobj_f36_i02_d05 0.930156906176818",
  "bbob-biobj_f36_i02_d10 0.793762012712450",
  "bbob-biobj_f36_i02_d20 0.853214386323061",
  "bbob-biobj_f36_i02_d40 0.813601612755896",
  "bbob-biobj_f36_i03_d02 0.808774401854655",
  "bbob-biobj_f36_i03_d03 0.842941654847907",
  "bbob-biobj_f36_i03_d05 0.627301358873446",
  "bbob-biobj_f36_i03_d10 0.756205125594118",
  "bbob-biobj_f36_i03_d20 0.787411418855015",
  "bbob-biobj_f36_i03_d40 0.828009661286235",
  "bbob-biobj_f36_i04_d02 0.913201943633561",
  "bbob-biobj_f36_i04_d03 0.970230637330057",
  "bbob-biobj_f36_i04_d05 0.914707866768181",
  "bbob-biobj_f36_i04_d10 0.850564650407951",
  "bbob-biobj_f36_i04_d20 0.765292066767800",
  "bbob-biobj_f36_i04_d40 0.832459999962591",
  "bbob-biobj_f36_i05_d02 0.771461442665405",
  "bbob-biobj_f36_i05_d03 0.852780214876783",
  "bbob-biobj_f36_i05_d05 0.954690339029068",
  "bbob-biobj_f36_i05_d10 0.808803658664174",
  "bbob-biobj_f36_i05_d20 0.865931157268564",
  "bbob-biobj_f36_i05_d40 0.817355203925663",
  "bbob-biobj_f36_i06_d02 0.933695763707815",
  "bbob-biobj_f36_i06_d03 0.801432031610929",
  "bbob-biobj_f36_i06_d05 0.830708343034619",
  "bbob-biobj_f36_i06_d10 0.753132543716327",
  "bbob-biobj_f36_i06_d20 0.886713265290704",
  "bbob-biobj_f36_i06_d40 0.845262138248134",
  "bbob-biobj_f36_i07_d02 0.937657704240961",
  "bbob-biobj_f36_i07_d03 0.733614515844247",
  "bbob-biobj_f36_i07_d05 0.892540846527257",
  "bbob-biobj_f36_i07_d10 0.726190913948427",
  "bbob-biobj_f36_i07_d20 0.877578776300129",
  "bbob-biobj_f36_i07_d40 0.792634847065283",
  "bbob-biobj_f36_i08_d02 0.773791972054986",
  "bbob-biobj_f36_i08_d03 0.582113784704423",
  "bbob-biobj_f36_i08_d05 0.935187367108437",
  "bbob-biobj_f36_i08_d10 0.876695021233218",
  "bbob-biobj_f36_i08_d20 0.824170206905267",
  "bbob-biobj_f36_i08_d40 0.865436409580512",
  "bbob-biobj_f36_i09_d02 0.891444356002966",
  "bbob-biobj_f36_i09_d03 0.625146222039320",
  "bbob-biobj_f36_i09_d05 0.820232336176047",
  "bbob-biobj_f36_i09_d10 0.861524351755032",
  "bbob-biobj_f36_i09_d20 0.829573948076514",
  "bbob-biobj_f36_i09_d40 0.829754021006423",
  "bbob-biobj_f36_i10_d02 0.934398020485006",
  "bbob-biobj_f36_i10_d03 0.705459873049955",
  "bbob-biobj_f36_i10_d05 0.920393161168462",
  "bbob-biobj_f36_i10_d10 0.835561034075361",
  "bbob-biobj_f36_i10_d20 0.809907143893169",
  "bbob-biobj_f36_i10_d40 0.790989826102430",
  "bbob-biobj_f37_i01_d02 0.887325713230270",
  "bbob-biobj_f37_i01_d03 0.807648927271962",
  "bbob-biobj_f37_i01_d05 0.835116542946276",
  "bbob-biobj_f37_i01_d10 0.782782389790906",
  "bbob-biobj_f37_i01_d20 0.841966768537659",
  "bbob-biobj_f37_i01_d40 0.788083322221278",
  "bbob-biobj_f37_i02_d02 0.737269283128202",
  "bbob-biobj_f37_i02_d03 0.835004353558839",
  "bbob-biobj_f37_i02_d05 0.929633551394840",
  "bbob-biobj_f37_i02_d10 0.807181719873087",
  "bbob-biobj_f37_i02_d20 0.794932415582117",
  "bbob-biobj_f37_i02_d40 0.798919821141479",
  "bbob-biobj_f37_i03_d02 0.802436084742524",
  "bbob-biobj_f37_i03_d03 0.779815898243700",
  "bbob-biobj_f37_i03_d05 0.841078751272964",
  "bbob-biobj_f37_i03_d10 0.835709362627074",
  "bbob-biobj_f37_i03_d20 0.808616248982069",
  "bbob-biobj_f37_i03_d40 0.794038185072231",
  "bbob-biobj_f37_i04_d02 0.827629900837309",
  "bbob-biobj_f37_i04_d03 0.927470901892586",
  "bbob-biobj_f37_i04_d05 0.767413924627218",
  "bbob-biobj_f37_i04_d10 0.783059860821889",
  "bbob-biobj_f37_i04_d20 0.834451557654742",
  "bbob-biobj_f37_i04_d40 0.771877921104696",
  "bbob-biobj_f37_i05_d02 0.833850868984802",
  "bbob-biobj_f37_i05_d03 0.841495086876068",
  "bbob-biobj_f37_i05_d05 0.837664278842414",
  "bbob-biobj_f37_i05_d10 0.756272500772225",
  "bbob-biobj_f37_i05_d20 0.759892814816738",
  "bbob-biobj_f37_i05_d40 0.806956545241705",
  "bbob-biobj_f37_i06_d02 0.966244887841858",
  "bbob-biobj_f37_i06_d03 0.919761036568247",
  "bbob-biobj_f37_i06_d05 0.870366097652817",
  "bbob-biobj_f37_i06_d10 0.802622600305857",
  "bbob-biobj_f37_i06_d20 0.787095683904922",
  "bbob-biobj_f37_i06_d40 0.822546942723380",
  "bbob-biobj_f37_i07_d02 0.865916213468063",
  "bbob-biobj_f37_i07_d03 0.887985521597400",
  "bbob-biobj_f37_i07_d05 0.852922047677374",
  "bbob-biobj_f37_i07_d10 0.768755196958351",
  "bbob-biobj_f37_i07_d20 0.755987827900339",
  "bbob-biobj_f37_i07_d40 0.760794475414749",
  "bbob-biobj_f37_i08_d02 0.954343889366429",
  "bbob-biobj_f37_i08_d03 0.779614728918267",
  "bbob-biobj_f37_i08_d05 0.866445305071004",
  "bbob-biobj_f37_i08_d10 0.822493279926153",
  "bbob-biobj_f37_i08_d20 0.853197609963602",
  "bbob-biobj_f37_i08_d40 0.733237656548855",
  "bbob-biobj_f37_i09_d02 0.651087460299773",
  "bbob-biobj_f37_i09_d03 0.971428844413651",
  "bbob-biobj_f37_i09_d05 0.874322221353586",
  "bbob-biobj_f37_i09_d10 0.933310217112306",
  "bbob-biobj_f37_i09_d20 0.824122184844244",
  "bbob-biobj_f37_i09_d40 0.802948224375232",
  "bbob-biobj_f37_i10_d02 0.846316671801481",
  "bbob-biobj_f37_i10_d03 0.783084272175719",
  "bbob-biobj_f37_i10_d05 0.772782281662687",
  "bbob-biobj_f37_i10_d10 0.826320652084208",
  "bbob-biobj_f37_i10_d20 0.840352163855009",
  "bbob-biobj_f37_i10_d40 0.728638839527935",
  "bbob-biobj_f38_i01_d02 0.877156626668269",
  "bbob-biobj_f38_i01_d03 0.866619456217389",
  "bbob-biobj_f38_i01_d05 0.915488889641595",
  "bbob-biobj_f38_i01_d10 0.841482318031030",
  "bbob-biobj_f38_i01_d20 0.857338087322617",
  "bbob-biobj_f38_i01_d40 0.823522621833064",
  "bbob-biobj_f38_i02_d02 0.906613082094068",
  "bbob-biobj_f38_i02_d03 0.793052232343624",
  "bbob-biobj_f38_i02_d05 0.872973219052256",
  "bbob-biobj_f38_i02_d10 0.886550516170649",
  "bbob-biobj_f38_i02_d20 0.836182460859613",
  "bbob-biobj_f38_i02_d40 0.846783720990037",
  "bbob-biobj_f38_i03_d02 0.700820137547856",
  "bbob-biobj_f38_i03_d03 0.753307482190734",
  "bbob-biobj_f38_i03_d05 0.890911636191755",
  "bbob-biobj_f38_i03_d10 0.875554141381291",
  "bbob-biobj_f38_i03_d20 0.831239834140523",
  "bbob-biobj_f38_i03_d40 0.884200336136909",
  "bbob-biobj_f38_i04_d02 0.629303687740648",
  "bbob-biobj_f38_i04_d03 0.886832778255523",
  "bbob-biobj_f38_i04_d05 0.807433718084827",
  "bbob-biobj_f38_i04_d10 0.886901305249582",
  "bbob-biobj_f38_i04_d20 0.898734518428570",
  "bbob-biobj_f38_i04_d40 0.860570031987039",
  "bbob-biobj_f38_i05_d02 0.802625415847500",
  "bbob-biobj_f38_i05_d03 0.930557240552683",
  "bbob-biobj_f38_i05_d05 0.906782638288660",
  "bbob-biobj_f38_i05_d10 0.873231647975663",
  "bbob-biobj_f38_i05_d20 0.878009332563302",
  "bbob-biobj_f38_i05_d40 0.853187105900135",
  "bbob-biobj_f38_i06_d02 0.853392565726335",
  "bbob-biobj_f38_i06_d03 0.903385341685126",
  "bbob-biobj_f38_i06_d05 0.797334587429967",
  "bbob-biobj_f38_i06_d10 0.828562611374134",
  "bbob-biobj_f38_i06_d20 0.898276774836316",
  "bbob-biobj_f38_i06_d40 0.855239203550405",
  "bbob-biobj_f38_i07_d02 0.737880353168575",
  "bbob-biobj_f38_i07_d03 0.949764351104347",
  "bbob-biobj_f38_i07_d05 0.865168292468046",
  "bbob-biobj_f38_i07_d10 0.906865989549259",
  "bbob-biobj_f38_i07_d20 0.869190970507406",
  "bbob-biobj_f38_i07_d40 0.844184254986462",
  "bbob-biobj_f38_i08_d02 0.935873070332340",
  "bbob-biobj_f38_i08_d03 0.838024245598931",
  "bbob-biobj_f38_i08_d05 0.898952467677852",
  "bbob-biobj_f38_i08_d10 0.872328382393700",
  "bbob-biobj_f38_i08_d20 0.853418746269193",
  "bbob-biobj_f38_i08_d40 0.776680480581279",
  "bbob-biobj_f38_i09_d02 0.629744163160619",
  "bbob-biobj_f38_i09_d03 0.942847439352580",
  "bbob-biobj_f38_i09_d05 0.838258676970166",
  "bbob-biobj_f38_i09_d10 0.766139777609307",
  "bbob-biobj_f38_i09_d20 0.845664544018467",
  "bbob-biobj_f38_i09_d40 0.825582404693079",
  "bbob-biobj_f38_i10_d02 0.861682224934748",
  "bbob-biobj_f38_i10_d03 0.907576125698610",
  "bbob-biobj_f38_i10_d05 0.947145965012375",
  "bbob-biobj_f38_i10_d10 0.933461200496374",
  "bbob-biobj_f38_i10_d20 0.896332331695701",
  "bbob-biobj_f38_i10_d40 0.870704202678407",
  "bbob-biobj_f39_i01_d02 0.904223494801458",
  "bbob-biobj_f39_i01_d03 0.853571617133106",
  "bbob-biobj_f39_i01_d05 0.891199716447439",
  "bbob-biobj_f39_i01_d10 0.980200073955946",
  "bbob-biobj_f39_i01_d20 0.867107099131664",
  "bbob-biobj_f39_i01_d40 0.840115962430494",
  "bbob-biobj_f39_i02_d02 0.983637077696188",
  "bbob-biobj_f39_i02_d03 0.963321106865497",
  "bbob-biobj_f39_i02_d05 0.979062405209409",
  "bbob-biobj_f39_i02_d10 0.937084496194365",
  "bbob-biobj_f39_i02_d20 0.927030925668967",
  "bbob-biobj_f39_i02_d40 0.872942659245748",
  "bbob-biobj_f39_i03_d02 0.979200985541726",
  "bbob-biobj_f39_i03_d03 0.842273782691279",
  "bbob-biobj_f39_i03_d05 0.903036610373556",
  "bbob-biobj_f39_i03_d10 0.857794499682875",
  "bbob-biobj_f39_i03_d20 0.868161430751310",
  "bbob-biobj_f39_i03_d40 0.881447097184690",
  "bbob-biobj_f39_i04_d02 0.978307294722389",
  "bbob-biobj_f39_i04_d03 0.956662509302352",
  "bbob-biobj_f39_i04_d05 0.968083139998477",
  "bbob-biobj_f39_i04_d10 0.906008720589640",
  "bbob-biobj_f39_i04_d20 0.876856275743194",
  "bbob-biobj_f39_i04_d40 0.887116495809417",
  "bbob-biobj_f39_i05_d02 0.986275115605065",
  "bbob-biobj_f39_i05_d03 0.958678930540330",
  "bbob-biobj_f39_i05_d05 0.921022072771728",
  "bbob-biobj_f39_i05_d10 0.884156852158941",
  "bbob-biobj_f39_i05_d20 0.868517819833135",
  "bbob-biobj_f39_i05_d40 0.866832889677387",
  "bbob-biobj_f39_i06_d02 0.996948457486112",
  "bbob-biobj_f39_i06_d03 0.953887245872590",
  "bbob-biobj_f39_i06_d05 0.844846258290417",
  "bbob-biobj_f39_i06_d10 0.891270192686639",
  "bbob-biobj_f39_i06_d20 0.889756133648198",
  "bbob-biobj_f39_i06_d40 0.886808632913449",
  "bbob-biobj_f39_i07_d02 0.945526844348775",
  "bbob-biobj_f39_i07_d03 0.990668954507648",
  "bbob-biobj_f39_i07_d05 0.942872265750309",
  "bbob-biobj_f39_i07_d10 0.844329327726334",
  "bbob-biobj_f39_i07_d20 0.871107914409235",
  "bbob-biobj_f39_i07_d40 0.851535146434555",
  "bbob-biobj_f39_i08_d02 0.784348726045595",
  "bbob-biobj_f39_i08_d03 0.787060415564065",
  "bbob-biobj_f39_i08_d05 0.992810204079008",
  "bbob-biobj_f39_i08_d10 0.810554149588323",
  "bbob-biobj_f39_i08_d20 0.849749605069192",
  "bbob-biobj_f39_i08_d40 0.852048361598631",
  "bbob-biobj_f39_i09_d02 0.995633482750447",
  "bbob-biobj_f39_i09_d03 0.929850476078642",
  "bbob-biobj_f39_i09_d05 0.880374822799933",
  "bbob-biobj_f39_i09_d10 0.992960401843265",
  "bbob-biobj_f39_i09_d20 0.902054192298907",
  "bbob-biobj_f39_i09_d40 0.874701425214485",
  "bbob-biobj_f39_i10_d02 0.733439797469165",
  "bbob-biobj_f39_i10_d03 0.810807877559873",
  "bbob-biobj_f39_i10_d05 0.906253128225098",
  "bbob-biobj_f39_i10_d10 0.864295076907232",
  "bbob-biobj_f39_i10_d20 0.870522053893872",
  "bbob-biobj_f39_i10_d40 0.862963349090832",
  "bbob-biobj_f40_i01_d02 0.800161246151599",
  "bbob-biobj_f40_i01_d03 0.914339524602664",
  "bbob-biobj_f40_i01_d05 0.903887232062397",
  "bbob-biobj_f40_i01_d10 0.728358823604654",
  "bbob-biobj_f40_i01_d20 0.603721990278000",
  "bbob-biobj_f40_i01_d40 0.454521176784640",
  "bbob-biobj_f40_i02_d02 0.814422882681879",
  "bbob-biobj_f40_i02_d03 0.684383727076659",
  "bbob-biobj_f40_i02_d05 0.823840701057669",
  "bbob-biobj_f40_i02_d10 0.744515559478456",
  "bbob-biobj_f40_i02_d20 0.644532905400418",
  "bbob-biobj_f40_i02_d40 0.475346750819245",
  "bbob-biobj_f40_i03_d02 0.841353974947936",
  "bbob-biobj_f40_i03_d03 0.938506394239854",
  "bbob-biobj_f40_i03_d05 0.809171650124029",
  "bbob-biobj_f40_i03_d10 0.729655032376510",
  "bbob-biobj_f40_i03_d20 0.526989228421848",
  "bbob-biobj_f40_i03_d40 0.594374630795347",
  "bbob-biobj_f40_i04_d02 0.545131553510016",
  "bbob-biobj_f40_i04_d03 0.820235713731847",
  "bbob-biobj_f40_i04_d05 0.914038330213304",
  "bbob-biobj_f40_i04_d10 0.707823424172015",
  "bbob-biobj_f40_i04_d20 0.748125315493907",
  "bbob-biobj_f40_i04_d40 0.554603485130928",
  "bbob-biobj_f40_i05_d02 0.710069762992855",
  "bbob-biobj_f40_i05_d03 0.761268399219076",
  "bbob-biobj_f40_i05_d05 0.946422296430662",
  "bbob-biobj_f40_i05_d10 0.738272073645916",
  "bbob-biobj_f40_i05_d20 0.618994719986882",
  "bbob-biobj_f40_i05_d40 0.524840592835769",
  "bbob-biobj_f40_i06_d02 0.940537942098922",
  "bbob-biobj_f40_i06_d03 0.965405293118817",
  "bbob-biobj_f40_i06_d05 0.873564661220253",
  "bbob-biobj_f40_i06_d10 0.712531441048768",
  "bbob-biobj_f40_i06_d20 0.837199760651909",
  "bbob-biobj_f40_i06_d40 0.558580661893026",
  "bbob-biobj_f40_i07_d02 0.911244402769808",
  "bbob-biobj_f40_i07_d03 0.831450960747717",
  "bbob-biobj_f40_i07_d05 0.906051359087307",
  "bbob-biobj_f40_i07_d10 0.785479527011245",
  "bbob-biobj_f40_i07_d20 0.631745363199700",
  "bbob-biobj_f40_i07_d40 0.601868801243822",
  "bbob-biobj_f40_i08_d02 0.931069806304358",
  "bbob-biobj_f40_i08_d03 0.960423548819111",
  "bbob-biobj_f40_i08_d05 0.865330316602925",
  "bbob-biobj_f40_i08_d10 0.809330752612949",
  "bbob-biobj_f40_i08_d20 0.646524483880005",
  "bbob-biobj_f40_i08_d40 0.591124207245141",
  "bbob-biobj_f40_i09_d02 0.929663089378316",
  "bbob-biobj_f40_i09_d03 0.706558047941276",
  "bbob-biobj_f40_i09_d05 0.846284534270260",
  "bbob-biobj_f40_i09_d10 0.873322034872999",
  "bbob-biobj_f40_i09_d20 0.678589733695060",
  "bbob-biobj_f40_i09_d40 0.543244841386772",
  "bbob-biobj_f40_i10_d02 0.978675048579028",
  "bbob-biobj_f40_i10_d03 0.932468687339708",
  "bbob-biobj_f40_i10_d05 0.837204852779286",
  "bbob-biobj_f40_i10_d10 0.834994341921121",
  "bbob-biobj_f40_i10_d20 0.720739247108300",
  "bbob-biobj_f40_i10_d40 0.548003839020667",
  "bbob-biobj_f41_i01_d02 0.822032636632435",
  "bbob-biobj_f41_i01_d03 0.885735589311779",
  "bbob-biobj_f41_i01_d05 0.975607309655536",
  "bbob-biobj_f41_i01_d10 0.976590634814697",
  "bbob-biobj_f41_i01_d20 0.971524548628439",
  "bbob-biobj_f41_i01_d40 0.965095266817710",
  "bbob-biobj_f41_i02_d02 0.914731719910576",
  "bbob-biobj_f41_i02_d03 0.923050990915773",
  "bbob-biobj_f41_i02_d05 0.924064779110605",
  "bbob-biobj_f41_i02_d10 0.931063386797241",
  "bbob-biobj_f41_i02_d20 0.981604969515710",
  "bbob-biobj_f41_i02_d40 0.951185152059880",
  "bbob-biobj_f41_i03_d02 0.853587389647327",
  "bbob-biobj_f41_i03_d03 0.914426287925777",
  "bbob-biobj_f41_i03_d05 0.974027750613517",
  "bbob-biobj_f41_i03_d10 0.871998139827258",
  "bbob-biobj_f41_i03_d20 0.941017702205533",
  "bbob-biobj_f41_i03_d40 0.954698190937639",
  "bbob-biobj_f41_i04_d02 0.512373025230438",
  "bbob-biobj_f41_i04_d03 0.695235231899295",
  "bbob-biobj_f41_i04_d05 0.889982461168974",
  "bbob-biobj_f41_i04_d10 0.951230829582214",
  "bbob-biobj_f41_i04_d20 0.913069182377331",
  "bbob-biobj_f41_i04_d40 0.980631033346320",
  "bbob-biobj_f41_i05_d02 0.821557917170107",
  "bbob-biobj_f41_i05_d03 0.822870004338797",
  "bbob-biobj_f41_i05_d05 0.876465115782367",
  "bbob-biobj_f41_i05_d10 0.966123487009661",
  "bbob-biobj_f41_i05_d20 0.924557082188807",
  "bbob-biobj_f41_i05_d40 0.951467783482356",
  "bbob-biobj_f41_i06_d02 0.767820247195123",
  "bbob-biobj_f41_i06_d03 0.896761247699564",
  "bbob-biobj_f41_i06_d05 0.792613675800985",
  "bbob-biobj_f41_i06_d10 0.893668404798527",
  "bbob-biobj_f41_i06_d20 0.963766756822212",
  "bbob-biobj_f41_i06_d40 0.957017816917340",
  "bbob-biobj_f41_i07_d02 0.935895666858192",
  "bbob-biobj_f41_i07_d03 0.915091613061218",
  "bbob-biobj_f41_i07_d05 0.922576398819569",
  "bbob-biobj_f41_i07_d10 0.850414058072923",
  "bbob-biobj_f41_i07_d20 0.901408374427671",
  "bbob-biobj_f41_i07_d40 0.977968591681356",
  "bbob-biobj_f41_i08_d02 0.948589453594646",
  "bbob-biobj_f41_i08_d03 0.818734428396040",
  "bbob-biobj_f41_i08_d05 0.967929065984017",
  "bbob-biobj_f41_i08_d10 0.823532159109795",
  "bbob-biobj_f41_i08_d20 0.958596440686672",
  "bbob-biobj_f41_i08_d40 0.945371539846042",
  "bbob-biobj_f41_i09_d02 0.784304762888349",
  "bbob-biobj_f41_i09_d03 0.946467961231456",
  "bbob-biobj_f41_i09_d05 0.824284297515830",
  "bbob-biobj_f41_i09_d10 0.938655860886718",
  "bbob-biobj_f41_i09_d20 0.980243032061870",
  "bbob-biobj_f41_i09_d40 0.949075055106026",
  "bbob-biobj_f41_i10_d02 0.636563178897354",
  "bbob-biobj_f41_i10_d03 0.839098772458601",
  "bbob-biobj_f41_i10_d05 0.911442757408109",
  "bbob-biobj_f41_i10_d10 0.875904631484907",
  "bbob-biobj_f41_i10_d20 0.976027130384361",
  "bbob-biobj_f41_i10_d40 0.952719059468594",
  "bbob-biobj_f42_i01_d02 0.948464765439568",
  "bbob-biobj_f42_i01_d03 0.964392985392458",
  "bbob-biobj_f42_i01_d05 0.935338815233537",
  "bbob-biobj_f42_i01_d10 0.942313119107756",
  "bbob-biobj_f42_i01_d20 0.963287312759796",
  "bbob-biobj_f42_i01_d40 0.965464047893837",
  "bbob-biobj_f42_i02_d02 0.938889926928718",
  "bbob-biobj_f42_i02_d03 0.896773013229354",
  "bbob-biobj_f42_i02_d05 0.918175512753701",
  "bbob-biobj_f42_i02_d10 0.932002927892709",
  "bbob-biobj_f42_i02_d20 0.932418928420582",
  "bbob-biobj_f42_i02_d40 0.929278847416208",
  "bbob-biobj_f42_i03_d02 0.813087200884499",
  "bbob-biobj_f42_i03_d03 0.936775808136482",
  "bbob-biobj_f42_i03_d05 0.860961914736147",
  "bbob-biobj_f42_i03_d10 0.911851048568690",
  "bbob-biobj_f42_i03_d20 0.939517430553683",
  "bbob-biobj_f42_i03_d40 0.937961190849771",
  "bbob-biobj_f42_i04_d02 0.882605831962261",
  "bbob-biobj_f42_i04_d03 0.956285623122566",
  "bbob-biobj_f42_i04_d05 0.929495427843329",
  "bbob-biobj_f42_i04_d10 0.870128731781709",
  "bbob-biobj_f42_i04_d20 0.951252475282381",
  "bbob-biobj_f42_i04_d40 0.917443861858564",
  "bbob-biobj_f42_i05_d02 0.980207491275489",
  "bbob-biobj_f42_i05_d03 0.933371934973166",
  "bbob-biobj_f42_i05_d05 0.970696599044961",
  "bbob-biobj_f42_i05_d10 0.913819853641673",
  "bbob-biobj_f42_i05_d20 0.955211173195455",
  "bbob-biobj_f42_i05_d40 0.952773557181925",
  "bbob-biobj_f42_i06_d02 0.948820967284496",
  "bbob-biobj_f42_i06_d03 0.886409022826839",
  "bbob-biobj_f42_i06_d05 0.921887986190143",
  "bbob-biobj_f42_i06_d10 0.944002279458475",
  "bbob-biobj_f42_i06_d20 0.948334654987018",
  "bbob-biobj_f42_i06_d40 0.922848340410259",
  "bbob-biobj_f42_i07_d02 0.973608788739378",
  "bbob-biobj_f42_i07_d03 0.911786597473834",
  "bbob-biobj_f42_i07_d05 0.953033195999899",
  "bbob-biobj_f42_i07_d10 0.894630196105259",
  "bbob-biobj_f42_i07_d20 0.965678396903687",
  "bbob-biobj_f42_i07_d40 0.938365462110054",
  "bbob-biobj_f42_i08_d02 0.917184749486420",
  "bbob-biobj_f42_i08_d03 0.965360550557004",
  "bbob-biobj_f42_i08_d05 0.943830127347422",
  "bbob-biobj_f42_i08_d10 0.951858228370568",
  "bbob-biobj_f42_i08_d20 0.955440555788915",
  "bbob-biobj_f42_i08_d40 0.904629808244290",
  "bbob-biobj_f42_i09_d02 0.966918705968638",
  "bbob-biobj_f42_i09_d03 0.964375591101980",
  "bbob-biobj_f42_i09_d05 0.912821333243006",
  "bbob-biobj_f42_i09_d10 0.926576853188507",
  "bbob-biobj_f42_i09_d20 0.961968306163202",
  "bbob-biobj_f42_i09_d40 0.903165659089824",
  "bbob-biobj_f42_i10_d02 0.968031865555333",
  "bbob-biobj_f42_i10_d03 0.957780257103254",
  "bbob-biobj_f42_i10_d05 0.918567004038035",
  "bbob-biobj_f42_i10_d10 0.945813733091541",
  "bbob-biobj_f42_i10_d20 0.966485720634660",
  "bbob-biobj_f42_i10_d40 0.869654175877484",
  "bbob-biobj_f43_i01_d02 0.806235783738648",
  "bbob-biobj_f43_i01_d03 0.961024116402951",
  "bbob-biobj_f43_i01_d05 0.898063277892934",
  "bbob-biobj_f43_i01_d10 0.993130440281062",
  "bbob-biobj_f43_i01_d20 0.938157853924559",
  "bbob-biobj_f43_i01_d40 0.983625753742276",
  "bbob-biobj_f43_i02_d02 0.928981100867030",
  "bbob-biobj_f43_i02_d03 0.900487190285260",
  "bbob-biobj_f43_i02_d05 0.953527657020681",
  "bbob-biobj_f43_i02_d10 0.943658984527069",
  "bbob-biobj_f43_i02_d20 0.955504681503142",
  "bbob-biobj_f43_i02_d40 0.975305093191710",
  "bbob-biobj_f43_i03_d02 0.703325488845851",
  "bbob-biobj_f43_i03_d03 0.738232261661174",
  "bbob-biobj_f43_i03_d05 0.988538904089657",
  "bbob-biobj_f43_i03_d10 0.982310546730188",
  "bbob-biobj_f43_i03_d20 0.971214544046469",
  "bbob-biobj_f43_i03_d40 0.943552964204884",
  "bbob-biobj_f43_i04_d02 0.955511249377773",
  "bbob-biobj_f43_i04_d03 0.898286408731989",
  "bbob-biobj_f43_i04_d05 0.945639505462955",
  "bbob-biobj_f43_i04_d10 0.962138222035042",
  "bbob-biobj_f43_i04_d20 0.965366859326638",
  "bbob-biobj_f43_i04_d40 0.966679911412569",
  "bbob-biobj_f43_i05_d02 0.527808634808338",
  "bbob-biobj_f43_i05_d03 0.967868516531061",
  "bbob-biobj_f43_i05_d05 0.902519492505864",
  "bbob-biobj_f43_i05_d10 0.968761193793143",
  "bbob-biobj_f43_i05_d20 0.978667811873788",
  "bbob-biobj_f43_i05_d40 0.965965120868660",
  "bbob-biobj_f43_i06_d02 0.913912796434573",
  "bbob-biobj_f43_i06_d03 0.911926671622942",
  "bbob-biobj_f43_i06_d05 0.968767671619481",
  "bbob-biobj_f43_i06_d10 0.930264654298930",
  "bbob-biobj_f43_i06_d20 0.959755090580140",
  "bbob-biobj_f43_i06_d40 0.916689640910848",
  "bbob-biobj_f43_i07_d02 0.901275709731305",
  "bbob-biobj_f43_i07_d03 0.881589817781074",
  "bbob-biobj_f43_i07_d05 0.959692328763792",
  "bbob-biobj_f43_i07_d10 0.984558957382070",
  "bbob-biobj_f43_i07_d20 0.967726668640305",
  "bbob-biobj_f43_i07_d40 0.948326276024638",
  "bbob-biobj_f43_i08_d02 0.701507552324083",
  "bbob-biobj_f43_i08_d03 0.967153417238951",
  "bbob-biobj_f43_i08_d05 0.941758666673685",
  "bbob-biobj_f43_i08_d10 0.986279003061626",
  "bbob-biobj_f43_i08_d20 0.974234017197797",
  "bbob-biobj_f43_i08_d40 0.955731037265340",
  "bbob-biobj_f43_i09_d02 0.862699790591577",
  "bbob-biobj_f43_i09_d03 0.862226281469927",
  "bbob-biobj_f43_i09_d05 0.989087524544532",
  "bbob-biobj_f43_i09_d10 0.978448520104203",
  "bbob-biobj_f43_i09_d20 0.971398821311664",
  "bbob-biobj_f43_i09_d40 0.853607336041518",
  "bbob-biobj_f43_i10_d02 0.848723462138902",
  "bbob-biobj_f43_i10_d03 0.819071040736239",
  "bbob-biobj_f43_i10_d05 0.931905908610244",
  "bbob-biobj_f43_i10_d10 0.958412523209058",
  "bbob-biobj_f43_i10_d20 0.981324831426623",
  "bbob-biobj_f43_i10_d40 0.964482582666160",
  "bbob-biobj_f44_i01_d02 0.990000542160782",
  "bbob-biobj_f44_i01_d03 0.989550390998165",
  "bbob-biobj_f44_i01_d05 0.993557132250125",
  "bbob-biobj_f44_i01_d10 0.977881273349498",
  "bbob-biobj_f44_i01_d20 0.989485461259344",
  "bbob-biobj_f44_i01_d40 0.992098877272136",
  "bbob-biobj_f44_i02_d02 0.996682440463554",
  "bbob-biobj_f44_i02_d03 0.744914592237804",
  "bbob-biobj_f44_i02_d05 0.951280292741821",
  "bbob-biobj_f44_i02_d10 0.975035300671498",
  "bbob-biobj_f44_i02_d20 0.984826262971341",
  "bbob-biobj_f44_i02_d40 0.981700455019443",
  "bbob-biobj_f44_i03_d02 0.982516595343871",
  "bbob-biobj_f44_i03_d03 0.988829949447386",
  "bbob-biobj_f44_i03_d05 0.958338343082007",
  "bbob-biobj_f44_i03_d10 0.978362251335445",
  "bbob-biobj_f44_i03_d20 0.977806363250877",
  "bbob-biobj_f44_i03_d40 0.981998615266170",
  "bbob-biobj_f44_i04_d02 0.975372267874977",
  "bbob-biobj_f44_i04_d03 0.927813089410043",
  "bbob-biobj_f44_i04_d05 0.965363937918500",
  "bbob-biobj_f44_i04_d10 0.987130916732076",
  "bbob-biobj_f44_i04_d20 0.984397834851993",
  "bbob-biobj_f44_i04_d40 0.980531891117153",
  "bbob-biobj_f44_i05_d02 0.992233423900571",
  "bbob-biobj_f44_i05_d03 0.972160548026576",
  "bbob-biobj_f44_i05_d05 0.966251845771374",
  "bbob-biobj_f44_i05_d10 0.972825160494269",
  "bbob-biobj_f44_i05_d20 0.974157531377812",
  "bbob-biobj_f44_i05_d40 0.981259160284495",
  "bbob-biobj_f44_i06_d02 0.955939797523231",
  "bbob-biobj_f44_i06_d03 0.938295348504735",
  "bbob-biobj_f44_i06_d05 0.926821904154388",
  "bbob-biobj_f44_i06_d10 0.970956651310300",
  "bbob-biobj_f44_i06_d20 0.981296651692079",
  "bbob-biobj_f44_i06_d40 0.989000723373714",
  "bbob-biobj_f44_i07_d02 0.982839357073339",
  "bbob-biobj_f44_i07_d03 0.968203446956965",
  "bbob-biobj_f44_i07_d05 0.996049574426312",
  "bbob-biobj_f44_i07_d10 0.961740169706816",
  "bbob-biobj_f44_i07_d20 0.936886232628911",
  "bbob-biobj_f44_i07_d40 0.987842905794241",
  "bbob-biobj_f44_i08_d02 0.881869872387817",
  "bbob-biobj_f44_i08_d03 0.950774552348263",
  "bbob-biobj_f44_i08_d05 0.982803359548557",
  "bbob-biobj_f44_i08_d10 0.976113895826915",
  "bbob-biobj_f44_i08_d20 0.985341041087579",
  "bbob-biobj_f44_i08_d40 0.981930570095902",
  "bbob-biobj_f44_i09_d02 0.974715179294778",
  "bbob-biobj_f44_i09_d03 0.993230798562239",
  "bbob-biobj_f44_i09_d05 0.952929260309795",
  "bbob-biobj_f44_i09_d10 0.966603115844480",
  "bbob-biobj_f44_i09_d20 0.996510158321955",
  "bbob-biobj_f44_i09_d40 0.991857743990371",
  "bbob-biobj_f44_i10_d02 0.742733624851043",
  "bbob-biobj_f44_i10_d03 0.788129186788466",
  "bbob-biobj_f44_i10_d05 0.997482012910640",
  "bbob-biobj_f44_i10_d10 0.974864988885511",
  "bbob-biobj_f44_i10_d20 0.961977716849415",
  "bbob-biobj_f44_i10_d40 0.984778763370486",
  "bbob-biobj_f45_i01_d02 0.951720478300025",
  "bbob-biobj_f45_i01_d03 0.828185907747463",
  "bbob-biobj_f45_i01_d05 0.917823772790742",
  "bbob-biobj_f45_i01_d10 0.962464906296622",
  "bbob-biobj_f45_i01_d20 0.916278583531648",
  "bbob-biobj_f45_i01_d40 0.734963567079671",
  "bbob-biobj_f45_i02_d02 0.797327510997148",
  "bbob-biobj_f45_i02_d03 0.976528231889640",
  "bbob-biobj_f45_i02_d05 0.892693577229959",
  "bbob-biobj_f45_i02_d10 0.937570463128827",
  "bbob-biobj_f45_i02_d20 0.852114827670767",
  "bbob-biobj_f45_i02_d40 0.841835435846432",
  "bbob-biobj_f45_i03_d02 0.729950018217370",
  "bbob-biobj_f45_i03_d03 0.953105240940278",
  "bbob-biobj_f45_i03_d05 0.953228358767808",
  "bbob-biobj_f45_i03_d10 0.939028992926788",
  "bbob-biobj_f45_i03_d20 0.945618407405583",
  "bbob-biobj_f45_i03_d40 0.870341204710171",
  "bbob-biobj_f45_i04_d02 0.962996536540628",
  "bbob-biobj_f45_i04_d03 0.954353868240249",
  "bbob-biobj_f45_i04_d05 0.946643482857920",
  "bbob-biobj_f45_i04_d10 0.932565865669521",
  "bbob-biobj_f45_i04_d20 0.929120624461619",
  "bbob-biobj_f45_i04_d40 0.920399420794777",
  "bbob-biobj_f45_i05_d02 0.924715066087407",
  "bbob-biobj_f45_i05_d03 0.942410755772471",
  "bbob-biobj_f45_i05_d05 0.909304053409049",
  "bbob-biobj_f45_i05_d10 0.875464727276631",
  "bbob-biobj_f45_i05_d20 0.913167131215061",
  "bbob-biobj_f45_i05_d40 0.851463837075312",
  "bbob-biobj_f45_i06_d02 0.904107320620183",
  "bbob-biobj_f45_i06_d03 0.926745813531225",
  "bbob-biobj_f45_i06_d05 0.908381292256803",
  "bbob-biobj_f45_i06_d10 0.974104258955199",
  "bbob-biobj_f45_i06_d20 0.838162094620567",
  "bbob-biobj_f45_i06_d40 0.844441637933453",
  "bbob-biobj_f45_i07_d02 0.436528957584281",
  "bbob-biobj_f45_i07_d03 0.967372004140929",
  "bbob-biobj_f45_i07_d05 0.948142132233281",
  "bbob-biobj_f45_i07_d10 0.959809149296999",
  "bbob-biobj_f45_i07_d20 0.868621742393358",
  "bbob-biobj_f45_i07_d40 0.800394785341641",
  "bbob-biobj_f45_i08_d02 0.780659318231633",
  "bbob-biobj_f45_i08_d03 0.986340760359067",
  "bbob-biobj_f45_i08_d05 0.971630237352987",
  "bbob-biobj_f45_i08_d10 0.895481269471874",
  "bbob-biobj_f45_i08_d20 0.904569145487528",
  "bbob-biobj_f45_i08_d40 0.935599722245083",
  "bbob-biobj_f45_i09_d02 0.717377978573964",
  "bbob-biobj_f45_i09_d03 0.983808114559429",
  "bbob-biobj_f45_i09_d05 0.661546908616817",
  "bbob-biobj_f45_i09_d10 0.941394391132751",
  "bbob-biobj_f45_i09_d20 0.879474099417815",
  "bbob-biobj_f45_i09_d40 0.858019798753395",
  "bbob-biobj_f45_i10_d02 0.963294936788939",
  "bbob-biobj_f45_i10_d03 0.911497187792734",
  "bbob-biobj_f45_i10_d05 0.888578532103580",
  "bbob-biobj_f45_i10_d10 0.921122750482157",
  "bbob-biobj_f45_i10_d20 0.919817170351647",
  "bbob-biobj_f45_i10_d40 0.849484621351490",
  "bbob-biobj_f46_i01_d02 0.761065158078483",
  "bbob-biobj_f46_i01_d03 0.903454101182485",
  "bbob-biobj_f46_i01_d05 0.947149328002746",
  "bbob-biobj_f46_i01_d10 0.916272715210604",
  "bbob-biobj_f46_i01_d20 0.942357928815136",
  "bbob-biobj_f46_i01_d40 0.849919241841063",
  "bbob-biobj_f46_i02_d02 0.848797238425957",
  "bbob-biobj_f46_i02_d03 0.962406656418374",
  "bbob-biobj_f46_i02_d05 0.902619923848109",
  "bbob-biobj_f46_i02_d10 0.906928366687421",
  "bbob-biobj_f46_i02_d20 0.965397086385498",
  "bbob-biobj_f46_i02_d40 0.868682770721595",
  "bbob-biobj_f46_i03_d02 0.924046817242088",
  "bbob-biobj_f46_i03_d03 0.899712790813471",
  "bbob-biobj_f46_i03_d05 0.893021851011440",
  "bbob-biobj_f46_i03_d10 0.891787108025778",
  "bbob-biobj_f46_i03_d20 0.923123910684272",
  "bbob-biobj_f46_i03_d40 0.886927414745710",
  "bbob-biobj_f46_i04_d02 0.934266788968986",
  "bbob-biobj_f46_i04_d03 0.892786792662613",
  "bbob-biobj_f46_i04_d05 0.909683052458784",
  "bbob-biobj_f46_i04_d10 0.906573356158082",
  "bbob-biobj_f46_i04_d20 0.940501260895213",
  "bbob-biobj_f46_i04_d40 0.944504408668998",
  "bbob-biobj_f46_i05_d02 0.925925506935702",
  "bbob-biobj_f46_i05_d03 0.846189938849445",
  "bbob-biobj_f46_i05_d05 0.942523699507068",
  "bbob-biobj_f46_i05_d10 0.898604788318887",
  "bbob-biobj_f46_i05_d20 0.905078417397695",
  "bbob-biobj_f46_i05_d40 0.920817665308800",
  "bbob-biobj_f46_i06_d02 0.860597324870132",
  "bbob-biobj_f46_i06_d03 0.908560135371942",
  "bbob-biobj_f46_i06_d05 0.898205575540547",
  "bbob-biobj_f46_i06_d10 0.943031520425248",
  "bbob-biobj_f46_i06_d20 0.932185791880173",
  "bbob-biobj_f46_i06_d40 0.894836698085377",
  "bbob-biobj_f46_i07_d02 0.839480891518788",
  "bbob-biobj_f46_i07_d03 0.924352421051578",
  "bbob-biobj_f46_i07_d05 0.885229249695310",
  "bbob-biobj_f46_i07_d10 0.932899016645717",
  "bbob-biobj_f46_i07_d20 0.929079747854824",
  "bbob-biobj_f46_i07_d40 0.837529460753000",
  "bbob-biobj_f46_i08_d02 0.936138999436885",
  "bbob-biobj_f46_i08_d03 0.934247711003668",
  "bbob-biobj_f46_i08_d05 0.918939976285668",
  "bbob-biobj_f46_i08_d10 0.937241161288771",
  "bbob-biobj_f46_i08_d20 0.897583888472514",
  "bbob-biobj_f46_i08_d40 0.870059626300228",
  "bbob-biobj_f46_i09_d02 0.883491787311581",
  "bbob-biobj_f46_i09_d03 0.960074190806033",
  "bbob-biobj_f46_i09_d05 0.904516729566851",
  "bbob-biobj_f46_i09_d10 0.913027260219081",
  "bbob-biobj_f46_i09_d20 0.916179229458673",
  "bbob-biobj_f46_i09_d40 0.828379770488438",
  "bbob-biobj_f46_i10_d02 0.881159321912029",
  "bbob-biobj_f46_i10_d03 0.911062675472851",
  "bbob-biobj_f46_i10_d05 0.939597032363025",
  "bbob-biobj_f46_i10_d10 0.936171859969650",
  "bbob-biobj_f46_i10_d20 0.953579426686093",
  "bbob-biobj_f46_i10_d40 0.863007419814745",
  "bbob-biobj_f47_i01_d02 0.712271521844314",
  "bbob-biobj_f47_i01_d03 0.868936887809205",
  "bbob-biobj_f47_i01_d05 0.942626543078939",
  "bbob-biobj_f47_i01_d10 0.956927722858880",
  "bbob-biobj_f47_i01_d20 0.953958051495361",
  "bbob-biobj_f47_i01_d40 0.897063620268465",
  "bbob-biobj_f47_i02_d02 0.939177831742199",
  "bbob-biobj_f47_i02_d03 0.954699108767530",
  "bbob-biobj_f47_i02_d05 0.930096709628886",
  "bbob-biobj_f47_i02_d10 0.905935402807431",
  "bbob-biobj_f47_i02_d20 0.953872180548399",
  "bbob-biobj_f47_i02_d40 0.896343561824198",
  "bbob-biobj_f47_i03_d02 0.739793093503317",
  "bbob-biobj_f47_i03_d03 0.961190108422431",
  "bbob-biobj_f47_i03_d05 0.976644922631438",
  "bbob-biobj_f47_i03_d10 0.947373226338933",
  "bbob-biobj_f47_i03_d20 0.944104864565450",
  "bbob-biobj_f47_i03_d40 0.909842805270730",
  "bbob-biobj_f47_i04_d02 0.779551093004040",
  "bbob-biobj_f47_i04_d03 0.939953689474088",
  "bbob-biobj_f47_i04_d05 0.922265797959458",
  "bbob-biobj_f47_i04_d10 0.953074812047075",
  "bbob-biobj_f47_i04_d20 0.946907776538079",
  "bbob-biobj_f47_i04_d40 0.958122737361373",
  "bbob-biobj_f47_i05_d02 0.944459623927879",
  "bbob-biobj_f47_i05_d03 0.910102071038915",
  "bbob-biobj_f47_i05_d05 0.870082492125726",
  "bbob-biobj_f47_i05_d10 0.976078696991155",
  "bbob-biobj_f47_i05_d20 0.918114403213010",
  "bbob-biobj_f47_i05_d40 0.959890019242779",
  "bbob-biobj_f47_i06_d02 0.844033597656672",
  "bbob-biobj_f47_i06_d03 0.968886162689234",
  "bbob-biobj_f47_i06_d05 0.912569180735884",
  "bbob-biobj_f47_i06_d10 0.955020327741348",
  "bbob-biobj_f47_i06_d20 0.959931736662928",
  "bbob-biobj_f47_i06_d40 0.863825720448060",
  "bbob-biobj_f47_i07_d02 0.868550984667518",
  "bbob-biobj_f47_i07_d03 0.828589441113759",
  "bbob-biobj_f47_i07_d05 0.965077074246411",
  "bbob-biobj_f47_i07_d10 0.979063661818018",
  "bbob-biobj_f47_i07_d20 0.989441498540136",
  "bbob-biobj_f47_i07_d40 0.888977656268294",
  "bbob-biobj_f47_i08_d02 0.958208327202276",
  "bbob-biobj_f47_i08_d03 0.913890286183394",
  "bbob-biobj_f47_i08_d05 0.907287541518252",
  "bbob-biobj_f47_i08_d10 0.950289934582031",
  "bbob-biobj_f47_i08_d20 0.956393436256001",
  "bbob-biobj_f47_i08_d40 0.865268951606849",
  "bbob-biobj_f47_i09_d02 0.889803641676326",
  "bbob-biobj_f47_i09_d03 0.960816560627750",
  "bbob-biobj_f47_i09_d05 0.940637892427626",
  "bbob-biobj_f47_i09_d10 0.936019590332102",
  "bbob-biobj_f47_i09_d20 0.938435983318132",
  "bbob-biobj_f47_i09_d40 0.842106041578083",
  "bbob-biobj_f47_i10_d02 0.645312273047757",
  "bbob-biobj_f47_i10_d03 0.896312911645273",
  "bbob-biobj_f47_i10_d05 0.902038741033290",
  "bbob-biobj_f47_i10_d10 0.951941478635096",
  "bbob-biobj_f47_i10_d20 0.952530563732015",
  "bbob-biobj_f47_i10_d40 0.813548779100876",
  "bbob-biobj_f48_i01_d02 0.848025765887428",
  "bbob-biobj_f48_i01_d03 0.973369992449986",
  "bbob-biobj_f48_i01_d05 0.970790088598286",
  "bbob-biobj_f48_i01_d10 0.974668055191848",
  "bbob-biobj_f48_i01_d20 0.972328623264621",
  "bbob-biobj_f48_i01_d40 0.967100206005242",
  "bbob-biobj_f48_i02_d02 0.994793668285615",
  "bbob-biobj_f48_i02_d03 0.870172001352784",
  "bbob-biobj_f48_i02_d05 0.991208400223761",
  "bbob-biobj_f48_i02_d10 0.987176940696567",
  "bbob-biobj_f48_i02_d20 0.972977939473435",
  "bbob-biobj_f48_i02_d40 0.950320787084731",
  "bbob-biobj_f48_i03_d02 0.974655294061295",
  "bbob-biobj_f48_i03_d03 0.980824615534915",
  "bbob-biobj_f48_i03_d05 0.992966701151600",
  "bbob-biobj_f48_i03_d10 0.979838599958335",
  "bbob-biobj_f48_i03_d20 0.972376025300435",
  "bbob-biobj_f48_i03_d40 0.963450647714541",
  "bbob-biobj_f48_i04_d02 0.991892533270107",
  "bbob-biobj_f48_i04_d03 0.977080995317503",
  "bbob-biobj_f48_i04_d05 0.979095136313054",
  "bbob-biobj_f48_i04_d10 0.970965597008950",
  "bbob-biobj_f48_i04_d20 0.984162743025940",
  "bbob-biobj_f48_i04_d40 0.961057143661875",
  "bbob-biobj_f48_i05_d02 0.988654636565983",
  "bbob-biobj_f48_i05_d03 0.980478422878336",
  "bbob-biobj_f48_i05_d05 0.973955007171175",
  "bbob-biobj_f48_i05_d10 0.987143043165127",
  "bbob-biobj_f48_i05_d20 0.963640714977916",
  "bbob-biobj_f48_i05_d40 0.973473305269489",
  "bbob-biobj_f48_i06_d02 0.389125092425358",
  "bbob-biobj_f48_i06_d03 0.992878211458856",
  "bbob-biobj_f48_i06_d05 0.964220661205509",
  "bbob-biobj_f48_i06_d10 0.982950597891665",
  "bbob-biobj_f48_i06_d20 0.969165048940038",
  "bbob-biobj_f48_i06_d40 0.938606560339671",
  "bbob-biobj_f48_i07_d02 0.986571537690663",
  "bbob-biobj_f48_i07_d03 0.984998898000340",
  "bbob-biobj_f48_i07_d05 0.987267880721301",
  "bbob-biobj_f48_i07_d10 0.981289929843155",
  "bbob-biobj_f48_i07_d20 0.993565979775790",
  "bbob-biobj_f48_i07_d40 0.923761387549727",
  "bbob-biobj_f48_i08_d02 0.954716094519366",
  "bbob-biobj_f48_i08_d03 0.953721237596988",
  "bbob-biobj_f48_i08_d05 0.968206356504483",
  "bbob-biobj_f48_i08_d10 0.974460701815274",
  "bbob-biobj_f48_i08_d20 0.986574637004161",
  "bbob-biobj_f48_i08_d40 0.933444093880587",
  "bbob-biobj_f48_i09_d02 0.957215799041695",
  "bbob-biobj_f48_i09_d03 0.952485035481476",
  "bbob-biobj_f48_i09_d05 0.971532847132109",
  "bbob-biobj_f48_i09_d10 0.991923965442916",
  "bbob-biobj_f48_i09_d20 0.978372242991204",
  "bbob-biobj_f48_i09_d40 0.921584491794025",
  "bbob-biobj_f48_i10_d02 0.993958970623558",
  "bbob-biobj_f48_i10_d03 0.983267975994625",
  "bbob-biobj_f48_i10_d05 0.968472704999628",
  "bbob-biobj_f48_i10_d10 0.994326886173225",
  "bbob-biobj_f48_i10_d20 0.987291099443830",
  "bbob-biobj_f48_i10_d40 0.914925901565964",
  "bbob-biobj_f49_i01_d02 0.926653565537568",
  "bbob-biobj_f49_i01_d03 0.957627076527394",
  "bbob-biobj_f49_i01_d05 0.956451000243124",
  "bbob-biobj_f49_i01_d10 0.924536477155924",
  "bbob-biobj_f49_i01_d20 0.882594082920566",
  "bbob-biobj_f49_i01_d40 0.805813845934873",
  "bbob-biobj_f49_i02_d02 0.961467341262251",
  "bbob-biobj_f49_i02_d03 0.946678903641252",
  "bbob-biobj_f49_i02_d05 0.955513299086273",
  "bbob-biobj_f49_i02_d10 0.849260821503450",
  "bbob-biobj_f49_i02_d20 0.839372263421586",
  "bbob-biobj_f49_i02_d40 0.749979167333866",
  "bbob-biobj_f49_i03_d02 0.848123281524619",
  "bbob-biobj_f49_i03_d03 0.942235598237515",
  "bbob-biobj_f49_i03_d05 0.964440386354545",
  "bbob-biobj_f49_i03_d10 0.927315375792253",
  "bbob-biobj_f49_i03_d20 0.867643492890123",
  "bbob-biobj_f49_i03_d40 0.770206936559363",
  "bbob-biobj_f49_i04_d02 0.934365739580161",
  "bbob-biobj_f49_i04_d03 0.927009796637233",
  "bbob-biobj_f49_i04_d05 0.968253100687801",
  "bbob-biobj_f49_i04_d10 0.878699507006033",
  "bbob-biobj_f49_i04_d20 0.872039319541837",
  "bbob-biobj_f49_i04_d40 0.903243998826101",
  "bbob-biobj_f49_i05_d02 0.928529327440040",
  "bbob-biobj_f49_i05_d03 0.899347718071056",
  "bbob-biobj_f49_i05_d05 0.951109258207328",
  "bbob-biobj_f49_i05_d10 0.875411107545713",
  "bbob-biobj_f49_i05_d20 0.839744829854665",
  "bbob-biobj_f49_i05_d40 0.838474270256212",
  "bbob-biobj_f49_i06_d02 0.786949353143911",
  "bbob-biobj_f49_i06_d03 0.938079446743783",
  "bbob-biobj_f49_i06_d05 0.887449562107843",
  "bbob-biobj_f49_i06_d10 0.958731685899055",
  "bbob-biobj_f49_i06_d20 0.915424880884173",
  "bbob-biobj_f49_i06_d40 0.744341733850764",
  "bbob-biobj_f49_i07_d02 0.725199175756735",
  "bbob-biobj_f49_i07_d03 0.976611362171364",
  "bbob-biobj_f49_i07_d05 0.954867125606800",
  "bbob-biobj_f49_i07_d10 0.948723375003504",
  "bbob-biobj_f49_i07_d20 0.851557196307848",
  "bbob-biobj_f49_i07_d40 0.751742480457290",
  "bbob-biobj_f49_i08_d02 0.959087812365026",
  "bbob-biobj_f49_i08_d03 0.979412441260057",
  "bbob-biobj_f49_i08_d05 0.874905490508062",
  "bbob-biobj_f49_i08_d10 0.932889781512857",
  "bbob-biobj_f49_i08_d20 0.908752919184370",
  "bbob-biobj_f49_i08_d40 0.682266052593558",
  "bbob-biobj_f49_i09_d02 0.970397812528866",
  "bbob-biobj_f49_i09_d03 0.933469142335403",
  "bbob-biobj_f49_i09_d05 0.979678466709802",
  "bbob-biobj_f49_i09_d10 0.914319670500001",
  "bbob-biobj_f49_i09_d20 0.809856761863152",
  "bbob-biobj_f49_i09_d40 0.756785375595970",
  "bbob-biobj_f49_i10_d02 0.914170511611676",
  "bbob-biobj_f49_i10_d03 0.956332706338877",
  "bbob-biobj_f49_i10_d05 0.977122011006727",
  "bbob-biobj_f49_i10_d10 0.969867306431658",
  "bbob-biobj_f49_i10_d20 0.890427598620254",
  "bbob-biobj_f49_i10_d40 0.841945443823362",
  "bbob-biobj_f50_i01_d02 0.908152951825516",
  "bbob-biobj_f50_i01_d03 0.969141636576372",
  "bbob-biobj_f50_i01_d05 0.971108504161016",
  "bbob-biobj_f50_i01_d10 0.935342802675479",
  "bbob-biobj_f50_i01_d20 0.967909054527284",
  "bbob-biobj_f50_i01_d40 0.963175365926732",
  "bbob-biobj_f50_i02_d02 0.951055945346928",
  "bbob-biobj_f50_i02_d03 0.967314095314753",
  "bbob-biobj_f50_i02_d05 0.954382281509502",
  "bbob-biobj_f50_i02_d10 0.960658627454367",
  "bbob-biobj_f50_i02_d20 0.968853344346544",
  "bbob-biobj_f50_i02_d40 0.961060188833321",
  "bbob-biobj_f50_i03_d02 0.921885271301822",
  "bbob-biobj_f50_i03_d03 0.912990840485907",
  "bbob-biobj_f50_i03_d05 0.961254450390365",
  "bbob-biobj_f50_i03_d10 0.934588375214548",
  "bbob-biobj_f50_i03_d20 0.958380193560993",
  "bbob-biobj_f50_i03_d40 0.929985993495247",
  "bbob-biobj_f50_i04_d02 0.908169132408745",
  "bbob-biobj_f50_i04_d03 0.826467384959616",
  "bbob-biobj_f50_i04_d05 0.972715407083074",
  "bbob-biobj_f50_i04_d10 0.979562594561223",
  "bbob-biobj_f50_i04_d20 0.986987999881198",
  "bbob-biobj_f50_i04_d40 0.948525163299196",
  "bbob-biobj_f50_i05_d02 0.808869730820323",
  "bbob-biobj_f50_i05_d03 0.990885009560591",
  "bbob-biobj_f50_i05_d05 0.994266708206341",
  "bbob-biobj_f50_i05_d10 0.974171284097960",
  "bbob-biobj_f50_i05_d20 0.981475413161458",
  "bbob-biobj_f50_i05_d40 0.974034468942063",
  "bbob-biobj_f50_i06_d02 0.917865473554499",
  "bbob-biobj_f50_i06_d03 0.981572384861056",
  "bbob-biobj_f50_i06_d05 0.949912067011369",
  "bbob-biobj_f50_i06_d10 0.925145696976871",
  "bbob-biobj_f50_i06_d20 0.983593193905868",
  "bbob-biobj_f50_i06_d40 0.945311175178354",
  "bbob-biobj_f50_i07_d02 0.900074558278857",
  "bbob-biobj_f50_i07_d03 0.917241784837036",
  "bbob-biobj_f50_i07_d05 0.955403892840706",
  "bbob-biobj_f50_i07_d10 0.915538765236820",
  "bbob-biobj_f50_i07_d20 0.964225658237811",
  "bbob-biobj_f50_i07_d40 0.895380500125037",
  "bbob-biobj_f50_i08_d02 0.947042097595082",
  "bbob-biobj_f50_i08_d03 0.903539600334574",
  "bbob-biobj_f50_i08_d05 0.920878731491008",
  "bbob-biobj_f50_i08_d10 0.950530846024484",
  "bbob-biobj_f50_i08_d20 0.963777388113597",
  "bbob-biobj_f50_i08_d40 0.915844256054755",
  "bbob-biobj_f50_i09_d02 0.864742017412034",
  "bbob-biobj_f50_i09_d03 0.833072953308811",
  "bbob-biobj_f50_i09_d05 0.969312547401378",
  "bbob-biobj_f50_i09_d10 0.956569954494617",
  "bbob-biobj_f50_i09_d20 0.976325266630043",
  "bbob-biobj_f50_i09_d40 0.917370839548581",
  "bbob-biobj_f50_i10_d02 0.900977792606790",
  "bbob-biobj_f50_i10_d03 0.786358327046805",
  "bbob-biobj_f50_i10_d05 0.838709408221832",
  "bbob-biobj_f50_i10_d10 0.964213559993532",
  "bbob-biobj_f50_i10_d20 0.977062606412686",
  "bbob-biobj_f50_i10_d40 0.919477061169469",
  "bbob-biobj_f51_i01_d02 0.940860537397499",
  "bbob-biobj_f51_i01_d03 0.868050298787463",
  "bbob-biobj_f51_i01_d05 0.986066284877236",
  "bbob-biobj_f51_i01_d10 0.979727485173911",
  "bbob-biobj_f51_i01_d20 0.987289470904479",
  "bbob-biobj_f51_i01_d40 0.961560752565228",
  "bbob-biobj_f51_i02_d02 0.909955742339090",
  "bbob-biobj_f51_i02_d03 0.990170991876764",
  "bbob-biobj_f51_i02_d05 0.980582713805348",
  "bbob-biobj_f51_i02_d10 0.990889288444959",
  "bbob-biobj_f51_i02_d20 0.993096262776686",
  "bbob-biobj_f51_i02_d40 0.971044927766033",
  "bbob-biobj_f51_i03_d02 0.941593182004721",
  "bbob-biobj_f51_i03_d03 0.933222688826134",
  "bbob-biobj_f51_i03_d05 0.972201564428101",
  "bbob-biobj_f51_i03_d10 0.968402510947740",
  "bbob-biobj_f51_i03_d20 0.983327760434045",
  "bbob-biobj_f51_i03_d40 0.976761206098427",
  "bbob-biobj_f51_i04_d02 0.986038924391174",
  "bbob-biobj_f51_i04_d03 0.947517578406138",
  "bbob-biobj_f51_i04_d05 0.983435485457280",
  "bbob-biobj_f51_i04_d10 0.997672798920720",
  "bbob-biobj_f51_i04_d20 0.988431280628644",
  "bbob-biobj_f51_i04_d40 0.987714677951497",
  "bbob-biobj_f51_i05_d02 0.971128838482503",
  "bbob-biobj_f51_i05_d03 0.983377571859782",
  "bbob-biobj_f51_i05_d05 0.994441829315396",
  "bbob-biobj_f51_i05_d10 0.965035354448681",
  "bbob-biobj_f51_i05_d20 0.983814056029309",
  "bbob-biobj_f51_i05_d40 0.980779908664141",
  "bbob-biobj_f51_i06_d02 0.911401665218032",
  "bbob-biobj_f51_i06_d03 0.972110908246713",
  "bbob-biobj_f51_i06_d05 0.980701255409435",
  "bbob-biobj_f51_i06_d10 0.980296412725697",
  "bbob-biobj_f51_i06_d20 0.978476692117109",
  "bbob-biobj_f51_i06_d40 0.976075236770952",
  "bbob-biobj_f51_i07_d02 0.993371169749195",
  "bbob-biobj_f51_i07_d03 0.992364256666899",
  "bbob-biobj_f51_i07_d05 0.981936194030115",
  "bbob-biobj_f51_i07_d10 0.989272822120050",
  "bbob-biobj_f51_i07_d20 0.977031890389840",
  "bbob-biobj_f51_i07_d40 0.966244464997337",
  "bbob-biobj_f51_i08_d02 0.985455569133304",
  "bbob-biobj_f51_i08_d03 0.760351775519736",
  "bbob-biobj_f51_i08_d05 0.957399683039877",
  "bbob-biobj_f51_i08_d10 0.994689881572957",
  "bbob-biobj_f51_i08_d20 0.984733505414459",
  "bbob-biobj_f51_i08_d40 0.962625068804169",
  "bbob-biobj_f51_i09_d02 0.942330495240780",
  "bbob-biobj_f51_i09_d03 0.993875040593506",
  "bbob-biobj_f51_i09_d05 0.974493263503671",
  "bbob-biobj_f51_i09_d10 0.980234695333247",
  "bbob-biobj_f51_i09_d20 0.989405834704547",
  "bbob-biobj_f51_i09_d40 0.902893968785303",
  "bbob-biobj_f51_i10_d02 0.969177500995289",
  "bbob-biobj_f51_i10_d03 0.981217877288816",
  "bbob-biobj_f51_i10_d05 0.956825308782293",
  "bbob-biobj_f51_i10_d10 0.975260891512735",
  "bbob-biobj_f51_i10_d20 0.978249764234613",
  "bbob-biobj_f51_i10_d40 0.972698880391288",
  "bbob-biobj_f52_i01_d02 0.947011512653082",
  "bbob-biobj_f52_i01_d03 0.887648497718161",
  "bbob-biobj_f52_i01_d05 0.961483017325633",
  "bbob-biobj_f52_i01_d10 0.877870505157748",
  "bbob-biobj_f52_i01_d20 0.946408601011959",
  "bbob-biobj_f52_i01_d40 0.959402130774812",
  "bbob-biobj_f52_i02_d02 0.786879894383408",
  "bbob-biobj_f52_i02_d03 0.975826628427212",
  "bbob-biobj_f52_i02_d05 0.978951060363674",
  "bbob-biobj_f52_i02_d10 0.908534195447815",
  "bbob-biobj_f52_i02_d20 0.951560419112876",
  "bbob-biobj_f52_i02_d40 0.747413942686227",
  "bbob-biobj_f52_i03_d02 0.931104786126480",
  "bbob-biobj_f52_i03_d03 0.944891102453365",
  "bbob-biobj_f52_i03_d05 0.981024778906876",
  "bbob-biobj_f52_i03_d10 0.927951108748541",
  "bbob-biobj_f52_i03_d20 0.887029220128141",
  "bbob-biobj_f52_i03_d40 0.929744956131365",
  "bbob-biobj_f52_i04_d02 0.719116577653030",
  "bbob-biobj_f52_i04_d03 0.963896489819758",
  "bbob-biobj_f52_i04_d05 0.990691026682001",
  "bbob-biobj_f52_i04_d10 0.991671616677661",
  "bbob-biobj_f52_i04_d20 0.921560853016768",
  "bbob-biobj_f52_i04_d40 0.847463228547690",
  "bbob-biobj_f52_i05_d02 0.781598663041217",
  "bbob-biobj_f52_i05_d03 0.970414169235087",
  "bbob-biobj_f52_i05_d05 0.985992027143371",
  "bbob-biobj_f52_i05_d10 0.880378120976968",
  "bbob-biobj_f52_i05_d20 0.936579438285795",
  "bbob-biobj_f52_i05_d40 0.908421756541198",
  "bbob-biobj_f52_i06_d02 0.642286646964650",
  "bbob-biobj_f52_i06_d03 0.884580608370422",
  "bbob-biobj_f52_i06_d05 0.965646602875562",
  "bbob-biobj_f52_i06_d10 0.781247664061266",
  "bbob-biobj_f52_i06_d20 0.846802757337477",
  "bbob-biobj_f52_i06_d40 0.859460759692182",
  "bbob-biobj_f52_i07_d02 0.921075710769415",
  "bbob-biobj_f52_i07_d03 0.981363520984398",
  "bbob-biobj_f52_i07_d05 0.934928152030306",
  "bbob-biobj_f52_i07_d10 0.943312785195896",
  "bbob-biobj_f52_i07_d20 0.929538608766954",
  "bbob-biobj_f52_i07_d40 0.860605783999671",
  "bbob-biobj_f52_i08_d02 0.993163244830518",
  "bbob-biobj_f52_i08_d03 0.935979252762347",
  "bbob-biobj_f52_i08_d05 0.968217041850876",
  "bbob-biobj_f52_i08_d10 0.901806757355025",
  "bbob-biobj_f52_i08_d20 0.890309663486097",
  "bbob-biobj_f52_i08_d40 0.910834164473842",
  "bbob-biobj_f52_i09_d02 0.975717536629409",
  "bbob-biobj_f52_i09_d03 0.946397644440578",
  "bbob-biobj_f52_i09_d05 0.969464143075948",
  "bbob-biobj_f52_i09_d10 0.932426134259369",
  "bbob-biobj_f52_i09_d20 0.986014464477212",
  "bbob-biobj_f52_i09_d40 0.830761409093344",
  "bbob-biobj_f52_i10_d02 0.873591452423777",
  "bbob-biobj_f52_i10_d03 0.975710559698658",
  "bbob-biobj_f52_i10_d05 0.965277547763757",
  "bbob-biobj_f52_i10_d10 0.944079671401443",
  "bbob-biobj_f52_i10_d20 0.921931372525067",
  "bbob-biobj_f52_i10_d40 0.833577857863187",
  "bbob-biobj_f53_i01_d02 0.997875679082022",
  "bbob-biobj_f53_i01_d03 0.999784713883165",
  "bbob-biobj_f53_i01_d05 0.999101042864677",
  "bbob-biobj_f53_i01_d10 0.997770078215682",
  "bbob-biobj_f53_i01_d20 0.993618426103621",
  "bbob-biobj_f53_i01_d40 0.997007591133070",
  "bbob-biobj_f53_i02_d02 0.975730723529639",
  "bbob-biobj_f53_i02_d03 0.981828782563799",
  "bbob-biobj_f53_i02_d05 0.998939988185228",
  "bbob-biobj_f53_i02_d10 0.995011747816985",
  "bbob-biobj_f53_i02_d20 0.987970981006518",
  "bbob-biobj_f53_i02_d40 0.997891707970316",
  "bbob-biobj_f53_i03_d02 0.975730722493676",
  "bbob-biobj_f53_i03_d03 0.981828812869113",
  "bbob-biobj_f53_i03_d05 0.997632018309794",
  "bbob-biobj_f53_i03_d10 0.999758269430491",
  "bbob-biobj_f53_i03_d20 0.998225088125355",
  "bbob-biobj_f53_i03_d40 0.999613973046764",
  "bbob-biobj_f53_i04_d02 0.997875679038962",
  "bbob-biobj_f53_i04_d03 0.999784713869230",
  "bbob-biobj_f53_i04_d05 0.999956638376693",
  "bbob-biobj_f53_i04_d10 0.998086536338641",
  "bbob-biobj_f53_i04_d20 0.998920776357780",
  "bbob-biobj_f53_i04_d40 0.998671059364190",
  "bbob-biobj_f53_i05_d02 0.706100997334071",
  "bbob-biobj_f53_i05_d03 0.998471496817153",
  "bbob-biobj_f53_i05_d05 0.998401673001044",
  "bbob-biobj_f53_i05_d10 0.999964504713894",
  "bbob-biobj_f53_i05_d20 0.995340213985215",
  "bbob-biobj_f53_i05_d40 0.992051098660953",
  "bbob-biobj_f53_i06_d02 0.975724821634887",
  "bbob-biobj_f53_i06_d03 0.981817938727238",
  "bbob-biobj_f53_i06_d05 0.999917387947711",
  "bbob-biobj_f53_i06_d10 0.999345583316699",
  "bbob-biobj_f53_i06_d20 0.991804521300823",
  "bbob-biobj_f53_i06_d40 0.994145338593300",
  "bbob-biobj_f53_i07_d02 0.975724694195531",
  "bbob-biobj_f53_i07_d03 0.981813830814595",
  "bbob-biobj_f53_i07_d05 0.997597927318593",
  "bbob-biobj_f53_i07_d10 0.991068298815472",
  "bbob-biobj_f53_i07_d20 0.994940642699376",
  "bbob-biobj_f53_i07_d40 0.997125746501427",
  "bbob-biobj_f53_i08_d02 0.975724655180667",
  "bbob-biobj_f53_i08_d03 0.981820944797427",
  "bbob-biobj_f53_i08_d05 0.986877776861968",
  "bbob-biobj_f53_i08_d10 0.992853938058744",
  "bbob-biobj_f53_i08_d20 0.991322204564829",
  "bbob-biobj_f53_i08_d40 0.993136892638962",
  "bbob-biobj_f53_i09_d02 0.997873443973818",
  "bbob-biobj_f53_i09_d03 0.999784485681346",
  "bbob-biobj_f53_i09_d05 0.999919434752277",
  "bbob-biobj_f53_i09_d10 0.999933777616880",
  "bbob-biobj_f53_i09_d20 0.997805428716612",
  "bbob-biobj_f53_i09_d40 0.989571516180620",
  "bbob-biobj_f53_i10_d02 0.705911911234661",
  "bbob-biobj_f53_i10_d03 0.998461587303110",
  "bbob-biobj_f53_i10_d05 0.996078454979063",
  "bbob-biobj_f53_i10_d10 0.999910197013833",
  "bbob-biobj_f53_i10_d20 0.998276887563269",
  "bbob-biobj_f53_i10_d40 0.996954813573633",
  "bbob-biobj_f54_i01_d02 0.943563863449301",
  "bbob-biobj_f54_i01_d03 0.975179056969319",
  "bbob-biobj_f54_i01_d05 0.986735380916194",
  "bbob-biobj_f54_i01_d10 0.971013682474010",
  "bbob-biobj_f54_i01_d20 0.975745421404581",
  "bbob-biobj_f54_i01_d40 0.940616895149468",
  "bbob-biobj_f54_i02_d02 0.929274293624055",
  "bbob-biobj_f54_i02_d03 0.971505408191091",
  "bbob-biobj_f54_i02_d05 0.988596144676046",
  "bbob-biobj_f54_i02_d10 0.979014651689302",
  "bbob-biobj_f54_i02_d20 0.978356607421568",
  "bbob-biobj_f54_i02_d40 0.944223762762190",
  "bbob-biobj_f54_i03_d02 0.991716129312491",
  "bbob-biobj_f54_i03_d03 0.991021292919039",
  "bbob-biobj_f54_i03_d05 0.990016831732568",
  "bbob-biobj_f54_i03_d10 0.971480686910145",
  "bbob-biobj_f54_i03_d20 0.985264349205171",
  "bbob-biobj_f54_i03_d40 0.948420705230625",
  "bbob-biobj_f54_i04_d02 0.984075279299642",
  "bbob-biobj_f54_i04_d03 0.985246593860014",
  "bbob-biobj_f54_i04_d05 0.992621281722102",
  "bbob-biobj_f54_i04_d10 0.987231811133686",
  "bbob-biobj_f54_i04_d20 0.980207864670661",
  "bbob-biobj_f54_i04_d40 0.913194060760829",
  "bbob-biobj_f54_i05_d02 0.977343738568073",
  "bbob-biobj_f54_i05_d03 0.956276745028416",
  "bbob-biobj_f54_i05_d05 0.991382280615900",
  "bbob-biobj_f54_i05_d10 0.963568490810405",
  "bbob-biobj_f54_i05_d20 0.984342844583395",
  "bbob-biobj_f54_i05_d40 0.912686198505554",
  "bbob-biobj_f54_i06_d02 0.737473437336051",
  "bbob-biobj_f54_i06_d03 0.999162952048421",
  "bbob-biobj_f54_i06_d05 0.989405825987285",
  "bbob-biobj_f54_i06_d10 0.978475567959749",
  "bbob-biobj_f54_i06_d20 0.975527289662498",
  "bbob-biobj_f54_i06_d40 0.971522360668396",
  "bbob-biobj_f54_i07_d02 0.905796815383532",
  "bbob-biobj_f54_i07_d03 0.998320771468793",
  "bbob-biobj_f54_i07_d05 0.978739902764936",
  "bbob-biobj_f54_i07_d10 0.967180269190981",
  "bbob-biobj_f54_i07_d20 0.964892503660819",
  "bbob-biobj_f54_i07_d40 0.952180291500368",
  "bbob-biobj_f54_i08_d02 0.989649762421448",
  "bbob-biobj_f54_i08_d03 0.980799989191962",
  "bbob-biobj_f54_i08_d05 0.978153370143813",
  "bbob-biobj_f54_i08_d10 0.979724729374434",
  "bbob-biobj_f54_i08_d20 0.988106140416881",
  "bbob-biobj_f54_i08_d40 0.965649966756634",
  "bbob-biobj_f54_i09_d02 0.948164312518049",
  "bbob-biobj_f54_i09_d03 0.982131499169035",
  "bbob-biobj_f54_i09_d05 0.994863896309283",
  "bbob-biobj_f54_i09_d10 0.983267634101101",
  "bbob-biobj_f54_i09_d20 0.977302646027685",
  "bbob-biobj_f54_i09_d40 0.961099977238744",
  "bbob-biobj_f54_i10_d02 0.974462537098040",
  "bbob-biobj_f54_i10_d03 0.993467856136654",
  "bbob-biobj_f54_i10_d05 0.979690555440993",
  "bbob-biobj_f54_i10_d10 0.978756363308017",
  "bbob-biobj_f54_i10_d20 0.973963755363511",
  "bbob-biobj_f54_i10_d40 0.808174749933682",
  "bbob-biobj_f55_i01_d02 0.994911896613141",
  "bbob-biobj_f55_i01_d03 0.969815971998309",
  "bbob-biobj_f55_i01_d05 0.966512720377832",
  "bbob-biobj_f55_i01_d10 0.911758869928201",
  "bbob-biobj_f55_i01_d20 0.590031829651283",
  "bbob-biobj_f55_i01_d40 0.718238704214367",
  "bbob-biobj_f55_i02_d02 0.936711247261198",
  "bbob-biobj_f55_i02_d03 0.932686554911415",
  "bbob-biobj_f55_i02_d05 0.982278794041039",
  "bbob-biobj_f55_i02_d10 0.933498280356903",
  "bbob-biobj_f55_i02_d20 0.739149453321999",
  "bbob-biobj_f55_i02_d40 0.523852099096600",
  "bbob-biobj_f55_i03_d02 0.979743453071692",
  "bbob-biobj_f55_i03_d03 0.960030204113943",
  "bbob-biobj_f55_i03_d05 0.964540824392663",
  "bbob-biobj_f55_i03_d10 0.953339798554894",
  "bbob-biobj_f55_i03_d20 0.629103942577891",
  "bbob-biobj_f55_i03_d40 0.261538385131529",
  "bbob-biobj_f55_i04_d02 0.931559634623346",
  "bbob-biobj_f55_i04_d03 0.978791005218763",
  "bbob-biobj_f55_i04_d05 0.987526210974901",
  "bbob-biobj_f55_i04_d10 0.853198355715606",
  "bbob-biobj_f55_i04_d20 0.673393030763764",
  "bbob-biobj_f55_i04_d40 0.386301095523013",
  "bbob-biobj_f55_i05_d02 0.890941438964542",
  "bbob-biobj_f55_i05_d03 0.890296554935423",
  "bbob-biobj_f55_i05_d05 0.972308467317385",
  "bbob-biobj_f55_i05_d10 0.932460303729344",
  "bbob-biobj_f55_i05_d20 0.831108627152462",
  "bbob-biobj_f55_i05_d40 0.613314768304081",
  "bbob-biobj_f55_i06_d02 0.878417777547335",
  "bbob-biobj_f55_i06_d03 0.907668201380734",
  "bbob-biobj_f55_i06_d05 0.960494096606577",
  "bbob-biobj_f55_i06_d10 0.859973719860841",
  "bbob-biobj_f55_i06_d20 0.632704252867494",
  "bbob-biobj_f55_i06_d40 0.318034425030417",
  "bbob-biobj_f55_i07_d02 0.938747642930046",
  "bbob-biobj_f55_i07_d03 0.971674798747285",
  "bbob-biobj_f55_i07_d05 0.986495846869529",
  "bbob-biobj_f55_i07_d10 0.938154359689089",
  "bbob-biobj_f55_i07_d20 0.737202106400754",
  "bbob-biobj_f55_i07_d40 0.409115232136645",
  "bbob-biobj_f55_i08_d02 0.901039061074305",
  "bbob-biobj_f55_i08_d03 0.984416828859526",
  "bbob-biobj_f55_i08_d05 0.985698579395567",
  "bbob-biobj_f55_i08_d10 0.852092110104912",
  "bbob-biobj_f55_i08_d20 0.695590004902784",
  "bbob-biobj_f55_i08_d40 0.263145726215256",
  "bbob-biobj_f55_i09_d02 0.869205310820279",
  "bbob-biobj_f55_i09_d03 0.938417264913280",
  "bbob-biobj_f55_i09_d05 0.954969728866659",
  "bbob-biobj_f55_i09_d10 0.910621949104183",
  "bbob-biobj_f55_i09_d20 0.753155231009609",
  "bbob-biobj_f55_i09_d40 0.444119724681479",
  "bbob-biobj_f55_i10_d02 0.851268210876801",
  "bbob-biobj_f55_i10_d03 0.983908789480269",
  "bbob-biobj_f55_i10_d05 0.960192981087286",
  "bbob-biobj_f55_i10_d10 0.892927198168716",
  "bbob-biobj_f55_i10_d20 0.732377974441418",
  "bbob-biobj_f55_i10_d40 0.270212188947950"
};
#line 19 "code-experiments/src/suite_biobj.c"

/**
 * @brief The array of triples biobj_instance - problem1_instance - problem2_instance connecting bi-objective
 * suite instances to the instances of the bbob suite.
 *
 * It should be updated with new instances when they are chosen.
 */
static const size_t suite_biobj_instances[][3] = {
    { 1, 2, 4 },
    { 2, 3, 5 },
    { 3, 7, 8 },
    { 4, 9, 10 },
    { 5, 11, 12 },
    { 6, 13, 14 },
    { 7, 15, 16 },
    { 8, 17, 18 },
    { 9, 19, 21 },
    { 10, 21, 22 }
};

/**
 * @brief The bbob-biobj suite data type.
 */
typedef struct {

  size_t **new_instances;    /**< @brief A matrix of new instances (equal in form to suite_biobj_instances)
                                   that needs to be used only when an instance that is not in
                                   suite_biobj_instances is being invoked. */

  size_t max_new_instances;  /**< @brief The maximal number of new instances. */

} suite_biobj_t;

static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                         const size_t number_of_functions,
                                         const size_t number_of_dimensions,
                                         const size_t *dimensions,
                                         const char *default_instances);
static void suite_biobj_free(void *stuff);
static size_t suite_biobj_get_new_instance(coco_suite_t *suite,
                                           const size_t instance,
                                           const size_t instance1,
                                           const size_t num_bbob_functions,
                                           const size_t *bbob_functions);

/**
 * @brief Sets the dimensions and default instances for the bbob-biobj suite.
 */
static coco_suite_t *suite_biobj_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = { 2, 3, 5, 10, 20, 40 };

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = coco_suite_allocate("bbob-biobj", 55, 6, dimensions, "year: 2016");

  return suite;
}

/**
 * @brief Sets the instances associated with years for the bbob-biobj suite.
 */
static const char *suite_biobj_get_instances_by_year(const int year) {

  if (year == 2016) {
    return "1-10";
  }
  else {
    coco_error("suite_biobj_get_instances_by_year(): year %d not defined for suite_biobj", year);
    return NULL;
  }
}

/**
 * @brief Returns the problem from the bbob-biobj suite that corresponds to the given parameters.
 *
 * Creates the bi-objective problem by constructing it from two single-objective problems from the bbob
 * suite. If the invoked instance number is not in suite_biobj_instances, the function uses the following
 * formula to construct a new appropriate instance:
 *
 *   problem1_instance = 2 * biobj_instance + 1
 *
 *   problem2_instance = problem1_instance + 1
 *
 * If needed, problem2_instance is increased (see also the explanation of suite_biobj_get_new_instance).
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_t *suite_biobj_get_problem(coco_suite_t *suite,
                                               const size_t function_idx,
                                               const size_t dimension_idx,
                                               const size_t instance_idx) {

  const size_t num_bbob_functions = 10;
  /* Functions from the bbob suite that are used to construct the bi-objective problem. */
  const size_t bbob_functions[] = { 1, 2, 6, 8, 13, 14, 15, 17, 20, 21 };

  coco_problem_t *problem1, *problem2, *problem = NULL;
  size_t function1_idx, function2_idx;
  size_t instance1 = 0, instance2 = 0;

  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];

  suite_biobj_t *data = (suite_biobj_t *) suite->data;
  size_t i, j;
  const size_t num_existing_instances = sizeof(suite_biobj_instances) / sizeof(suite_biobj_instances[0]);
  int instance_found = 0;

  double *smallest_values_of_interest = coco_allocate_vector_with_value(dimension, -100);
  double *largest_values_of_interest = coco_allocate_vector_with_value(dimension, 100);

  /* A "magic" formula to compute the BBOB function index from the bi-objective function index */
  function1_idx = num_bbob_functions
      - coco_double_to_size_t(
          floor(-0.5 + sqrt(0.25 + 2.0 * (double) (suite->number_of_functions - function_idx - 1)))) - 1;
  function2_idx = function_idx - (function1_idx * num_bbob_functions) +
      (function1_idx * (function1_idx + 1)) / 2;

  /* First search for instance in suite_biobj_instances */
  for (i = 0; i < num_existing_instances; i++) {
    if (suite_biobj_instances[i][0] == instance) {
      /* The instance has been found in suite_biobj_instances */
      instance1 = suite_biobj_instances[i][1];
      instance2 = suite_biobj_instances[i][2];
      instance_found = 1;
      break;
    }
  }

  if ((!instance_found) && (data)) {
    /* Next, search for instance in new_instances */
    for (i = 0; i < data->max_new_instances; i++) {
      if (data->new_instances[i][0] == 0)
        break;
      if (data->new_instances[i][0] == instance) {
        /* The instance has been found in new_instances */
        instance1 = data->new_instances[i][1];
        instance2 = data->new_instances[i][2];
        instance_found = 1;
        break;
      }
    }
  }

  if (!instance_found) {
    /* Finally, if the instance is not found, create a new one */

    if (!data) {
      /* Allocate space needed for saving new instances */
      data = (suite_biobj_t *) coco_allocate_memory(sizeof(*data));

      /* Most often the actual number of new instances will be lower than max_new_instances, because
       * some of them are already in suite_biobj_instances. However, in order to avoid iterating over
       * suite_biobj_instances, the allocation uses max_new_instances. */
      data->max_new_instances = suite->number_of_instances;

      data->new_instances = (size_t **) coco_allocate_memory(data->max_new_instances * sizeof(size_t *));
      for (i = 0; i < data->max_new_instances; i++) {
        data->new_instances[i] = (size_t *) malloc(3 * sizeof(size_t));
        for (j = 0; j < 3; j++) {
          data->new_instances[i][j] = 0;
        }
      }
      suite->data_free_function = suite_biobj_free;
      suite->data = data;
    }

    /* A simple formula to set the first instance */
    instance1 = 2 * instance + 1;
    instance2 = suite_biobj_get_new_instance(suite, instance, instance1, num_bbob_functions, bbob_functions);
  }

  problem1 = coco_get_bbob_problem(bbob_functions[function1_idx], dimension, instance1);
  problem2 = coco_get_bbob_problem(bbob_functions[function2_idx], dimension, instance2);

  problem = coco_problem_stacked_allocate(problem1, problem2, smallest_values_of_interest, largest_values_of_interest);

  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(suite, function_idx, dimension_idx, instance_idx);

  /* Use the standard stacked problem_id as problem_name and construct a new suite-specific problem_id */
  coco_problem_set_name(problem, problem->problem_id);
  coco_problem_set_id(problem, "bbob-biobj_f%02lu_i%02lu_d%02lu", (unsigned long) function,
  		(unsigned long) instance, (unsigned long) dimension);

  /* Construct problem type */
  coco_problem_set_type(problem, "%s_%s", problem1->problem_type, problem2->problem_type);

  coco_free_memory(smallest_values_of_interest);
  coco_free_memory(largest_values_of_interest);

  return problem;
}

/**
 * @brief Computes the instance number of the second problem/objective so that the resulting bi-objective
 * problem has more than a single optimal solution.
 *
 * Starts by setting instance2 = instance1 + 1 and increases this number until an appropriate instance has
 * been found (or until a maximum number of tries has been reached, in which case it throws a coco_error).
 * An appropriate instance is the one for which the resulting bi-objective problem (in any considered
 * dimension) has the ideal and nadir points apart enough in the objective space and the extreme optimal
 * points apart enough in the decision space. When the instance has been found, it is output through
 * coco_warning, so that the user can see it and eventually manually add it to suite_biobj_instances.
 */
static size_t suite_biobj_get_new_instance(coco_suite_t *suite,
                                           const size_t instance,
                                           const size_t instance1,
                                           const size_t num_bbob_functions,
                                           const size_t *bbob_functions) {
  size_t instance2 = 0;
  size_t num_tries = 0;
  const size_t max_tries = 1000;
  const double apart_enough = 1e-4;
  int appropriate_instance_found = 0, break_search, warning_produced = 0;
  coco_problem_t *problem1, *problem2, *problem = NULL;
  size_t d, f1, f2, i;
  size_t function1, function2, dimension;
  double norm;
  double *smallest_values_of_interest, *largest_values_of_interest;

  suite_biobj_t *data;
  assert(suite->data);
  data = (suite_biobj_t *) suite->data;

  while ((!appropriate_instance_found) && (num_tries < max_tries)) {
    num_tries++;
    instance2 = instance1 + num_tries;
    break_search = 0;

    /* An instance is "appropriate" if the ideal and nadir points in the objective space and the two
     * extreme optimal points in the decisions space are apart enough for all problems (all dimensions
     * and function combinations); therefore iterate over all dimensions and function combinations  */
    for (f1 = 0; (f1 < num_bbob_functions) && !break_search; f1++) {
      function1 = bbob_functions[f1];
      for (f2 = f1; (f2 < num_bbob_functions) && !break_search; f2++) {
        function2 = bbob_functions[f2];
        for (d = 0; (d < suite->number_of_dimensions) && !break_search; d++) {
          dimension = suite->dimensions[d];

          if (dimension == 0) {
            if (!warning_produced)
              coco_warning("suite_biobj_get_new_instance(): remove filtering of dimensions to get generally acceptable instances!");
            warning_produced = 1;
            continue;
          }

          problem1 = coco_get_bbob_problem(function1, dimension, instance1);
          problem2 = coco_get_bbob_problem(function2, dimension, instance2);
          if (problem) {
            coco_problem_stacked_free(problem);
            problem = NULL;
          }

          /* Set smallest and largest values of interest to some value (not important which, it just needs to be a
           * vector of doubles of the right dimension) */
          smallest_values_of_interest = coco_allocate_vector_with_value(dimension, -100);
          largest_values_of_interest = coco_allocate_vector_with_value(dimension, 100);
          problem = coco_problem_stacked_allocate(problem1, problem2, smallest_values_of_interest,
          		largest_values_of_interest);
          coco_free_memory(smallest_values_of_interest);
          coco_free_memory(largest_values_of_interest);

          /* Check whether the ideal and nadir points are too close in the objective space */
          norm = mo_get_norm(problem->best_value, problem->nadir_value, 2);
          if (norm < 1e-1) { /* TODO How to set this value in a sensible manner? */
            coco_debug(
                "suite_biobj_get_new_instance(): The ideal and nadir points of %s are too close in the objective space",
                problem->problem_id);
            coco_debug("norm = %e, ideal = %e\t%e, nadir = %e\t%e", norm, problem->best_value[0],
                problem->best_value[1], problem->nadir_value[0], problem->nadir_value[1]);
            break_search = 1;
          }

          /* Check whether the extreme optimal points are too close in the decision space */
          norm = mo_get_norm(problem1->best_parameter, problem2->best_parameter, problem->number_of_variables);
          if (norm < apart_enough) {
            coco_debug(
                "suite_biobj_get_new_instance(): The extreme points of %s are too close in the decision space",
                problem->problem_id);
            coco_debug("norm = %e", norm);
            break_search = 1;
          }
        }
      }
    }
    /* Clean up */
    if (problem) {
      coco_problem_stacked_free(problem);
      problem = NULL;
    }

    if (break_search) {
      /* The search was broken, continue with next instance2 */
      continue;
    } else {
      /* An appropriate instance was found */
      appropriate_instance_found = 1;
      coco_info("suite_biobj_set_new_instance(): Instance %lu created from instances %lu and %lu",
      		(unsigned long) instance, (unsigned long) instance1, (unsigned long) instance2);

      /* Save the instance to new_instances */
      for (i = 0; i < data->max_new_instances; i++) {
        if (data->new_instances[i][0] == 0) {
          data->new_instances[i][0] = instance;
          data->new_instances[i][1] = instance1;
          data->new_instances[i][2] = instance2;
          break;
        };
      }
    }
  }

  if (!appropriate_instance_found) {
    coco_error("suite_biobj_get_new_instance(): Could not find suitable instance %lu in %lu tries",
    		(unsigned long) instance, (unsigned long) num_tries);
    return 0; /* Never reached */
  }

  return instance2;
}

/**
 * @brief  Frees the memory of the given bi-objective suite.
 */
static void suite_biobj_free(void *stuff) {

  suite_biobj_t *data;
  size_t i;

  assert(stuff != NULL);
  data = (suite_biobj_t *) stuff;

  if (data->new_instances) {
    for (i = 0; i < data->max_new_instances; i++) {
      if (data->new_instances[i]) {
        coco_free_memory(data->new_instances[i]);
        data->new_instances[i] = NULL;
      }
    }
  }
  coco_free_memory(data->new_instances);
  data->new_instances = NULL;
}

/**
 * @brief Returns the best known value for indicator_name matching the given key if the key is found, and
 * throws a coco_error otherwise.
 */
static double suite_biobj_get_best_value(const char *indicator_name, const char *key) {

  size_t i, count;
  double best_value = 0;
  char *curr_key;

  if (strcmp(indicator_name, "hyp") == 0) {

    curr_key = coco_allocate_string(COCO_PATH_MAX);
    count = sizeof(suite_biobj_best_values_hyp) / sizeof(char *);
    for (i = 0; i < count; i++) {
      sscanf(suite_biobj_best_values_hyp[i], "%s %lf", curr_key, &best_value);
      if (strcmp(curr_key, key) == 0) {
        coco_free_memory(curr_key);
        return best_value;
      }
    }

    coco_free_memory(curr_key);
    coco_warning("suite_biobj_get_best_value(): best value of %s could not be found; set to 1.0", key);
    return 1.0;

  } else {
    coco_error("suite_biobj_get_best_value(): indicator %s not supported", indicator_name);
    return 0; /* Never reached */
  }

  coco_error("suite_biobj_get_best_value(): unexpected exception");
  return 0; /* Never reached */
}
#line 20 "code-experiments/src/coco_suite.c"
#line 1 "code-experiments/src/suite_toy.c"
/**
 * @file suite_toy.c
 * @brief Implementation of a toy suite containing 6 noiseless "basic" single-objective functions in 5
 * dimensions.
 */

#line 8 "code-experiments/src/suite_toy.c"
#line 9 "code-experiments/src/suite_toy.c"
#line 10 "code-experiments/src/suite_toy.c"
#line 11 "code-experiments/src/suite_toy.c"
#line 12 "code-experiments/src/suite_toy.c"
#line 13 "code-experiments/src/suite_toy.c"
#line 14 "code-experiments/src/suite_toy.c"

static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                         const size_t number_of_functions,
                                         const size_t number_of_dimensions,
                                         const size_t *dimensions,
                                         const char *default_instances);

/**
 * @brief Sets the dimensions and default instances for the toy suite.
 */
static coco_suite_t *suite_toy_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = { 2, 3, 5, 10, 20 };

  suite = coco_suite_allocate("toy", 6, 3, dimensions, "instances:1");

  return suite;
}

/**
 * @brief Returns the problem from the toy suite that corresponds to the given parameters.
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_t *suite_toy_get_problem(coco_suite_t *suite,
                                             const size_t function_idx,
                                             const size_t dimension_idx,
                                             const size_t instance_idx) {


  coco_problem_t *problem = NULL;

  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];

  if (function == 1) {
    problem = f_sphere_allocate(dimension);
  } else if (function == 2) {
    problem = f_ellipsoid_allocate(dimension);
  } else if (function == 3) {
    problem = f_rastrigin_allocate(dimension);
  } else if (function == 4) {
    problem = f_bueche_rastrigin_allocate(dimension);
  } else if (function == 5) {
    double xopt[40] = { 5.0 };
    problem = f_linear_slope_allocate(dimension, xopt);
  } else if (function == 6) {
    problem = f_rosenbrock_allocate(dimension);
  } else {
    coco_error("suite_toy_get_problem(): function %lu does not exist in this suite", (unsigned long) function);
    return NULL; /* Never reached */
  }

  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(suite, function_idx, dimension_idx, instance_idx);

  return problem;
}
#line 21 "code-experiments/src/coco_suite.c"
#line 1 "code-experiments/src/suite_largescale.c"
/**
 * @file suite_largescale.c
 * @brief Implementation of the bbob large-scale suite containing 1 function in 6 large dimensions.
 */

#line 7 "code-experiments/src/suite_largescale.c"

#line 9 "code-experiments/src/suite_largescale.c"

static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                         const size_t number_of_functions,
                                         const size_t number_of_dimensions,
                                         const size_t *dimensions,
                                         const char *default_instances);

/**
 * @brief Sets the dimensions and default instances for the bbob large-scale suite.
 */
static coco_suite_t *suite_largescale_initialize(void) {
  
  coco_suite_t *suite;
  /*const size_t dimensions[] = { 8, 16, 32, 64, 128, 256,512,1024};*/
  const size_t dimensions[] = { 40, 80, 160, 320, 640, 1280};
  suite = coco_suite_allocate("bbob-largescale", 1, 6, dimensions, "instances:1-15");
  return suite;
}

/**
 * @brief Creates and returns a large-scale problem without needing the actual large-scale suite.
 */
static coco_problem_t *coco_get_largescale_problem(const size_t function,
                                                   const size_t dimension,
                                                   const size_t instance) {
  coco_problem_t *problem = NULL;

  const char *problem_id_template = "bbob_f%03lu_i%02lu_d%02lu";
  const char *problem_name_template = "BBOB suite problem f%lu instance %lu in %luD";

  const long rseed = (long) (function + 10000 * instance);
  /*const long rseed_3 = (long) (3 + 10000 * instance);*/
  /*const long rseed_17 = (long) (17 + 10000 * instance);*/
  if (function == 1) {
    problem = f_ellipsoid_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
        problem_id_template, problem_name_template);
  } else {
    coco_error("coco_get_largescale_problem(): cannot retrieve problem f%lu instance %lu in %luD",
    		(unsigned long) function, (unsigned long) instance, (unsigned long) dimension);
    return NULL; /* Never reached */
  }

  return problem;
}

/**
 * @brief Returns the problem from the bbob large-scale suite that corresponds to the given parameters.
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_t *suite_largescale_get_problem(coco_suite_t *suite,
                                                    const size_t function_idx,
                                                    const size_t dimension_idx,
                                                    const size_t instance_idx) {
  
  coco_problem_t *problem = NULL;
  
  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];
  
  problem = coco_get_largescale_problem(function, dimension, instance);
  
  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(suite, function_idx, dimension_idx, instance_idx);
  
  return problem;
}
#line 22 "code-experiments/src/coco_suite.c"

/** @brief The maximum number of different instances in a suite. */
#define COCO_MAX_INSTANCES 1000

/**
 * @brief Calls the initializer of the given suite.
 *
 * @note This function needs to be updated when a new suite is added to COCO.
 */
static coco_suite_t *coco_suite_intialize(const char *suite_name) {

  coco_suite_t *suite;

  if (strcmp(suite_name, "toy") == 0) {
    suite = suite_toy_initialize();
  } else if (strcmp(suite_name, "bbob") == 0) {
    suite = suite_bbob_initialize();
  } else if (strcmp(suite_name, "bbob-biobj") == 0) {
    suite = suite_biobj_initialize();
  } else if (strcmp(suite_name, "bbob-largescale") == 0) {
    suite = suite_largescale_initialize();
  }
  else {
    coco_error("coco_suite_intialize(): unknown problem suite");
    return NULL;
  }

  return suite;
}

/**
 * @brief Calls the function that sets the instanced by year for the given suite.
 *
 * @note This function needs to be updated when a new suite is added to COCO.
 */
static const char *coco_suite_get_instances_by_year(const coco_suite_t *suite, const int year) {
  const char *year_string;

  if (strcmp(suite->suite_name, "bbob") == 0) {
    year_string = suite_bbob_get_instances_by_year(year);
  } else if (strcmp(suite->suite_name, "bbob-biobj") == 0) {
    year_string = suite_biobj_get_instances_by_year(year);
  } else {
    coco_error("coco_suite_get_instances_by_year(): suite '%s' has no years defined", suite->suite_name);
    return NULL;
  }

  return year_string;
}

/**
 * @brief Calls the function that returns the problem corresponding to the given suite, function index,
 * dimension index and instance index. If the indices don't correspond to a problem because of suite
 * filtering, it returns NULL.
 *
 * @note This function needs to be updated when a new suite is added to COCO.
 */
static coco_problem_t *coco_suite_get_problem_from_indices(coco_suite_t *suite,
                                                           const size_t function_idx,
                                                           const size_t dimension_idx,
                                                           const size_t instance_idx) {

  coco_problem_t *problem;

  if ((suite->functions[function_idx] == 0) ||
      (suite->dimensions[dimension_idx] == 0) ||
	  (suite->instances[instance_idx] == 0)) {
	  return NULL;
  }

  if (strcmp(suite->suite_name, "toy") == 0) {
    problem = suite_toy_get_problem(suite, function_idx, dimension_idx, instance_idx);
  } else if (strcmp(suite->suite_name, "bbob") == 0) {
    problem = suite_bbob_get_problem(suite, function_idx, dimension_idx, instance_idx);
  } else if (strcmp(suite->suite_name, "bbob-biobj") == 0) {
    problem = suite_biobj_get_problem(suite, function_idx, dimension_idx, instance_idx);
  } else if (strcmp(suite->suite_name, "bbob-largescale") == 0) {
    problem = suite_largescale_get_problem(suite, function_idx, dimension_idx, instance_idx);
  } else {
    coco_error("coco_suite_get_problem_from_indices(): unknown problem suite");
    return NULL;
  }

  coco_problem_set_suite(problem, suite);

  return problem;
}

/**
 * @note: While a suite can contain multiple problems with equal function, dimension and instance, this
 * function always returns the first problem in the suite with the given function, dimension and instance
 * values. If the given values don't correspond to a problem, the function returns NULL.
 */
coco_problem_t *coco_suite_get_problem_by_function_dimension_instance(coco_suite_t *suite,
                                                                      const size_t function,
                                                                      const size_t dimension,
                                                                      const size_t instance) {

  size_t i;
  int function_idx, dimension_idx, instance_idx;
  int found;

  found = 0;
  for (i = 0; i < suite->number_of_functions; i++) {
    if (suite->functions[i] == function) {
      function_idx = (int) i;
      found = 1;
      break;
    }
  }
  if (!found)
    return NULL;

  found = 0;
  for (i = 0; i < suite->number_of_dimensions; i++) {
    if (suite->dimensions[i] == dimension) {
      dimension_idx = (int) i;
      found = 1;
      break;
    }
  }
  if (!found)
    return NULL;

  found = 0;
  for (i = 0; i < suite->number_of_instances; i++) {
    if (suite->instances[i] == instance) {
      instance_idx = (int) i;
      found = 1;
      break;
    }
  }
  if (!found)
    return NULL;

  return coco_suite_get_problem_from_indices(suite, (size_t) function_idx, (size_t) dimension_idx, (size_t) instance_idx);
}


/**
 * @brief Allocates the space for a coco_suite_t instance.
 *
 * This function sets the functions and dimensions contained in the suite, while the instances are set by
 * the function coco_suite_set_instance.
 */
static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                         const size_t number_of_functions,
                                         const size_t number_of_dimensions,
                                         const size_t *dimensions,
                                         const char *default_instances) {

  coco_suite_t *suite;
  size_t i;

  suite = (coco_suite_t *) coco_allocate_memory(sizeof(*suite));

  suite->suite_name = coco_strdup(suite_name);

  suite->number_of_dimensions = number_of_dimensions;
  assert(number_of_dimensions > 0);
  suite->dimensions = coco_allocate_vector_size_t(suite->number_of_dimensions);
  for (i = 0; i < suite->number_of_dimensions; i++) {
    suite->dimensions[i] = dimensions[i];
  }

  suite->number_of_functions = number_of_functions;
  assert(number_of_functions > 0);
  suite->functions = coco_allocate_vector_size_t(suite->number_of_functions);
  for (i = 0; i < suite->number_of_functions; i++) {
    suite->functions[i] = i + 1;
  }

  assert(strlen(default_instances) > 0);
  suite->default_instances = coco_strdup(default_instances);

  /* Will be set to the first valid dimension index before the constructor ends */
  suite->current_dimension_idx = -1;
  /* Will be set to the first valid function index before the constructor ends  */
  suite->current_function_idx = -1;

  suite->current_instance_idx = -1;
  suite->current_problem = NULL;

  /* To be set in coco_suite_set_instance() */
  suite->number_of_instances = 0;
  suite->instances = NULL;

  /* To be set in particular suites if needed */
  suite->data = NULL;
  suite->data_free_function = NULL;

  return suite;
}

/**
 * @brief Sets the suite instance to the given instance_numbers.
 */
static void coco_suite_set_instance(coco_suite_t *suite,
                                    const size_t *instance_numbers) {

  size_t i;

  if (!instance_numbers) {
    coco_error("coco_suite_set_instance(): no instance given");
    return;
  }

  suite->number_of_instances = coco_count_numbers(instance_numbers, COCO_MAX_INSTANCES, "suite instance numbers");
  suite->instances = coco_allocate_vector_size_t(suite->number_of_instances);
  for (i = 0; i < suite->number_of_instances; i++) {
    suite->instances[i] = instance_numbers[i];
  }

}

/**
 * @brief Filters the given items w.r.t. the given indices (starting from 1).
 *
 * Sets items[i - 1] to 0 for every i that cannot be found in indices (this function performs the conversion
 * from user-friendly indices starting from 1 to C-friendly indices starting from 0).
 */
static void coco_suite_filter_indices(size_t *items, const size_t number_of_items, const size_t *indices, const char *name) {

  size_t i, j;
  size_t count = coco_count_numbers(indices, COCO_MAX_INSTANCES, name);
  int found;

  for (i = 1; i <= number_of_items; i++) {
    found = 0;
    for (j = 0; j < count; j++) {
      if (i == indices[j]) {
        found = 1;
        break;
      }
    }
    if (!found)
      items[i - 1] = 0;
  }

}

/**
 * @brief Filters dimensions w.r.t. the given dimension_numbers.
 *
 * Sets suite->dimensions[i] to 0 for every dimension value that cannot be found in dimension_numbers.
 */
static void coco_suite_filter_dimensions(coco_suite_t *suite, const size_t *dimension_numbers) {

  size_t i, j;
  size_t count = coco_count_numbers(dimension_numbers, COCO_MAX_INSTANCES, "dimensions");
  int found;

  for (i = 0; i < suite->number_of_dimensions; i++) {
    found = 0;
    for (j = 0; j < count; j++) {
      if (suite->dimensions[i] == dimension_numbers[j])
        found = 1;
    }
    if (!found)
      suite->dimensions[i] = 0;
  }

}

/**
 * @param suite The given suite.
 * @param function_idx The index of the function in question (starting from 0).
 *
 * @return The function number in position function_idx in the suite. If the function has been filtered out
 * through suite_options in the coco_suite function, the result is 0.
 */
size_t coco_suite_get_function_from_function_index(const coco_suite_t *suite, const size_t function_idx) {

  if (function_idx >= suite->number_of_functions) {
    coco_error("coco_suite_get_function_from_function_index(): function index exceeding the number of functions in the suite");
    return 0; /* Never reached*/
  }

 return suite->functions[function_idx];
}

/**
 * @param suite The given suite.
 * @param dimension_idx The index of the dimension in question (starting from 0).
 *
 * @return The dimension number in position dimension_idx in the suite. If the dimension has been filtered out
 * through suite_options in the coco_suite function, the result is 0.
 */
size_t coco_suite_get_dimension_from_dimension_index(const coco_suite_t *suite, const size_t dimension_idx) {

  if (dimension_idx >= suite->number_of_dimensions) {
    coco_error("coco_suite_get_dimension_from_dimension_index(): dimensions index exceeding the number of dimensions in the suite");
    return 0; /* Never reached*/
  }

 return suite->dimensions[dimension_idx];
}

/**
 * @param suite The given suite.
 * @param instance_idx The index of the instance in question (starting from 0).
 *
 * @return The instance number in position instance_idx in the suite. If the instance has been filtered out
 * through suite_options in the coco_suite function, the result is 0.
 */
size_t coco_suite_get_instance_from_instance_index(const coco_suite_t *suite, const size_t instance_idx) {

  if (instance_idx >= suite->number_of_instances) {
    coco_error("coco_suite_get_instance_from_instance_index(): instance index exceeding the number of instances in the suite");
    return 0; /* Never reached*/
  }

 return suite->functions[instance_idx];
}

void coco_suite_free(coco_suite_t *suite) {

  if (suite != NULL) {

    if (suite->suite_name) {
      coco_free_memory(suite->suite_name);
      suite->suite_name = NULL;
    }
    if (suite->dimensions) {
      coco_free_memory(suite->dimensions);
      suite->dimensions = NULL;
    }
    if (suite->functions) {
      coco_free_memory(suite->functions);
      suite->functions = NULL;
    }
    if (suite->instances) {
      coco_free_memory(suite->instances);
      suite->instances = NULL;
    }
    if (suite->default_instances) {
      coco_free_memory(suite->default_instances);
      suite->default_instances = NULL;
    }

    if (suite->current_problem) {
      coco_problem_free(suite->current_problem);
      suite->current_problem = NULL;
    }

    if (suite->data != NULL) {
      if (suite->data_free_function != NULL) {
        suite->data_free_function(suite->data);
      }
      coco_free_memory(suite->data);
      suite->data = NULL;
    }

    coco_free_memory(suite);
    suite = NULL;
  }
}

/**
 * Note that the problem_index depends on the number of instances a suite is defined with.
 *
 * @param suite The given suite.
 * @param problem_index The index of the problem to be returned.
 *
 * @return The problem of the suite defined by problem_index (NULL if this problem has been filtered out
 * from the suite).
 */
coco_problem_t *coco_suite_get_problem(coco_suite_t *suite, const size_t problem_index) {

  size_t function_idx = 0, instance_idx = 0, dimension_idx = 0;
  coco_suite_decode_problem_index(suite, problem_index, &function_idx, &dimension_idx, &instance_idx);

  return coco_suite_get_problem_from_indices(suite, function_idx, dimension_idx, instance_idx);
}

/**
 * The number of problems in the suite is computed as a product of the number of instances, number of
 * functions and number of dimensions and therefore doesn't account for any filtering done through the
 * suite_options parameter of the coco_suite function.
 *
 * @param suite The given suite.
 *
 * @return The number of problems in the suite.
 */
size_t coco_suite_get_number_of_problems(const coco_suite_t *suite) {
  return (suite->number_of_instances * suite->number_of_functions * suite->number_of_dimensions);
}


/**
 * @brief Returns the instances read from either a "year: YEAR" or "instances: NUMBERS" string.
 *
 * If both "year" and "instances" are given, the second is ignored (and a warning is raised). See the
 * coco_suite function for more information about the required format.
 */
static size_t *coco_suite_get_instance_indices(const coco_suite_t *suite, const char *suite_instance) {

  int year = -1;
  char *instances = NULL;
  const char *year_string = NULL;
  long year_found, instances_found;
  int parce_year = 1, parce_instances = 1;
  size_t *result = NULL;

  if (suite_instance == NULL)
    return NULL;

  year_found = coco_strfind(suite_instance, "year");
  instances_found = coco_strfind(suite_instance, "instances");

  if ((year_found < 0) && (instances_found < 0))
    return NULL;

  if ((year_found > 0) && (instances_found > 0)) {
    if (year_found < instances_found) {
      parce_instances = 0;
      coco_warning("coco_suite_get_instance_indices(): 'instances' suite option ignored because it follows 'year'");
    }
    else {
      parce_year = 0;
      coco_warning("coco_suite_get_instance_indices(): 'year' suite option ignored because it follows 'instances'");
    }
  }

  if ((year_found >= 0) && (parce_year == 1)) {
    if (coco_options_read_int(suite_instance, "year", &(year)) != 0) {
      year_string = coco_suite_get_instances_by_year(suite, year);
      result = coco_string_parse_ranges(year_string, 1, 0, "instances", COCO_MAX_INSTANCES);
    } else {
      coco_warning("coco_suite_get_instance_indices(): problems parsing the 'year' suite_instance option, ignored");
    }
  }

  instances = coco_allocate_string(COCO_MAX_INSTANCES);
  if ((instances_found >= 0) && (parce_instances == 1)) {
    if (coco_options_read_values(suite_instance, "instances", instances) > 0) {
      result = coco_string_parse_ranges(instances, 1, 0, "instances", COCO_MAX_INSTANCES);
    } else {
      coco_warning("coco_suite_get_instance_indices(): problems parsing the 'instance' suite_instance option, ignored");
    }
  }
  coco_free_memory(instances);

  return result;
}

/**
 * @brief Iterates through the items from the current_item_id position on in search for the next positive
 * item.
 *
 * If such an item is found, current_item_id points to this item and the method returns 1. If such an
 * item cannot be found, current_item_id points to the first positive item and the method returns 0.
 */
static int coco_suite_is_next_item_found(const size_t *items, const size_t number_of_items, long *current_item_id) {

  if ((*current_item_id) != number_of_items - 1)  {
    /* Not the last item, iterate through items */
    do {
      (*current_item_id)++;
    } while (((*current_item_id) < number_of_items - 1) && (items[*current_item_id] == 0));

    assert((*current_item_id) < number_of_items);
    if (items[*current_item_id] != 0) {
      /* Next item is found, return true */
      return 1;
    }
  }

  /* Next item cannot be found, move to the first good item and return false */
  *current_item_id = -1;
  do {
    (*current_item_id)++;
  } while ((*current_item_id < number_of_items - 1) && (items[*current_item_id] == 0));
  if (items[*current_item_id] == 0)
    coco_error("coco_suite_is_next_item_found(): the chosen suite has no valid (positive) items");
  return 0;
}

/**
 * @brief Iterates through the instances of the given suite from the current_instance_idx position on in
 * search for the next positive instance.
 *
 * If such an instance is found, current_instance_idx points to this instance and the method returns 1. If
 * such an instance cannot be found, current_instance_idx points to the first positive instance and the
 * method returns 0.
 */
static int coco_suite_is_next_instance_found(coco_suite_t *suite) {

  return coco_suite_is_next_item_found(suite->instances, suite->number_of_instances,
      &suite->current_instance_idx);
}

/**
 * @brief Iterates through the functions of the given suite from the current_function_idx position on in
 * search for the next positive function.
 *
 * If such a function is found, current_function_idx points to this function and the method returns 1. If
 * such a function cannot be found, current_function_idx points to the first positive function,
 * current_instance_idx points to the first positive instance and the method returns 0.
 */
static int coco_suite_is_next_function_found(coco_suite_t *suite) {

  int result = coco_suite_is_next_item_found(suite->functions, suite->number_of_functions,
      &suite->current_function_idx);
  if (!result) {
    /* Reset the instances */
    suite->current_instance_idx = -1;
    coco_suite_is_next_instance_found(suite);
  }
  return result;
}

/**
 * @brief Iterates through the dimensions of the given suite from the current_dimension_idx position on in
 * search for the next positive dimension.
 *
 * If such a dimension is found, current_dimension_idx points to this dimension and the method returns 1. If
 * such a dimension cannot be found, current_dimension_idx points to the first positive dimension and the
 * method returns 0.
 */
static int coco_suite_is_next_dimension_found(coco_suite_t *suite) {

  return coco_suite_is_next_item_found(suite->dimensions, suite->number_of_dimensions,
      &suite->current_dimension_idx);
}

/**
 * Currently, four suites are supported:
 * - "bbob" contains 24 <a href="http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf">
 * single-objective functions</a> in 6 dimensions (2, 3, 5, 10, 20, 40)
 * - "bbob-biobj" contains 55 <a href="http://numbbo.github.io/coco-doc/bbob-biobj/functions">bi-objective
 * functions</a> in 6 dimensions (2, 3, 5, 10, 20, 40)
 * - "bbob-largescale" contains 24 <a href="http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf">
 * single-objective functions</a> in 6 large dimensions (40, 80, 160, 320, 640, 1280)
 * - "toy" contains 6 <a href="http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf">
 * single-objective functions</a> in 5 dimensions (2, 3, 5, 10, 20)
 *
 * Only the suite_name parameter needs to be non-empty. The suite_instance and suite_options can be "" or
 * NULL. In this case, default values are taken (default instances of a suite are those used in the last year
 * and the suite is not filtered by default).
 *
 * @param suite_name A string containing the name of the suite. Currently supported suite names are "bbob",
 * "bbob-biobj", "bbob-largescale" and "toy".
 * @param suite_instance A string used for defining the suite instances. Two ways are supported:
 * - "year: YEAR", where YEAR is the year of the BBOB workshop, includes the instances (to be) used in that
 * year's workshop;
 * - "instances: VALUES", where VALUES are instance numbers from 1 on written as a comma-separated list or a
 * range m-n.
 * @param suite_options A string of pairs "key: value" used to filter the suite (especially useful for
 * parallelizing the experiments). Supported options:
 * - "dimensions: LIST", where LIST is the list of dimensions to keep in the suite (range-style syntax is
 * not allowed here),
 * - "dimension_indices: VALUES", where VALUES is a list or a range of dimension indices (starting from 1) to keep
 * in the suite, and
 * - "function_indices: VALUES", where VALUES is a list or a range of function indices (starting from 1) to keep
 * in the suite, and
 * - "instance_indices: VALUES", where VALUES is a list or a range of instance indices (starting from 1) to keep
 * in the suite.
 *
 * @return The constructed suite object.
 */
coco_suite_t *coco_suite(const char *suite_name, const char *suite_instance, const char *suite_options) {

  coco_suite_t *suite;
  size_t *instances;
  char *option_string = NULL;
  char *ptr;
  size_t *indices = NULL;
  size_t *dimensions = NULL;
  long dim_found, dim_idx_found;
  int parce_dim = 1, parce_dim_idx = 1;

  coco_option_keys_t *known_option_keys, *given_option_keys, *redundant_option_keys;

  /* Sets the valid keys for suite options and suite instance */
  const char *known_keys_o[] = { "dimensions", "dimension_indices", "function_indices", "instance_indices" };
  const char *known_keys_i[] = { "year", "instances" };

  /* Initialize the suite */
  suite = coco_suite_intialize(suite_name);

  /* Set the instance */
  if ((!suite_instance) || (strlen(suite_instance) == 0))
    instances = coco_suite_get_instance_indices(suite, suite->default_instances);
  else {
    instances = coco_suite_get_instance_indices(suite, suite_instance);

    if (!instances) {
      /* Something wrong in the suite_instance string, use default instead */
      instances = coco_suite_get_instance_indices(suite, suite->default_instances);
    }

    /* Check for redundant option keys for suite instance */
    known_option_keys = coco_option_keys_allocate(sizeof(known_keys_i) / sizeof(char *), known_keys_i);
    given_option_keys = coco_option_keys(suite_instance);

    if (given_option_keys) {
      redundant_option_keys = coco_option_keys_get_redundant(known_option_keys, given_option_keys);

      if ((redundant_option_keys != NULL) && (redundant_option_keys->count > 0)) {
        /* Warn the user that some of given options are being ignored and output the valid options */
        char *output_redundant = coco_option_keys_get_output_string(redundant_option_keys,
            "coco_suite(): Some keys in suite instance were ignored:\n");
        char *output_valid = coco_option_keys_get_output_string(known_option_keys,
            "Valid keys for suite instance are:\n");
        coco_warning("%s%s", output_redundant, output_valid);
        coco_free_memory(output_redundant);
        coco_free_memory(output_valid);
      }

      coco_option_keys_free(given_option_keys);
      coco_option_keys_free(redundant_option_keys);
    }
    coco_option_keys_free(known_option_keys);
  }
  coco_suite_set_instance(suite, instances);
  coco_free_memory(instances);

  /* Apply filter if any given by the suite_options */
  if ((suite_options) && (strlen(suite_options) > 0)) {
    option_string = coco_allocate_string(COCO_PATH_MAX);
    if (coco_options_read_values(suite_options, "function_indices", option_string) > 0) {
      indices = coco_string_parse_ranges(option_string, 1, suite->number_of_functions, "function_indices", COCO_MAX_INSTANCES);
      if (indices != NULL) {
        coco_suite_filter_indices(suite->functions, suite->number_of_functions, indices, "function_indices");
        coco_free_memory(indices);
      }
    }
    coco_free_memory(option_string);

    option_string = coco_allocate_string(COCO_PATH_MAX);
    if (coco_options_read_values(suite_options, "instance_indices", option_string) > 0) {
      indices = coco_string_parse_ranges(option_string, 1, suite->number_of_instances, "instance_indices", COCO_MAX_INSTANCES);
      if (indices != NULL) {
        coco_suite_filter_indices(suite->instances, suite->number_of_instances, indices, "instance_indices");
        coco_free_memory(indices);
      }
    }
    coco_free_memory(option_string);

    dim_found = coco_strfind(suite_options, "dimensions");
    dim_idx_found = coco_strfind(suite_options, "dimension_indices");

    if ((dim_found > 0) && (dim_idx_found > 0)) {
      if (dim_found < dim_idx_found) {
        parce_dim_idx = 0;
        coco_warning("coco_suite(): 'dimension_indices' suite option ignored because it follows 'dimensions'");
      }
      else {
        parce_dim = 0;
        coco_warning("coco_suite(): 'dimensions' suite option ignored because it follows 'dimension_indices'");
      }
    }

    option_string = coco_allocate_string(COCO_PATH_MAX);
    if ((dim_idx_found >= 0) && (parce_dim_idx == 1)
        && (coco_options_read_values(suite_options, "dimension_indices", option_string) > 0)) {
      indices = coco_string_parse_ranges(option_string, 1, suite->number_of_dimensions, "dimension_indices",
          COCO_MAX_INSTANCES);
      if (indices != NULL) {
        coco_suite_filter_indices(suite->dimensions, suite->number_of_dimensions, indices, "dimension_indices");
        coco_free_memory(indices);
      }
    }
    coco_free_memory(option_string);

    option_string = coco_allocate_string(COCO_PATH_MAX);
    if ((dim_found >= 0) && (parce_dim == 1)
        && (coco_options_read_values(suite_options, "dimensions", option_string) > 0)) {
      ptr = option_string;
      /* Check for disallowed characters */
      while (*ptr != '\0') {
        if ((*ptr != ',') && !isdigit((unsigned char )*ptr)) {
          coco_warning("coco_suite(): 'dimensions' suite option ignored because of disallowed characters");
          return NULL;
        } else
          ptr++;
      }
      dimensions = coco_string_parse_ranges(option_string, suite->dimensions[0],
          suite->dimensions[suite->number_of_dimensions - 1], "dimensions", COCO_MAX_INSTANCES);
      if (dimensions != NULL) {
        coco_suite_filter_dimensions(suite, dimensions);
        coco_free_memory(dimensions);
      }
    }
    coco_free_memory(option_string);

    /* Check for redundant option keys for suite options */
    known_option_keys = coco_option_keys_allocate(sizeof(known_keys_o) / sizeof(char *), known_keys_o);
    given_option_keys = coco_option_keys(suite_options);

    if (given_option_keys) {
      redundant_option_keys = coco_option_keys_get_redundant(known_option_keys, given_option_keys);

      if ((redundant_option_keys != NULL) && (redundant_option_keys->count > 0)) {
        /* Warn the user that some of given options are being ignored and output the valid options */
        char *output_redundant = coco_option_keys_get_output_string(redundant_option_keys,
            "coco_suite(): Some keys in suite options were ignored:\n");
        char *output_valid = coco_option_keys_get_output_string(known_option_keys,
            "Valid keys for suite options are:\n");
        coco_warning("%s%s", output_redundant, output_valid);
        coco_free_memory(output_redundant);
        coco_free_memory(output_valid);
      }

      coco_option_keys_free(given_option_keys);
      coco_option_keys_free(redundant_option_keys);
    }
    coco_option_keys_free(known_option_keys);

  }

  /* Check that there are enough dimensions, functions and instances left */
  if ((suite->number_of_dimensions < 1)
      || (suite->number_of_functions < 1)
      || (suite->number_of_instances < 1)) {
    coco_error("coco_suite(): the suite does not contain at least one dimension, function and instance");
    return NULL;
  }

  /* Set the starting values of the current indices in such a way, that when the instance_idx is incremented,
   * this results in a valid problem */
  coco_suite_is_next_function_found(suite);
  coco_suite_is_next_dimension_found(suite);

  return suite;
}

/**
 * Iterates through the suite first by instances, then by functions and finally by dimensions.
 * The instances/functions/dimensions that have been filtered out using the suite_options of the coco_suite
 * function are skipped. Outputs some information regarding the current place in the iteration. The returned
 * problem is wrapped with the observer. If the observer is NULL, the returned problem is unobserved.
 *
 * @param suite The given suite.
 * @param observer The observer used to wrap the problem. If NULL, the problem is returned unobserved.
 *
 * @returns The next problem of the suite or NULL if there is no next problem left.
 */
coco_problem_t *coco_suite_get_next_problem(coco_suite_t *suite, coco_observer_t *observer) {

  size_t function_idx;
  size_t dimension_idx;
  size_t instance_idx;
  coco_problem_t *problem;

  long previous_function_idx;
  long previous_dimension_idx;
  long previous_instance_idx;

  assert(suite != NULL);

  previous_function_idx = suite->current_function_idx;
  previous_dimension_idx = suite->current_dimension_idx;
  previous_instance_idx = suite->current_instance_idx;

  /* Iterate through the suite by instances, then functions and lastly dimensions in search for the next
   * problem. Note that these functions set the values of suite fields current_instance_idx,
   * current_function_idx and current_dimension_idx. */
  if (!coco_suite_is_next_instance_found(suite)
      && !coco_suite_is_next_function_found(suite)
      && !coco_suite_is_next_dimension_found(suite)) {
    coco_info_partial("done\n");
    return NULL;
  }

  if (suite->current_problem) {
    coco_problem_free(suite->current_problem);
  }

  assert(suite->current_function_idx >= 0);
  assert(suite->current_dimension_idx >= 0);
  assert(suite->current_instance_idx >= 0);

  function_idx = (size_t) suite->current_function_idx;
  dimension_idx = (size_t) suite->current_dimension_idx;
  instance_idx = (size_t) suite->current_instance_idx;

  problem = coco_suite_get_problem_from_indices(suite, function_idx, dimension_idx, instance_idx);
  if (observer != NULL)
    problem = coco_problem_add_observer(problem, observer);
  suite->current_problem = problem;

  /* Output information regarding the current place in the iteration */
  if (coco_log_level >= COCO_INFO) {
    if (((long) dimension_idx != previous_dimension_idx) || (previous_instance_idx < 0)) {
      /* A new dimension started */
      char *time_string = coco_current_time_get_string();
      if (dimension_idx > 0)
        coco_info_partial("done\n");
      else
        coco_info_partial("\n");
      coco_info_partial("COCO INFO: %s, d=%lu, running: f%02lu", time_string,
      		(unsigned long) suite->dimensions[dimension_idx], (unsigned long) suite->functions[function_idx]);
      coco_free_memory(time_string);
    }
    else if ((long) function_idx != previous_function_idx){
      /* A new function started */
      coco_info_partial("f%02lu", (unsigned long) suite->functions[function_idx]);
    }
    /* One dot for each instance */
    coco_info_partial(".", suite->instances[instance_idx]);
  }

  return problem;
}

/* See coco.h for more information on encoding and decoding problem index */

/**
 * @param suite The suite.
 * @param function_idx Index of the function (starting with 0).
 * @param dimension_idx Index of the dimension (starting with 0).
 * @param instance_idx Index of the insatnce (starting with 0).
 *
 * @return The problem index in the suite computed from function_idx, dimension_idx and instance_idx.
 */
size_t coco_suite_encode_problem_index(const coco_suite_t *suite,
                                       const size_t function_idx,
                                       const size_t dimension_idx,
                                       const size_t instance_idx) {

  return instance_idx + (function_idx * suite->number_of_instances) +
      (dimension_idx * suite->number_of_instances * suite->number_of_functions);

}

/**
 * @param suite The suite.
 * @param problem_index Index of the problem in the suite (starting with 0).
 * @param function_idx Pointer to the index of the function, which is set by this function.
 * @param dimension_idx Pointer to the index of the dimension, which is set by this function.
 * @param instance_idx Pointer to the index of the instance, which is set by this function.
 */
void coco_suite_decode_problem_index(const coco_suite_t *suite,
                                     const size_t problem_index,
                                     size_t *function_idx,
                                     size_t *dimension_idx,
                                     size_t *instance_idx) {

  if (problem_index > (suite->number_of_instances * suite->number_of_functions * suite->number_of_dimensions) - 1) {
    coco_warning("coco_suite_decode_problem_index(): problem_index too large");
    function_idx = 0;
    instance_idx = 0;
    dimension_idx = 0;
    return;
  }

  *instance_idx = problem_index % suite->number_of_instances;
  *function_idx = (problem_index / suite->number_of_instances) % suite->number_of_functions;
  *dimension_idx = problem_index / (suite->number_of_instances * suite->number_of_functions);

}
#line 1 "code-experiments/src/coco_observer.c"
/**
 * @file coco_observer.c
 * @brief Definitions of functions regarding COCO observers.
 */

#line 7 "code-experiments/src/coco_observer.c"
#line 8 "code-experiments/src/coco_observer.c"
#include <limits.h>
#include <float.h>
#include <math.h>

/**
 * @brief The type for triggers based on target values.
 *
 * The target values that trigger logging are at every 10**(exponent/number_of_triggers) from positive
 * infinity down to precision, at 0, and from -precision on with step -10**(exponent/number_of_triggers) until
 * negative infinity.
 */
typedef struct {

  int exponent;               /**< @brief Value used to compare with the previously hit target. */
  double value;               /**< @brief Value of the currently hit target. */
  size_t number_of_triggers;  /**< @brief Number of target triggers between 10**i and 10**(i+1) for any i. */
  double precision;           /**< @brief Minimal precision of interest. */

} coco_observer_targets_t;

/**
 * @brief The type for triggers based on numbers of evaluations.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
typedef struct {

  /* First trigger */
  size_t value1;              /**< @brief The next value for the first trigger. */
  size_t exponent1;           /**< @brief Exponent used to compute the first trigger. */
  size_t number_of_triggers;  /**< @brief Number of target triggers between 10**i and 10**(i+1) for any i. */

  /* Second trigger */
  size_t value2;              /**< @brief The next value for the second trigger. */
  size_t exponent2;           /**< @brief Exponent used to compute the second trigger. */
  size_t *base_evaluations;   /**< @brief The base evaluation numbers used to compute the actual evaluation
                                   numbers that trigger logging. */
  size_t base_count;          /**< @brief The number of base evaluations. */
  size_t base_index;          /**< @brief The next index of the base evaluations. */
  size_t dimension;           /**< @brief Dimension used in the calculation of the first trigger. */

} coco_observer_evaluations_t;

/**
 * @brief The maximum number of evaluations to trigger logging.
 *
 * @note This is not the maximal evaluation number to be logged, but the maximal number of times logging is
 * triggered by the number of evaluations.
 */
#define COCO_MAX_EVALS_TO_LOG 1000

/***********************************************************************************************************/

/**
 * @name Methods regarding triggers based on target values
 */
/**@{*/

/**
 * @brief Creates and returns a structure containing information on targets.
 *
 * @param number_of_targets The number of targets between 10**(i/n) and 10**((i+1)/n) for each i.
 * @param precision Minimal precision of interest.
 */
static coco_observer_targets_t *coco_observer_targets(const size_t number_of_targets,
                                                      const double precision) {

  coco_observer_targets_t *targets = (coco_observer_targets_t *) coco_allocate_memory(sizeof(*targets));
  targets->exponent = INT_MAX;
  targets->value = DBL_MAX;
  targets->number_of_triggers = number_of_targets;
  targets->precision = precision;

  return targets;
}

/**
 * @brief Computes and returns whether the given value should trigger logging.
 */
static int coco_observer_targets_trigger(coco_observer_targets_t *targets, const double given_value) {

  int update_performed = 0;

  const double number_of_targets_double = (double) (long) targets->number_of_triggers;

  double verified_value = 0;
  int last_exponent = 0;
  int current_exponent = 0;
  int adjusted_exponent = 0;

  assert(targets != NULL);

  /* The given_value is positive or zero */
  if (given_value >= 0) {

  	if (given_value == 0) {
  		/* If zero, use even smaller value than precision */
  		verified_value = targets->precision / 10.0;
  	} else if (given_value < targets->precision) {
      /* If close to zero, use precision instead of the given_value*/
      verified_value = targets->precision;
    } else {
      verified_value = given_value;
    }

    current_exponent = (int) (ceil(log10(verified_value) * number_of_targets_double));

    /* If this is the first time the update was called, set the last_exponent to some value greater than the
     * current exponent */
    if (last_exponent == INT_MAX) {
      last_exponent = current_exponent + 1;
    } else {
      last_exponent = targets->exponent;
    }

    if (current_exponent < last_exponent) {
      /* Update the target information */
      targets->exponent = current_exponent;
      if (given_value == 0)
      	targets->value = 0;
      else
      	targets->value = pow(10, (double) current_exponent / number_of_targets_double);
      update_performed = 1;
    }
  }
  /* The given_value is negative, therefore adjustments need to be made */
  else {

    /* If close to zero, use precision instead of the given_value*/
    if (given_value > -targets->precision) {
      verified_value = targets->precision;
    } else {
      verified_value = -given_value;
    }

    /* Adjustment: use floor instead of ceil! */
    current_exponent = (int) (floor(log10(verified_value) * number_of_targets_double));

    /* If this is the first time the update was called, set the last_exponent to some value greater than the
     * current exponent */
    if (last_exponent == INT_MAX) {
      last_exponent = current_exponent + 1;
    } else {
      last_exponent = targets->exponent;
    }

    /* Compute the adjusted_exponent in such a way, that it is always diminishing in value. The adjusted
     * exponent can only be used to verify if a new target has been hit. To compute the actual target
     * value, the current_exponent needs to be used. */
    adjusted_exponent = 2 * (int) (ceil(log10(targets->precision / 10.0) * number_of_targets_double))
        - current_exponent - 1;

    if (adjusted_exponent < last_exponent) {
      /* Update the target information */
      targets->exponent = adjusted_exponent;
      targets->value = - pow(10, (double) current_exponent / number_of_targets_double);
      update_performed = 1;
    }
  }

  return update_performed;
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding triggers based on numbers of evaluations.
 */
/**@{*/

/**
 * @brief Creates and returns a structure containing information on triggers based on evaluation numbers.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 *
 * @note The coco_observer_evaluations_t object instances need to be freed using the
 * coco_observer_evaluations_free function!
 *
 * @param base_evaluations Evaluation numbers formatted as a string, which are used as the base to compute
 * the second trigger. For example, if base_evaluations = "1,2,5", the logger will be triggered by
 * evaluations dim*1, dim*2, dim*5, 10*dim*1, 10*dim*2, 10*dim*5, 100*dim*1, 100*dim*2, 100*dim*5, ...
 */
static coco_observer_evaluations_t *coco_observer_evaluations(const char *base_evaluations,
                                                              const size_t dimension) {

  coco_observer_evaluations_t *evaluations = (coco_observer_evaluations_t *) coco_allocate_memory(
      sizeof(*evaluations));

  /* First trigger */
  evaluations->value1 = 1;
  evaluations->exponent1 = 0;
  evaluations->number_of_triggers = 20;

  /* Second trigger */
  evaluations->base_evaluations = coco_string_parse_ranges(base_evaluations, 1, 0, "base_evaluations",
      COCO_MAX_EVALS_TO_LOG);
  evaluations->dimension = dimension;
  evaluations->base_count = coco_count_numbers(evaluations->base_evaluations, COCO_MAX_EVALS_TO_LOG,
      "base_evaluations");
  evaluations->base_index = 0;
  evaluations->value2 = dimension * evaluations->base_evaluations[0];
  evaluations->exponent2 = 0;

  return evaluations;
}

/**
 * @brief Computes and returns whether the given evaluation number triggers the first condition of the
 * logging based on the number of evaluations.
 *
 * The second condition is:
 * evaluation_number == 10**(exponent1/number_of_triggers)
 */
static int coco_observer_evaluations_trigger_first(coco_observer_evaluations_t *evaluations,
                                                   const size_t evaluation_number) {

  assert(evaluations != NULL);

  if (evaluation_number == evaluations->value1) {
    /* Compute the next value for the first trigger */
    while (coco_double_to_size_t(floor(pow(10, (double) evaluations->exponent1 / (double) evaluations->number_of_triggers)) <= evaluations->value1)) {
      evaluations->exponent1++;
    }
    evaluations->value1 = coco_double_to_size_t(floor(pow(10, (double) evaluations->exponent1 / (double) evaluations->number_of_triggers)));
    return 1;
  }
  return 0;
}

/**
 * @brief Computes and returns whether the given evaluation number triggers the second condition of the
 * logging based on the number of evaluations.
 *
 * The second condition is:
 * evaluation_number == base_evaluation[base_index] * dimension * (10**exponent2)
 */
static int coco_observer_evaluations_trigger_second(coco_observer_evaluations_t *evaluations,
                                                    const size_t evaluation_number) {

  assert(evaluations != NULL);

  if (evaluation_number == evaluations->value2) {
    /* Compute the next value for the second trigger */
    if (evaluations->base_index < evaluations->base_count - 1) {
      evaluations->base_index++;
    } else {
      evaluations->base_index = 0;
      evaluations->exponent2++;
    }
    evaluations->value2 = coco_double_to_size_t(pow(10, (double) evaluations->exponent2)
        * (double) (long) evaluations->dimension
        * (double) (long) evaluations->base_evaluations[evaluations->base_index]);
    return 1;
  }
  return 0;
}

/**
 * @brief Returns 1 if any of the two triggers based on the number of evaluations equal 1 and 0 otherwise.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
static int coco_observer_evaluations_trigger(coco_observer_evaluations_t *evaluations,
                                             const size_t evaluation_number) {

  /* Both functions need to be called so that both triggers are correctly updated */
  int first = coco_observer_evaluations_trigger_first(evaluations, evaluation_number);
  int second = coco_observer_evaluations_trigger_second(evaluations, evaluation_number);

  return (first + second > 0) ? 1: 0;
}

/**
 * @brief Frees the given evaluations object.
 */
static void coco_observer_evaluations_free(coco_observer_evaluations_t *evaluations) {

  assert(evaluations != NULL);
  coco_free_memory(evaluations->base_evaluations);
  coco_free_memory(evaluations);
}

/**@}*/

/***********************************************************************************************************/

/**
 * @brief Allocates memory for a coco_observer_t instance.
 */
static coco_observer_t *coco_observer_allocate(const char *result_folder,
                                               const char *observer_name,
                                               const char *algorithm_name,
                                               const char *algorithm_info,
                                               const size_t number_target_triggers,
                                               const double target_precision,
                                               const size_t number_evaluation_triggers,
                                               const char *base_evaluation_triggers,
                                               const int precision_x,
                                               const int precision_f) {

  coco_observer_t *observer;
  observer = (coco_observer_t *) coco_allocate_memory(sizeof(*observer));
  /* Initialize fields to sane/safe defaults */
  observer->result_folder = coco_strdup(result_folder);
  observer->observer_name = coco_strdup(observer_name);
  observer->algorithm_name = coco_strdup(algorithm_name);
  observer->algorithm_info = coco_strdup(algorithm_info);
  observer->number_target_triggers = number_target_triggers;
  observer->target_precision = target_precision;
  observer->number_evaluation_triggers = number_evaluation_triggers;
  observer->base_evaluation_triggers = coco_strdup(base_evaluation_triggers);
  observer->precision_x = precision_x;
  observer->precision_f = precision_f;
  observer->data = NULL;
  observer->data_free_function = NULL;
  observer->logger_allocate_function = NULL;
  observer->logger_free_function = NULL;
  observer->is_active = 1;
  return observer;
}

void coco_observer_free(coco_observer_t *observer) {

  if (observer != NULL) {
    observer->is_active = 0;
    if (observer->observer_name != NULL)
      coco_free_memory(observer->observer_name);
    if (observer->result_folder != NULL)
      coco_free_memory(observer->result_folder);
    if (observer->algorithm_name != NULL)
      coco_free_memory(observer->algorithm_name);
    if (observer->algorithm_info != NULL)
      coco_free_memory(observer->algorithm_info);

    if (observer->base_evaluation_triggers != NULL)
      coco_free_memory(observer->base_evaluation_triggers);

    if (observer->data != NULL) {
      if (observer->data_free_function != NULL) {
        observer->data_free_function(observer->data);
      }
      coco_free_memory(observer->data);
      observer->data = NULL;
    }

    observer->logger_allocate_function = NULL;
    observer->logger_free_function = NULL;

    coco_free_memory(observer);
    observer = NULL;
  }
}

#line 1 "code-experiments/src/logger_bbob.c"
/**
 * @file logger_bbob.c
 * @brief Implementation of the bbob logger.
 *
 * Logs the performance of a single-objective optimizer on noisy or noiseless problems.
 * It produces four kinds of files:
 * - The "info" files ...
 * - The "dat" files ...
 * - The "tdat" files ...
 * - The "rdat" files ...
 */

/* TODO: Document this file in doxygen style! */

#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#line 23 "code-experiments/src/logger_bbob.c"

#line 25 "code-experiments/src/logger_bbob.c"
#line 26 "code-experiments/src/logger_bbob.c"
#line 27 "code-experiments/src/logger_bbob.c"
#line 1 "code-experiments/src/observer_bbob.c"
/**
 * @file observer_bbob.c
 * @brief Implementation of the bbob observer.
 */

#line 7 "code-experiments/src/observer_bbob.c"
#line 8 "code-experiments/src/observer_bbob.c"

static coco_problem_t *logger_bbob(coco_observer_t *observer, coco_problem_t *problem);
static void logger_bbob_free(void *logger);

/**
 * @brief The bbob observer data type.
 */
typedef struct {
  /* TODO: Can be used to store variables that need to be accessible during one run (i.e. for multiple
   * problems). For example, the following global variables from logger_bbob.c could be stored here: */
  size_t current_dim;
  size_t current_fun_id;
  /* ... and others */
} observer_bbob_data_t;

/**
 * @brief Initializes the bbob observer.
 */
static void observer_bbob(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys) {

  observer->logger_allocate_function = logger_bbob;
  observer->logger_free_function = logger_bbob_free;
  observer->data_free_function = NULL;
  observer->data = NULL;

  *option_keys = NULL;

  (void) options; /* To silence the compiler */
}
#line 28 "code-experiments/src/logger_bbob.c"

/*static const size_t bbob_nbpts_nbevals = 20; Wassim: tentative, are now observer options with these default values*/
/*static const size_t bbob_nbpts_fval = 5;*/
static size_t bbob_current_dim = 0;
static size_t bbob_current_funId = 0;
static size_t bbob_infoFile_firstInstance = 0;
char bbob_infoFile_firstInstance_char[3];
/* a possible solution: have a list of dims that are already in the file, if the ones we're about to log
 * is != bbob_current_dim and the funId is currend_funId, create a new .info file with as suffix the
 * number of the first instance */
static const int bbob_number_of_dimensions = 6;
static size_t bbob_dimensions_in_current_infoFile[6] = { 0, 0, 0, 0, 0, 0 }; /* TODO should use dimensions from the suite */

/* The current_... mechanism fails if several problems are open.
 * For the time being this should lead to an error.
 *
 * A possible solution: bbob_logger_is_open becomes a reference
 * counter and as long as another logger is open, always a new info
 * file is generated.
 * TODO: Shouldn't the new way of handling observers already fix this?
 */
static int bbob_logger_is_open = 0; /* this could become lock-list of .info files */

/* TODO: add possibility of adding a prefix to the index files (easy to do through observer options) */

/**
 * @brief The bbob logger data type.
 */
typedef struct {
  coco_observer_t *observer;
  int is_initialized;
  /*char *path;// relative path to the data folder. //Wassim: now fetched from the observer */
  /*const char *alg_name; the alg name, for now, temporarily the same as the path. Wassim: Now in the observer */
  FILE *index_file; /* index file */
  FILE *fdata_file; /* function value aligned data file */
  FILE *tdata_file; /* number of function evaluations aligned data file */
  FILE *rdata_file; /* restart info data file */
  size_t number_of_evaluations;
  double best_fvalue;
  double last_fvalue;
  short written_last_eval; /* allows writing the the data of the final fun eval in the .tdat file if not already written by the t_trigger*/
  double *best_solution;
  /* The following are to only pass data as a parameter in the free function. The
   * interface should probably be the same for all free functions so passing the
   * problem as a second parameter is not an option even though we need info
   * form it.*/
  size_t function_id; /*TODO: consider changing name*/
  size_t instance_id;
  size_t number_of_variables;
  double optimal_fvalue;

  coco_observer_targets_t *targets;          /**< @brief Triggers based on target values. */
  coco_observer_evaluations_t *evaluations;  /**< @brief Triggers based on the number of evaluations. */

} logger_bbob_data_t;

static const char *bbob_file_header_str = "%% function evaluation | "
    "noise-free fitness - Fopt (%13.12e) | "
    "best noise-free fitness - Fopt | "
    "measured fitness | "
    "best measured fitness | "
    "x1 | "
    "x2...\n";

/**
 * adds a formated line to a data file
 */
static void logger_bbob_write_data(FILE *target_file,
                                   size_t number_of_evaluations,
                                   double fvalue,
                                   double best_fvalue,
                                   double best_value,
                                   const double *x,
                                   size_t number_of_variables) {
  /* for some reason, it's %.0f in the old code instead of the 10.9e
   * in the documentation
   */
  fprintf(target_file, "%lu %+10.9e %+10.9e %+10.9e %+10.9e", (unsigned long) number_of_evaluations,
  		fvalue - best_value, best_fvalue - best_value, fvalue, best_fvalue);
  if (number_of_variables < 22) {
    size_t i;
    for (i = 0; i < number_of_variables; i++) {
      fprintf(target_file, " %+5.4e", x[i]);
    }
  }
  fprintf(target_file, "\n");
}

/**
 * Error when trying to create the file "path"
 */
static void logger_bbob_error_io(FILE *path, int errnum) {
  const char *error_format = "Error opening file: %s\n ";
  coco_error(error_format, strerror(errnum), path);
}

/**
 * Creates the data files or simply opens it
 */

/*
 calling sequence:
 logger_bbob_open_dataFile(&(logger->fdata_file), logger->observer->output_folder, dataFile_path,
 ".dat");
 */

static void logger_bbob_open_dataFile(FILE **target_file,
                                      const char *path,
                                      const char *dataFile_path,
                                      const char *file_extension) {
  char file_path[COCO_PATH_MAX] = { 0 };
  char relative_filePath[COCO_PATH_MAX] = { 0 };
  int errnum;
  strncpy(relative_filePath, dataFile_path,
  COCO_PATH_MAX - strlen(relative_filePath) - 1);
  strncat(relative_filePath, file_extension,
  COCO_PATH_MAX - strlen(relative_filePath) - 1);
  coco_join_path(file_path, sizeof(file_path), path, relative_filePath, NULL);
  if (*target_file == NULL) {
    *target_file = fopen(file_path, "a+");
    errnum = errno;
    if (*target_file == NULL) {
      logger_bbob_error_io(*target_file, errnum);
    }
  }
}

/*
static void logger_bbob_open_dataFile(FILE **target_file,
                                      const char *path,
                                      const char *dataFile_path,
                                      const char *file_extension) {
  char file_path[COCO_PATH_MAX] = { 0 };
  char relative_filePath[COCO_PATH_MAX] = { 0 };
  int errnum;
  strncpy(relative_filePath, dataFile_path,
  COCO_PATH_MAX - strlen(relative_filePath) - 1);
  strncat(relative_filePath, file_extension,
  COCO_PATH_MAX - strlen(relative_filePath) - 1);
  coco_join_path(file_path, sizeof(file_path), path, relative_filePath, NULL);
  if (*target_file == NULL) {
    *target_file = fopen(file_path, "a+");
    errnum = errno;
    if (*target_file == NULL) {
      _bbob_logger_error_io(*target_file, errnum);
    }
  }
}*/

/**
 * Creates the index file fileName_prefix+problem_id+file_extension in
 * folde_path
 */
static void logger_bbob_openIndexFile(logger_bbob_data_t *logger,
                                      const char *folder_path,
                                      const char *indexFile_prefix,
                                      const char *function_id,
                                      const char *dataFile_path,
                                      const char *suite_name) {
  /* to add the instance number TODO: this should be done outside to avoid redoing this for the .*dat files */
  char used_dataFile_path[COCO_PATH_MAX] = { 0 };
  int errnum, newLine; /* newLine is at 1 if we need a new line in the info file */
  char function_id_char[3]; /* TODO: consider adding them to logger */
  char file_name[COCO_PATH_MAX] = { 0 };
  char file_path[COCO_PATH_MAX] = { 0 };
  FILE **target_file;
  FILE *tmp_file;
  strncpy(used_dataFile_path, dataFile_path, COCO_PATH_MAX - strlen(used_dataFile_path) - 1);
  if (bbob_infoFile_firstInstance == 0) {
    bbob_infoFile_firstInstance = logger->instance_id;
  }
  sprintf(function_id_char, "%lu", (unsigned long) logger->function_id);
  sprintf(bbob_infoFile_firstInstance_char, "%lu", (unsigned long) bbob_infoFile_firstInstance);
  target_file = &(logger->index_file);
  tmp_file = NULL; /* to check whether the file already exists. Don't want to use target_file */
  strncpy(file_name, indexFile_prefix, COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, "_f", COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, function_id_char, COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, "_i", COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, bbob_infoFile_firstInstance_char, COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, ".info", COCO_PATH_MAX - strlen(file_name) - 1);
  coco_join_path(file_path, sizeof(file_path), folder_path, file_name, NULL);
  if (*target_file == NULL) {
    tmp_file = fopen(file_path, "r"); /* to check for existence */
    if ((tmp_file) && (bbob_current_dim == logger->number_of_variables)
        && (bbob_current_funId == logger->function_id)) {
        /* new instance of current funId and current dim */
      newLine = 0;
      *target_file = fopen(file_path, "a+");
      if (*target_file == NULL) {
        errnum = errno;
        logger_bbob_error_io(*target_file, errnum);
      }
      fclose(tmp_file);
    } else { /* either file doesn't exist (new funId) or new Dim */
      /* check that the dim was not already present earlier in the file, if so, create a new info file */
      if (bbob_current_dim != logger->number_of_variables) {
        int i, j;
        for (i = 0;
            i < bbob_number_of_dimensions && bbob_dimensions_in_current_infoFile[i] != 0
                && bbob_dimensions_in_current_infoFile[i] != logger->number_of_variables; i++) {
          ; /* checks whether dimension already present in the current infoFile */
        }
        if (i < bbob_number_of_dimensions && bbob_dimensions_in_current_infoFile[i] == 0) {
          /* new dimension seen for the first time */
          bbob_dimensions_in_current_infoFile[i] = logger->number_of_variables;
          newLine = 1;
        } else {
          if (i < bbob_number_of_dimensions) { /* dimension already present, need to create a new file */
            newLine = 0;
            file_path[strlen(file_path) - strlen(bbob_infoFile_firstInstance_char) - 7] = 0; /* truncate the instance part */
            bbob_infoFile_firstInstance = logger->instance_id;
            sprintf(bbob_infoFile_firstInstance_char, "%lu", (unsigned long) bbob_infoFile_firstInstance);
            strncat(file_path, "_i", COCO_PATH_MAX - strlen(file_name) - 1);
            strncat(file_path, bbob_infoFile_firstInstance_char, COCO_PATH_MAX - strlen(file_name) - 1);
            strncat(file_path, ".info", COCO_PATH_MAX - strlen(file_name) - 1);
          } else {/*we have all dimensions*/
            newLine = 1;
          }
          for (j = 0; j < bbob_number_of_dimensions; j++) { /* new info file, reinitialize list of dims */
            bbob_dimensions_in_current_infoFile[j] = 0;
          }
          bbob_dimensions_in_current_infoFile[i] = logger->number_of_variables;
        }
      } else {
        if ( bbob_current_funId != logger->function_id ) {
          /*new function in the same file */
          newLine = 1;
        }
      }
      *target_file = fopen(file_path, "a+"); /* in any case, we append */
      if (*target_file == NULL) {
        errnum = errno;
        logger_bbob_error_io(*target_file, errnum);
      }
      if (tmp_file) { /* File already exists, new dim so just a new line. Also, close the tmp_file */
        if (newLine) {
          fprintf(*target_file, "\n");
        }
        fclose(tmp_file);
      }

      fprintf(*target_file,
          "funcId = %d, DIM = %lu, Precision = %.3e, algId = '%s', coco_version = '%s'\n",
          (int) strtol(function_id, NULL, 10), (unsigned long) logger->number_of_variables,
          pow(10, -8), logger->observer->algorithm_name, coco_version);
      /* TODO: Use this once suite can be read by the postprocessing
      fprintf(*target_file,
          "suite = '%s', funcId = %d, DIM = %lu, Precision = %.3e, algId = '%s', coco_version = '%s'\n",
          suite_name, (int) strtol(function_id, NULL, 10), (unsigned long) logger->number_of_variables,
          pow(10, -8), logger->observer->algorithm_name, coco_version);
      */
      fprintf(*target_file, "%%\n");
      strncat(used_dataFile_path, "_i", COCO_PATH_MAX - strlen(used_dataFile_path) - 1);
      strncat(used_dataFile_path, bbob_infoFile_firstInstance_char,
      COCO_PATH_MAX - strlen(used_dataFile_path) - 1);
      fprintf(*target_file, "%s.dat", used_dataFile_path); /* dataFile_path does not have the extension */
      bbob_current_dim = logger->number_of_variables;
      bbob_current_funId = logger->function_id;
    }
  }
}

/**
 * Generates the different files and folder needed by the logger to store the
 * data if these don't already exist
 */
static void logger_bbob_initialize(logger_bbob_data_t *logger, coco_problem_t *inner_problem) {
  /*
   Creates/opens the data and index files
   */
  char dataFile_path[COCO_PATH_MAX] = { 0 }; /* relative path to the .dat file from where the .info file is */
  char folder_path[COCO_PATH_MAX] = { 0 };
  char *tmpc_funId; /* serves to extract the function id as a char *. There should be a better way of doing this! */
  char *tmpc_dim; /* serves to extract the dimension as a char *. There should be a better way of doing this! */
  char indexFile_prefix[10] = "bbobexp"; /* TODO (minor): make the prefix bbobexp a parameter that the user can modify */
  size_t str_length_funId, str_length_dim;
  
  str_length_funId = coco_double_to_size_t(bbob2009_fmax(1, ceil(log10((double) coco_problem_get_suite_dep_function(inner_problem)))));
  str_length_dim = coco_double_to_size_t(bbob2009_fmax(1, ceil(log10((double) inner_problem->number_of_variables))));
  tmpc_funId = coco_allocate_string(str_length_funId);
  tmpc_dim = coco_allocate_string(str_length_dim);

  assert(logger != NULL);
  assert(inner_problem != NULL);
  assert(inner_problem->problem_id != NULL);

  sprintf(tmpc_funId, "%lu", (unsigned long) coco_problem_get_suite_dep_function(inner_problem));
  sprintf(tmpc_dim, "%lu", (unsigned long) inner_problem->number_of_variables);

  /* prepare paths and names */
  strncpy(dataFile_path, "data_f", COCO_PATH_MAX);
  strncat(dataFile_path, tmpc_funId,
  COCO_PATH_MAX - strlen(dataFile_path) - 1);
  coco_join_path(folder_path, sizeof(folder_path), logger->observer->result_folder, dataFile_path,
  NULL);
  coco_create_directory(folder_path);
  strncat(dataFile_path, "/bbobexp_f",
  COCO_PATH_MAX - strlen(dataFile_path) - 1);
  strncat(dataFile_path, tmpc_funId,
  COCO_PATH_MAX - strlen(dataFile_path) - 1);
  strncat(dataFile_path, "_DIM", COCO_PATH_MAX - strlen(dataFile_path) - 1);
  strncat(dataFile_path, tmpc_dim, COCO_PATH_MAX - strlen(dataFile_path) - 1);

  /* index/info file */
  assert(coco_problem_get_suite(inner_problem));
  logger_bbob_openIndexFile(logger, logger->observer->result_folder, indexFile_prefix, tmpc_funId,
      dataFile_path, coco_problem_get_suite(inner_problem)->suite_name);
  fprintf(logger->index_file, ", %lu", (unsigned long) coco_problem_get_suite_dep_instance(inner_problem));
  /* data files */
  /* TODO: definitely improvable but works for now */
  strncat(dataFile_path, "_i", COCO_PATH_MAX - strlen(dataFile_path) - 1);
  strncat(dataFile_path, bbob_infoFile_firstInstance_char,
  COCO_PATH_MAX - strlen(dataFile_path) - 1);
  logger_bbob_open_dataFile(&(logger->fdata_file), logger->observer->result_folder, dataFile_path, ".dat");
  fprintf(logger->fdata_file, bbob_file_header_str, logger->optimal_fvalue);

  logger_bbob_open_dataFile(&(logger->tdata_file), logger->observer->result_folder, dataFile_path, ".tdat");
  fprintf(logger->tdata_file, bbob_file_header_str, logger->optimal_fvalue);

  logger_bbob_open_dataFile(&(logger->rdata_file), logger->observer->result_folder, dataFile_path, ".rdat");
  fprintf(logger->rdata_file, bbob_file_header_str, logger->optimal_fvalue);
  logger->is_initialized = 1;
  coco_free_memory(tmpc_dim);
  coco_free_memory(tmpc_funId);
}

/**
 * Layer added to the transformed-problem evaluate_function by the logger
 */
static void logger_bbob_evaluate(coco_problem_t *problem, const double *x, double *y) {
  logger_bbob_data_t *logger = (logger_bbob_data_t *) coco_problem_transformed_get_data(problem);
  coco_problem_t * inner_problem = coco_problem_transformed_get_inner_problem(problem);

  size_t i;

  if (!logger->is_initialized) {
    logger_bbob_initialize(logger, inner_problem);
  }
  if ((coco_log_level >= COCO_DEBUG) && logger->number_of_evaluations == 0) {
    coco_debug("%4lu: ", (unsigned long) inner_problem->suite_dep_index);
    coco_debug("on problem %s ... ", coco_problem_get_id(inner_problem));
  }
  coco_evaluate_function(inner_problem, x, y);
  logger->last_fvalue = y[0];
  logger->written_last_eval = 0;
  if (logger->number_of_evaluations == 0 || y[0] < logger->best_fvalue) {
    logger->best_fvalue = y[0];
    for (i = 0; i < problem->number_of_variables; i++)
      logger->best_solution[i] = x[i];
  }
  logger->number_of_evaluations++;

  /* Add sanity check for optimal f value */
  assert(y[0] + 1e-13 >= logger->optimal_fvalue);

  /* Add a line in the .dat file for each logging target reached. */
    if (coco_observer_targets_trigger(logger->targets, y[0] - logger->optimal_fvalue)) {

    logger_bbob_write_data(logger->fdata_file, logger->number_of_evaluations, y[0], logger->best_fvalue,
        logger->optimal_fvalue, x, problem->number_of_variables);
  }

  /* Add a line in the .tdat file each time an fevals trigger is reached.*/
  if (coco_observer_evaluations_trigger(logger->evaluations, logger->number_of_evaluations)) {
    logger->written_last_eval = 1;
    logger_bbob_write_data(logger->tdata_file, logger->number_of_evaluations, y[0], logger->best_fvalue,
        logger->optimal_fvalue, x, problem->number_of_variables);
  }

  /* Flush output so that impatient users can see progress. */
  fflush(logger->fdata_file);
}

/**
 * Also serves as a finalize run method so. Must be called at the end
 * of Each run to correctly fill the index file
 *
 * TODO: make sure it is called at the end of each run or move the
 * writing into files to another function
 */
static void logger_bbob_free(void *stuff) {
  /* TODO: do all the "non simply freeing" stuff in another function
   * that can have problem as input
   */
  logger_bbob_data_t *logger = (logger_bbob_data_t *) stuff;

  if ((coco_log_level >= COCO_DEBUG) && logger && logger->number_of_evaluations > 0) {
    coco_debug("best f=%e after %lu fevals (done observing)\n", logger->best_fvalue,
    		(unsigned long) logger->number_of_evaluations);
  }
  if (logger->index_file != NULL) {
    fprintf(logger->index_file, ":%lu|%.1e", (unsigned long) logger->number_of_evaluations,
        logger->best_fvalue - logger->optimal_fvalue);
    fclose(logger->index_file);
    logger->index_file = NULL;
  }
  if (logger->fdata_file != NULL) {
    fclose(logger->fdata_file);
    logger->fdata_file = NULL;
  }
  if (logger->tdata_file != NULL) {
    /* TODO: make sure it handles restarts well. i.e., it writes
     * at the end of a single run, not all the runs on a given
     * instance. Maybe start with forcing it to generate a new
     * "instance" of problem for each restart in the beginning
     */
    if (!logger->written_last_eval) {
      logger_bbob_write_data(logger->tdata_file, logger->number_of_evaluations, logger->last_fvalue,
          logger->best_fvalue, logger->optimal_fvalue, logger->best_solution, logger->number_of_variables);
    }
    fclose(logger->tdata_file);
    logger->tdata_file = NULL;
  }

  if (logger->rdata_file != NULL) {
    fclose(logger->rdata_file);
    logger->rdata_file = NULL;
  }

  if (logger->best_solution != NULL) {
    coco_free_memory(logger->best_solution);
    logger->best_solution = NULL;
  }

  if (logger->targets != NULL){
    coco_free_memory(logger->targets);
    logger->targets = NULL;
  }

  if (logger->evaluations != NULL){
    coco_observer_evaluations_free(logger->evaluations);
    logger->evaluations = NULL;
  }

  bbob_logger_is_open = 0;
}

static coco_problem_t *logger_bbob(coco_observer_t *observer, coco_problem_t *inner_problem) {
  logger_bbob_data_t *logger_bbob;
  coco_problem_t *problem;

  logger_bbob = (logger_bbob_data_t *) coco_allocate_memory(sizeof(*logger_bbob));
  logger_bbob->observer = observer;

  if (inner_problem->number_of_objectives != 1) {
    coco_warning("logger_bbob(): The bbob logger shouldn't be used to log a problem with %d objectives",
        inner_problem->number_of_objectives);
  }

  if (bbob_logger_is_open)
    coco_error("The current bbob_logger (observer) must be closed before a new one is opened");
  /* This is the name of the folder which happens to be the algName */
  /*logger->path = coco_strdup(observer->output_folder);*/
  logger_bbob->index_file = NULL;
  logger_bbob->fdata_file = NULL;
  logger_bbob->tdata_file = NULL;
  logger_bbob->rdata_file = NULL;
  logger_bbob->number_of_variables = inner_problem->number_of_variables;
  if (inner_problem->best_value == NULL) {
    /* coco_error("Optimal f value must be defined for each problem in order for the logger to work properly"); */
    /* Setting the value to 0 results in the assertion y>=optimal_fvalue being susceptible to failure */
    coco_warning("undefined optimal f value. Set to 0");
    logger_bbob->optimal_fvalue = 0;
  } else {
    logger_bbob->optimal_fvalue = *(inner_problem->best_value);
  }

  logger_bbob->number_of_evaluations = 0;
  logger_bbob->best_solution = coco_allocate_vector(inner_problem->number_of_variables);
  /* TODO: the following inits are just to be in the safe side and
   * should eventually be removed. Some fields of the bbob_logger struct
   * might be useless
   */
  logger_bbob->function_id = coco_problem_get_suite_dep_function(inner_problem);
  logger_bbob->instance_id = coco_problem_get_suite_dep_instance(inner_problem);
  logger_bbob->written_last_eval = 1;
  logger_bbob->last_fvalue = DBL_MAX;
  logger_bbob->is_initialized = 0;

  /* Initialize triggers based on target values and number of evaluations */
  logger_bbob->targets = coco_observer_targets(observer->number_target_triggers, observer->target_precision);
  logger_bbob->evaluations = coco_observer_evaluations(observer->base_evaluation_triggers, inner_problem->number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, logger_bbob, logger_bbob_free, observer->observer_name);

  problem->evaluate_function = logger_bbob_evaluate;
  bbob_logger_is_open = 1;
  return problem;
}

#line 370 "code-experiments/src/coco_observer.c"
#line 1 "code-experiments/src/logger_biobj.c"
/**
 * @file logger_biobj.c
 * @brief Implementation of the bbob-biobj logger.
 *
 * Logs the values of the implemented indicators and archives nondominated solutions.
 * Produces four kinds of files:
 * - The "info" files contain high-level information on the performed experiment. One .info file is created
 * for each problem group (and indicator type) and contains information on all the problems in that problem
 * group (and indicator type).
 * - The "dat" files contain function evaluations, indicator values and target hits for every target hit as
 * well as for the last evaluation. One .dat file is created for each problem function and dimension (and
 * indicator type) and contains information for all instances of that problem (and indicator type).
 * - The "tdat" files contain function evaluation and indicator values for every predefined evaluation
 * number as well as for the last evaluation. One .tdat file is created for each problem function and
 * dimension (and indicator type) and contains information for all instances of that problem (and indicator
 * type).
 * - The "adat" files are archive files that contain function evaluations, 2 objectives and dim variables
 * for every nondominated solution. Whether these files are created, at what point in time the logger writes
 * nondominated solutions to the archive and whether the decision variables are output or not depends on
 * the values of log_nondom_mode and log_nondom_mode. See the bi-objective observer constructor
 * observer_biobj() for more information. One .adat file is created for each problem function, dimension
 * and instance.
 *
 * @note Whenever in this file a ROI is mentioned, it means the (normalized) region of interest in the
 * objective space. The non-normalized ROI is a rectangle with the ideal and nadir points as its two
 * opposite vertices, while the normalized ROI is the square [0, 1]^2. If not specifically mentioned, the
 * normalized ROI is assumed.
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#line 35 "code-experiments/src/logger_biobj.c"
#line 36 "code-experiments/src/logger_biobj.c"

#line 38 "code-experiments/src/logger_biobj.c"
#line 39 "code-experiments/src/logger_biobj.c"
#line 40 "code-experiments/src/logger_biobj.c"
#line 1 "code-experiments/src/mo_avl_tree.c"
/*****************************************************************************

 avl.c - Source code for libavl

 Copyright (c) 1998  Michael H. Buselli <cosine@cosine.org>
 Copyright (c) 2000-2009  Wessel Dankers <wsl@fruit.je>

 This file is part of libavl.

 libavl is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation, either version 3 of
 the License, or (at your option) any later version.

 libavl is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU General Public License
 and a copy of the GNU Lesser General Public License along with
 libavl.  If not, see <http://www.gnu.org/licenses/>.

 Augmented AVL-tree. Original by Michael H. Buselli <cosine@cosine.org>.

 Modified by Wessel Dankers <wsl@fruit.je> to add a bunch of bloat
 to the source code, change the interface and replace a few bugs.
 Mail him if you find any new bugs.

 Renamed and additionally modified by BOBBies to fit the COCO platform.

 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

/* In order to easily (un)comment unused functions */
#define AVL_TREE_COMMENT_UNUSED 1

/* We need either depths, counts or both (the latter being the default) */
#if !defined(AVL_DEPTH) && !defined(AVL_COUNT)
#define AVL_DEPTH
#define AVL_COUNT
#endif

/* User supplied function to compare two items like strcmp() does.
 * For example: compare(a,b) will return:
 *   -1  if a < b
 *    0  if a = b
 *    1  if a > b
 */
typedef int (*avl_compare_t)(const void *a, const void *b, void *userdata);

/* User supplied function to delete an item when a node is free()d.
 * If NULL, the item is not free()d.
 */
typedef void (*avl_free_t)(void *item, void *userdata);

#define AVL_CMP(a,b) ((a) < (b) ? -1 : (a) != (b))

#if defined(AVL_COUNT) && defined(AVL_DEPTH)
#define AVL_NODE_INITIALIZER(item) { 0, 0, 0, 0, 0, (item), 0, 0 }
#else
#define AVL_NODE_INITIALIZER(item) { 0, 0, 0, 0, 0, (item), 0 }
#endif

typedef struct avl_node {
  struct avl_node *prev;
  struct avl_node *next;
  struct avl_node *parent;
  struct avl_node *left;
  struct avl_node *right;
  void *item;
#ifdef AVL_COUNT
  unsigned long count;
#endif
#ifdef AVL_DEPTH
  unsigned char depth;
#endif
} avl_node_t;

#define AVL_TREE_INITIALIZER(cmp, free) { 0, 0, 0, (cmp), (free), {0}, 0, 0 }

typedef struct avl_tree {
  avl_node_t *top;
  avl_node_t *head;
  avl_node_t *tail;
  avl_compare_t cmpitem;
  avl_free_t freeitem;
  void *userdata;
  struct avl_allocator *allocator;
} avl_tree_t;

#define AVL_ALLOCATOR_INITIALIZER(alloc, dealloc) { (alloc), (dealloc) }

typedef avl_node_t *(*avl_allocate_t)(struct avl_allocator *);
typedef void (*avl_deallocate_t)(struct avl_allocator *, avl_node_t *);

typedef struct avl_allocator {
  avl_allocate_t allocate;
  avl_deallocate_t deallocate;
} avl_allocator_t;

static void avl_rebalance(avl_tree_t *, avl_node_t *);
static avl_node_t *avl_node_insert_after(avl_tree_t *avltree, avl_node_t *node, avl_node_t *newnode);

#ifdef AVL_COUNT
#define AVL_NODE_COUNT(n)  ((n) ? (n)->count : 0)
#define AVL_L_COUNT(n)     (AVL_NODE_COUNT((n)->left))
#define AVL_R_COUNT(n)     (AVL_NODE_COUNT((n)->right))
#define AVL_CALC_COUNT(n)  (AVL_L_COUNT(n) + AVL_R_COUNT(n) + 1)
#endif

#ifdef AVL_DEPTH
#define AVL_NODE_DEPTH(n)  ((n) ? (n)->depth : 0)
#define AVL_L_DEPTH(n)     (AVL_NODE_DEPTH((n)->left))
#define AVL_R_DEPTH(n)     (AVL_NODE_DEPTH((n)->right))
#define AVL_CALC_DEPTH(n)  ((unsigned char)((AVL_L_DEPTH(n) > AVL_R_DEPTH(n) ? AVL_L_DEPTH(n) : AVL_R_DEPTH(n)) + 1))
#endif

const avl_node_t avl_node_0 = { 0, 0, 0, 0, 0, 0, 0, 0 };
const avl_tree_t avl_tree_0 = { 0, 0, 0, 0, 0, 0, 0 };
const avl_allocator_t avl_allocator_0 = { 0, 0 };

#define AVL_CONST_NODE(x) ((avl_node_t *)(x))
#define AVL_CONST_ITEM(x) ((void *)(x))

static int avl_check_balance(avl_node_t *avlnode) {
#ifdef AVL_DEPTH
  int d;
  d = AVL_R_DEPTH(avlnode) - AVL_L_DEPTH(avlnode);
  return d < -1 ? -1 : d > 1;
#else
  /*  int d;
   *  d = ffs(AVL_R_COUNT(avl_node)) - ffs(AVL_L_COUNT(avl_node));
   *  d = d < -1 ? -1 : d > 1;
   */
#ifdef AVL_COUNT
  int pl, r;

  pl = ffs(AVL_L_COUNT(avlnode));
  r = AVL_R_COUNT(avlnode);

  if (r >> pl + 1)
  return 1;
  if (pl < 2 || r >> pl - 2)
  return 0;
  return -1;
#else
#error No balancing possible.
#endif
#endif
}

#ifdef AVL_COUNT
static unsigned long avl_count(const avl_tree_t *avltree) {
  if (!avltree)
    return 0;
  return AVL_NODE_COUNT(avltree->top);
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_at(const avl_tree_t *avltree, unsigned long index) {
  avl_node_t *avlnode;
  unsigned long c;

  if (!avltree)
    return NULL;

  avlnode = avltree->top;

  while (avlnode) {
    c = AVL_L_COUNT(avlnode);

    if (index < c) {
      avlnode = avlnode->left;
    } else if (index > c) {
      avlnode = avlnode->right;
      index -= c + 1;
    } else {
      return avlnode;
    }
  }
  return NULL;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static unsigned long avl_index(const avl_node_t *avlnode) {
  avl_node_t *next;
  unsigned long c;

  if (!avlnode)
    return 0;

  c = AVL_L_COUNT(avlnode);

  while ((next = avlnode->parent)) {
    if (avlnode == next->right)
      c += AVL_L_COUNT(next) + 1;
    avlnode = next;
  }

  return c;
}
#endif
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static const avl_node_t *avl_search_leftmost_equal(const avl_tree_t *tree, const avl_node_t *node,
    const void *item) {
  avl_compare_t cmp = tree->cmpitem;
  void *userdata = tree->userdata;
  const avl_node_t *r = node;

  for (;;) {
    for (;;) {
      node = node->left;
      if (!node)
        return r;
      if (cmp(item, node->item, userdata))
        break;
      r = node;
    }
    for (;;) {
      node = node->right;
      if (!node)
        return r;
      if (!cmp(item, node->item, userdata))
        break;
    }
    r = node;
  }

  return NULL; /* To silence the compiler */

}
#endif

static const avl_node_t *avl_search_rightmost_equal(const avl_tree_t *tree,
                                                    const avl_node_t *node,
                                                    const void *item) {
  avl_compare_t cmp = tree->cmpitem;
  void *userdata = tree->userdata;
  const avl_node_t *r = node;

  for (;;) {
    for (;;) {
      node = node->right;
      if (!node)
        return r;
      if (cmp(item, node->item, userdata))
        break;
      r = node;
    }
    for (;;) {
      node = node->left;
      if (!node)
        return r;
      if (!cmp(item, node->item, userdata))
        break;
    }
    r = node;
  }

  return NULL; /* To silence the compiler */
}

/* Searches for an item, returning either some exact
 * match, or (if no exact match could be found) the first (leftmost)
 * of the nodes that have an item larger than the search item.
 * If exact is not NULL, *exact will be set to:
 *    0  if the returned node is unequal or NULL
 *    1  if the returned node is equal
 * Returns NULL if no equal or larger element could be found.
 * O(lg n) */
#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_search_leftish(const avl_tree_t *tree, const void *item, int *exact) {
  avl_node_t *node;
  avl_compare_t cmp;
  void *userdata;
  int c;

  if (!exact)
    exact = &c;

  if (!tree)
    return *exact = 0, (avl_node_t *) NULL;

  node = tree->top;
  if (!node)
    return *exact = 0, (avl_node_t *) NULL;

  cmp = tree->cmpitem;
  userdata = tree->userdata;

  for (;;) {
    c = cmp(item, node->item, userdata);

    if (c < 0) {
      if (node->left)
        node = node->left;
      else
        return *exact = 0, node;
    } else if (c > 0) {
      if (node->right)
        node = node->right;
      else
        return *exact = 0, node->next;
    } else {
      return *exact = 1, node;
    }
  }

  return NULL; /* To silence the compiler */

}
#endif

/* Searches for an item, returning either some exact
 * match, or (if no exact match could be found) the last (rightmost)
 * of the nodes that have an item smaller than the search item.
 * If exact is not NULL, *exact will be set to:
 *    0  if the returned node is unequal or NULL
 *    1  if the returned node is equal
 * Returns NULL if no equal or smaller element could be found.
 * O(lg n) */
static avl_node_t *avl_search_rightish(const avl_tree_t *tree, const void *item, int *exact) {
  avl_node_t *node;
  avl_compare_t cmp;
  void *userdata;
  int c;

  if (!exact)
    exact = &c;

  if (!tree)
    return *exact = 0, (avl_node_t *) NULL;

  node = tree->top;
  if (!node)
    return *exact = 0, (avl_node_t *) NULL;

  cmp = tree->cmpitem;
  userdata = tree->userdata;

  for (;;) {
    c = cmp(item, node->item, userdata);

    if (c < 0) {
      if (node->left)
        node = node->left;
      else
        return *exact = 0, node->prev;
    } else if (c > 0) {
      if (node->right)
        node = node->right;
      else
        return *exact = 0, node;
    } else {
      return *exact = 1, node;
    }
  }

  return NULL; /* To silence the compiler */
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_search_left(const avl_tree_t *tree, const void *item, int *exact) {
  avl_node_t *node;
  int c;

  if (!exact)
    exact = &c;

  if (!tree)
    return *exact = 0, (avl_node_t *) NULL;

  node = avl_search_leftish(tree, item, exact);
  if (*exact)
    return AVL_CONST_NODE(avl_search_leftmost_equal(tree, node, item));

  return AVL_CONST_NODE(node);
}
#endif

/* Searches for an item, returning either the last (rightmost) exact
 * match, or (if no exact match could be found) the last (rightmost)
 * of the nodes that have an item smaller than the search item.
 * If exact is not NULL, *exact will be set to:
 *    0  if the returned node is inequal or NULL
 *    1  if the returned node is equal
 * Returns NULL if no equal or smaller element could be found.
 * O(lg n) */
static avl_node_t *avl_item_search_right(const avl_tree_t *tree, const void *item, int *exact) {
  const avl_node_t *node;
  int c;

  if (!exact)
    exact = &c;

  node = avl_search_rightish(tree, item, exact);
  if (*exact)
    return AVL_CONST_NODE(avl_search_rightmost_equal(tree, node, item));

  return AVL_CONST_NODE(node);
}

/* Searches for the item in the tree and returns a matching node if found
 * or NULL if not.
 * O(lg n) */
static avl_node_t *avl_item_search(const avl_tree_t *avltree, const void *item) {
  int c;
  avl_node_t *n;
  n = avl_search_rightish(avltree, item, &c);
  return c ? n : NULL;
}

/* Initializes a new tree for elements that will be ordered using
 * the supplied strcmp()-like function.
 * Returns the value of avltree (even if it's NULL).
 * O(1) */
static avl_tree_t *avl_tree_init(avl_tree_t *avltree, avl_compare_t cmp, avl_free_t free) {
  if (avltree) {
    avltree->head = NULL;
    avltree->tail = NULL;
    avltree->top = NULL;
    avltree->cmpitem = cmp;
    avltree->freeitem = free;
    avltree->userdata = NULL;
    avltree->allocator = NULL;
  }
  return avltree;
}

/* Allocates and initializes a new tree for elements that will be
 * ordered using the supplied strcmp()-like function.
 * Returns NULL if memory could not be allocated.
 * O(1) */
static avl_tree_t *avl_tree_construct(avl_compare_t cmp, avl_free_t free) {
  return avl_tree_init((avl_tree_t *) malloc(sizeof(avl_tree_t)), cmp, free);
}

/* Reinitializes the tree structure for reuse. Nothing is free()d.
 * Compare and free functions are left alone.
 * Returns the value of avltree (even if it's NULL).
 * O(1) */
static avl_tree_t *avl_tree_clear(avl_tree_t *avltree) {
  if (avltree)
    avltree->top = avltree->head = avltree->tail = NULL;
  return avltree;
}

static void avl_node_free(avl_tree_t *avltree, avl_node_t *node) {
  avl_allocator_t *allocator;
  avl_deallocate_t deallocate;

  allocator = avltree->allocator;
  if (allocator) {
    deallocate = allocator->deallocate;
    if (deallocate)
      deallocate(allocator, node);
  } else {
    free(node);
  }
}

/* Free()s all nodes in the tree but leaves the tree itself.
 * If the tree's free is not NULL it will be invoked on every item.
 * Returns the value of avltree (even if it's NULL).
 * O(n) */
static avl_tree_t *avl_tree_purge(avl_tree_t *avltree) {
  avl_node_t *node, *next;
  avl_free_t func;
  avl_allocator_t *allocator;
  avl_deallocate_t deallocate;
  void *userdata;

  if (!avltree)
    return NULL;

  userdata = avltree->userdata;

  func = avltree->freeitem;
  allocator = avltree->allocator;
  deallocate = allocator ? allocator->deallocate : (avl_deallocate_t) NULL;

  for (node = avltree->head; node; node = next) {
    next = node->next;
    if (func)
      func(node->item, userdata);
    if (allocator) {
      if (deallocate)
        deallocate(allocator, node);
    } else {
      free(node);
    }
  }

  return avl_tree_clear(avltree);
}

/* Frees the entire tree efficiently. Nodes will be free()d.
 * If the tree's free is not NULL it will be invoked on every item.
 * O(n) */
static void avl_tree_destruct(avl_tree_t *avltree) {
  if (!avltree)
    return;
  (void) avl_tree_purge(avltree);
  free(avltree);
}

static void avl_node_clear(avl_node_t *newnode) {
  newnode->left = newnode->right = NULL;
#   ifdef AVL_COUNT
  newnode->count = 1;
#   endif
#   ifdef AVL_DEPTH
  newnode->depth = 1;
#   endif
}

/* Initializes memory for use as a node.
 * Returns the value of avlnode (even if it's NULL).
 * O(1) */
static avl_node_t *avl_node_init(avl_node_t *newnode, const void *item) {
  if (newnode)
    newnode->item = AVL_CONST_ITEM(item);
  return newnode;
}

/* Allocates and initializes memory for use as a node.
 * Returns the value of avlnode (or NULL if the allocation failed).
 * O(1) */
static avl_node_t *avl_alloc(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;
  avl_allocator_t *allocator = avltree ? avltree->allocator : (avl_allocator_t *) NULL;
  avl_allocate_t allocate;
  if (allocator) {
    allocate = allocator->allocate;
    if (allocator) {
      newnode = allocate(allocator);
    } else {
      errno = ENOSYS;
      newnode = NULL;
    }
  } else {
    newnode = (avl_node_t *) malloc(sizeof *newnode);
  }
  return avl_node_init(newnode, item);
}

/* Insert a node in an empty tree. If avl_node is NULL, the tree will be
 * cleared and ready for re-use.
 * If the tree is not empty, the old nodes are left dangling.
 * O(1) */
static avl_node_t *avl_insert_top(avl_tree_t *avltree, avl_node_t *newnode) {
  avl_node_clear(newnode);
  newnode->prev = newnode->next = newnode->parent = NULL;
  avltree->head = avltree->tail = avltree->top = newnode;
  return newnode;
}

/* Insert a node before another node. Returns the new node.
 * If old is NULL, the item is appended to the tree.
 * O(lg n) */
static avl_node_t *avl_node_insert_before(avl_tree_t *avltree, avl_node_t *node, avl_node_t *newnode) {
  if (!avltree || !newnode)
    return NULL;

  if (!node)
    return
        avltree->tail ?
            avl_node_insert_after(avltree, avltree->tail, newnode) : avl_insert_top(avltree, newnode);

  if (node->left)
    return avl_node_insert_after(avltree, node->prev, newnode);

  avl_node_clear(newnode);

  newnode->next = node;
  newnode->parent = node;

  newnode->prev = node->prev;
  if (node->prev)
    node->prev->next = newnode;
  else
    avltree->head = newnode;
  node->prev = newnode;

  node->left = newnode;
  avl_rebalance(avltree, node);
  return newnode;
}

/* Insert a node after another node. Returns the new node.
 * If old is NULL, the item is prepended to the tree.
 * O(lg n) */
static avl_node_t *avl_node_insert_after(avl_tree_t *avltree, avl_node_t *node, avl_node_t *newnode) {
  if (!avltree || !newnode)
    return NULL;

  if (!node)
    return
        avltree->head ?
            avl_node_insert_before(avltree, avltree->head, newnode) : avl_insert_top(avltree, newnode);

  if (node->right)
    return avl_node_insert_before(avltree, node->next, newnode);

  avl_node_clear(newnode);

  newnode->prev = node;
  newnode->parent = node;

  newnode->next = node->next;
  if (node->next)
    node->next->prev = newnode;
  else
    avltree->tail = newnode;
  node->next = newnode;

  node->right = newnode;
  avl_rebalance(avltree, node);
  return newnode;
}

/* Insert a node into the tree and return it.
 * Returns NULL if an equal node is already in the tree.
 * O(lg n) */
static avl_node_t *avl_node_insert(avl_tree_t *avltree, avl_node_t *newnode) {
  avl_node_t *node;
  int c;

  node = avl_search_rightish(avltree, newnode->item, &c);
  return c ? NULL : avl_node_insert_after(avltree, node, newnode);
}


#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_insert_left(avl_tree_t *avltree, avl_node_t *newnode) {
  return avl_node_insert_before(avltree, avl_item_search_left(avltree, newnode->item, NULL), newnode);
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_insert_right(avl_tree_t *avltree, avl_node_t *newnode) {
  return avl_node_insert_after(avltree, avl_item_search_right(avltree, newnode->item, NULL), newnode);
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_insert_somewhere(avl_tree_t *avltree, avl_node_t *newnode) {
  return avl_node_insert_after(avltree, avl_search_rightish(avltree, newnode->item, NULL), newnode);
}
#endif

/* Insert an item into the tree and return the new node.
 * Returns NULL and sets errno if memory for the new node could not be
 * allocated or if the node is already in the tree (EEXIST).
 * O(lg n) */
static avl_node_t *avl_item_insert(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree)
    return errno = EFAULT, (avl_node_t *) NULL;

  newnode = avl_alloc(avltree, item);
  if (newnode) {
    if (avl_node_insert(avltree, newnode))
      return newnode;
    avl_node_free(avltree, newnode);
    errno = EEXIST;
  }
  return NULL;
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_somewhere(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree)
    return errno = EFAULT, (avl_node_t *) NULL;

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_somewhere(avltree, newnode);
  return NULL;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_before(avl_tree_t *avltree, avl_node_t *node, const void *item) {
  avl_node_t *newnode;

  if (!avltree)
    return errno = EFAULT, (avl_node_t *) NULL;

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_before(avltree, node, newnode);
  return NULL;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_after(avl_tree_t *avltree, avl_node_t *node, const void *item) {
  avl_node_t *newnode;

  if (!avltree)
    return errno = EFAULT, (avl_node_t *) NULL;

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_after(avltree, node, newnode);
  return NULL;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_left(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree)
    return errno = EFAULT, (avl_node_t *) NULL;

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_left(avltree, newnode);
  return NULL;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_right(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree)
    return errno = EFAULT, (avl_node_t *) NULL;

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_right(avltree, newnode);
  return NULL;
}
#endif

/* Deletes a node from the tree.
 * Returns the value of the node (even if it's NULL).
 * The item will NOT be free()d regardless of the tree's free handler.
 * This function comes in handy if you need to update the search key.
 * O(lg n) */
static avl_node_t *avl_node_unlink(avl_tree_t *avltree, avl_node_t *avlnode) {
  avl_node_t *parent;
  avl_node_t **superparent;
  avl_node_t *subst, *left, *right;
  avl_node_t *balnode;

  if (!avltree || !avlnode)
    return NULL;

  if (avlnode->prev)
    avlnode->prev->next = avlnode->next;
  else
    avltree->head = avlnode->next;

  if (avlnode->next)
    avlnode->next->prev = avlnode->prev;
  else
    avltree->tail = avlnode->prev;

  parent = avlnode->parent;

  superparent = parent ? avlnode == parent->left ? &parent->left : &parent->right : &avltree->top;

  left = avlnode->left;
  right = avlnode->right;
  if (!left) {
    *superparent = right;
    if (right)
      right->parent = parent;
    balnode = parent;
  } else if (!right) {
    *superparent = left;
    left->parent = parent;
    balnode = parent;
  } else {
    subst = avlnode->prev;
    if (subst == left) {
      balnode = subst;
    } else {
      balnode = subst->parent;
      balnode->right = subst->left;
      if (balnode->right)
        balnode->right->parent = balnode;
      subst->left = left;
      left->parent = subst;
    }
    subst->right = right;
    subst->parent = parent;
    right->parent = subst;
    *superparent = subst;
  }

  avl_rebalance(avltree, balnode);

  return avlnode;
}

/* Deletes a node from the tree. Returns immediately if the node is NULL.
 * If the tree's free is not NULL, it is invoked on the item.
 * If it is, returns the item. In all other cases returns NULL.
 * O(lg n) */
static void *avl_node_delete(avl_tree_t *avltree, avl_node_t *avlnode) {
  void *item = NULL;
  if (avlnode) {
    item = avlnode->item;
    (void) avl_node_unlink(avltree, avlnode);
    if (avltree->freeitem)
      avltree->freeitem(item, avltree->userdata);
    avl_node_free(avltree, avlnode);
  }
  return item;
}

/* Searches for an item in the tree and deletes it if found.
 * If the tree's free is not NULL, it is invoked on the item.
 * If it is, returns the item. In all other cases returns NULL.
 * O(lg n) */
static void *avl_item_delete(avl_tree_t *avltree, const void *item) {
  return avl_node_delete(avltree, avl_item_search(avltree, item));
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_fixup(avl_tree_t *avltree, avl_node_t *newnode) {
  avl_node_t *oldnode = NULL, *node;

  if (!avltree || !newnode)
    return NULL;

  node = newnode->prev;
  if (node) {
    oldnode = node->next;
    node->next = newnode;
  } else {
    avltree->head = newnode;
  }

  node = newnode->next;
  if (node) {
    oldnode = node->prev;
    node->prev = newnode;
  } else {
    avltree->tail = newnode;
  }

  node = newnode->parent;
  if (node) {
    if (node->left == oldnode)
      node->left = newnode;
    else
      node->right = newnode;
  } else {
    oldnode = avltree->top;
    avltree->top = newnode;
  }

  return oldnode;
}
#endif

/**
 * avl_rebalance:
 * Rebalances the tree if one side becomes too heavy.  This function
 * assumes that both subtrees are AVL-trees with consistent data.  The
 * function has the additional side effect of recalculating the count of
 * the tree at this node.  It should be noted that at the return of this
 * function, if a rebalance takes place, the top of this subtree is no
 * longer going to be the same node.
 */
static void avl_rebalance(avl_tree_t *avltree, avl_node_t *avlnode) {
  avl_node_t *child;
  avl_node_t *gchild;
  avl_node_t *parent;
  avl_node_t **superparent;

  parent = avlnode;

  while (avlnode) {
    parent = avlnode->parent;

    superparent = parent ? avlnode == parent->left ? &parent->left : &parent->right : &avltree->top;

    switch (avl_check_balance(avlnode)) {
    case -1:
      child = avlnode->left;
#           ifdef AVL_DEPTH
      if (AVL_L_DEPTH(child) >= AVL_R_DEPTH(child)) {
#           else
#           ifdef AVL_COUNT
        if (AVL_L_COUNT(child) >= AVL_R_COUNT(child)) {
#           else
#           error No balancing possible.
#           endif
#           endif
        avlnode->left = child->right;
        if (avlnode->left)
          avlnode->left->parent = avlnode;
        child->right = avlnode;
        avlnode->parent = child;
        *superparent = child;
        child->parent = parent;
#               ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
#               endif
#               ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
#               endif
      } else {
        gchild = child->right;
        avlnode->left = gchild->right;
        if (avlnode->left)
          avlnode->left->parent = avlnode;
        child->right = gchild->left;
        if (child->right)
          child->right->parent = child;
        gchild->right = avlnode;
        if (gchild->right)
          gchild->right->parent = gchild;
        gchild->left = child;
        if (gchild->left)
          gchild->left->parent = gchild;
        *superparent = gchild;
        gchild->parent = parent;
#               ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
        gchild->count = AVL_CALC_COUNT(gchild);
#               endif
#               ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
        gchild->depth = AVL_CALC_DEPTH(gchild);
#               endif
      }
      break;
    case 1:
      child = avlnode->right;
#           ifdef AVL_DEPTH
      if (AVL_R_DEPTH(child) >= AVL_L_DEPTH(child)) {
#           else
#           ifdef AVL_COUNT
        if (AVL_R_COUNT(child) >= AVL_L_COUNT(child)) {
#           else
#           error No balancing possible.
#           endif
#           endif
        avlnode->right = child->left;
        if (avlnode->right)
          avlnode->right->parent = avlnode;
        child->left = avlnode;
        avlnode->parent = child;
        *superparent = child;
        child->parent = parent;
#               ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
#               endif
#               ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
#               endif
      } else {
        gchild = child->left;
        avlnode->right = gchild->left;
        if (avlnode->right)
          avlnode->right->parent = avlnode;
        child->left = gchild->right;
        if (child->left)
          child->left->parent = child;
        gchild->left = avlnode;
        if (gchild->left)
          gchild->left->parent = gchild;
        gchild->right = child;
        if (gchild->right)
          gchild->right->parent = gchild;
        *superparent = gchild;
        gchild->parent = parent;
#               ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
        gchild->count = AVL_CALC_COUNT(gchild);
#               endif
#               ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
        gchild->depth = AVL_CALC_DEPTH(gchild);
#               endif
      }
      break;
    default:
#           ifdef AVL_COUNT
      avlnode->count = AVL_CALC_COUNT(avlnode);
#           endif
#           ifdef AVL_DEPTH
      avlnode->depth = AVL_CALC_DEPTH(avlnode);
#           endif
    }
    avlnode = parent;
  }
}
#line 41 "code-experiments/src/logger_biobj.c"
#line 1 "code-experiments/src/observer_biobj.c"
/**
 * @file observer_biobj.c
 * @brief Implementation of the bbob-biobj observer.
 */

#line 7 "code-experiments/src/observer_biobj.c"
#line 8 "code-experiments/src/observer_biobj.c"

#line 10 "code-experiments/src/observer_biobj.c"
#line 11 "code-experiments/src/observer_biobj.c"

/** @brief Enum for denoting the way in which the nondominated solutions are treated. */
typedef enum {
  LOG_NONDOM_NONE, LOG_NONDOM_FINAL, LOG_NONDOM_ALL, LOG_NONDOM_READ
} observer_biobj_log_nondom_e;

/** @brief Enum for denoting the when the decision variables are logged. */
typedef enum {
  LOG_VARS_NEVER, LOG_VARS_LOW_DIM, LOG_VARS_ALWAYS
} observer_biobj_log_vars_e;

/**
 * @brief The bbob-biobj observer data type.
 */
typedef struct {
  observer_biobj_log_nondom_e log_nondom_mode; /**< @brief Handling of the nondominated solutions. */
  observer_biobj_log_vars_e log_vars_mode;     /**< @brief When the decision variables are logged. */

  int compute_indicators;                      /**< @brief Whether to compute indicators. */
  int produce_all_data;                        /**< @brief Whether to produce all data. */

  long previous_function;                      /**< @brief Function of the previous logged problem. */
  long previous_dimension;                     /**< @brief Dimension of the previous logged problem */

} observer_biobj_data_t;

static coco_problem_t *logger_biobj(coco_observer_t *observer, coco_problem_t *problem);
static void logger_biobj_free(void *logger);

/**
 * @brief Initializes the bi-objective observer.
 *
 * Possible options:
 *
 * - "log_nondominated: STRING" determines how the nondominated solutions are handled. STRING can take on the
 * values "none" (don't log nondominated solutions), "final" (log only the final nondominated solutions),
 * "all" (log every solution that is nondominated at creation time) and "read" (the nondominated solutions
 * are not logged, but are passed to the logger as input - this is a functionality needed in pre-processing
 * of the data). The default value is "all".
 *
 * - "log_decision_variables: STRING" determines whether the decision variables are to be logged in addition
 * to the objective variables in the output of nondominated solutions. STRING can take on the values "none"
 * (don't output decision variables), "low_dim"(output decision variables only for dimensions lower or equal
 * to 5) and "all" (output all decision variables). The default value is "low_dim".
 *
 * - "compute_indicators: VALUE" determines whether to compute and output performance indicators (1) or not
 * (0). The default value is 1.
 *
 * - "produce_all_data: VALUE" determines whether to produce all data required for the workshop. If set to 1,
 * it overwrites some other options and is equivalent to setting "log_nondominated: all",
 * "log_decision_variables: low_dim" and "compute_indicators: 1". If set to 0, it does not change the values
 * of the other options. The default value is 0.
 */
static void observer_biobj(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys) {

  observer_biobj_data_t *observer_biobj;
  char string_value[COCO_PATH_MAX];

  /* Sets the valid keys for bbob-biobj observer options
   * IMPORTANT: This list should be up-to-date with the code and the documentation */
  const char *known_keys[] = { "log_nondominated", "log_decision_variables", "compute_indicators",
      "produce_all_data" };
  *option_keys = coco_option_keys_allocate(sizeof(known_keys) / sizeof(char *), known_keys);

  observer_biobj = (observer_biobj_data_t *) coco_allocate_memory(sizeof(*observer_biobj));

  observer_biobj->log_nondom_mode = LOG_NONDOM_ALL;
  if (coco_options_read_string(options, "log_nondominated", string_value) > 0) {
    if (strcmp(string_value, "none") == 0)
      observer_biobj->log_nondom_mode = LOG_NONDOM_NONE;
    else if (strcmp(string_value, "final") == 0)
      observer_biobj->log_nondom_mode = LOG_NONDOM_FINAL;
    else if (strcmp(string_value, "all") == 0)
      observer_biobj->log_nondom_mode = LOG_NONDOM_ALL;
    else if (strcmp(string_value, "read") == 0)
      observer_biobj->log_nondom_mode = LOG_NONDOM_READ;
  }

  observer_biobj->log_vars_mode = LOG_VARS_LOW_DIM;
  if (coco_options_read_string(options, "log_decision_variables", string_value) > 0) {
    if (strcmp(string_value, "none") == 0)
      observer_biobj->log_vars_mode = LOG_VARS_NEVER;
    else if (strcmp(string_value, "all") == 0)
      observer_biobj->log_vars_mode = LOG_VARS_ALWAYS;
    else if (strcmp(string_value, "low_dim") == 0)
      observer_biobj->log_vars_mode = LOG_VARS_LOW_DIM;
  }

  if (coco_options_read_int(options, "compute_indicators", &(observer_biobj->compute_indicators)) == 0)
    observer_biobj->compute_indicators = 1;

  if (coco_options_read_int(options, "produce_all_data", &(observer_biobj->produce_all_data)) == 0)
    observer_biobj->produce_all_data = 0;

  if (observer_biobj->produce_all_data) {
    observer_biobj->compute_indicators = 1;
    observer_biobj->log_nondom_mode = LOG_NONDOM_ALL;
    observer_biobj->log_vars_mode = LOG_VARS_LOW_DIM;
  }

  if (observer_biobj->compute_indicators) {
    observer_biobj->previous_function = -1;
    observer_biobj->previous_dimension = -1;
  }

  observer->logger_allocate_function = logger_biobj;
  observer->logger_free_function = logger_biobj_free;
  observer->data_free_function = NULL;
  observer->data = observer_biobj;

  if ((observer_biobj->log_nondom_mode == LOG_NONDOM_NONE) && (!observer_biobj->compute_indicators)) {
    /* No logging required */
    observer->is_active = 0;
  }
}
#line 42 "code-experiments/src/logger_biobj.c"

#line 44 "code-experiments/src/logger_biobj.c"

/** @brief Number of implemented indicators */
#define LOGGER_BIOBJ_NUMBER_OF_INDICATORS 1

/** @brief Names of implemented indicators
 *
 * "hyp" stands for the hypervolume indicator.
 * */
const char *logger_biobj_indicators[LOGGER_BIOBJ_NUMBER_OF_INDICATORS] = { "hyp" };

/**
 * @brief The indicator type.
 *
 * <B> The hypervolume indicator ("hyp") </B>
 *
 * The hypervolume indicator measures the volume of the portion of the ROI in the objective space that is
 * dominated by the current Pareto front approximation. Instead of logging the hypervolume indicator value,
 * this implementation logs the difference between the best know hypervolume indicator (a value stored in
 * best_value) and the hypervolume indicator of the current Pareto front approximation (current_value). The
 * current_value equals 0 if no solution is located in the ROI. In order to be able to register the
 * performance of an optimizer even before the ROI is reached, an additional value is computed when no
 * solutions are located inside the ROI. This value is stored in additional_penalty and equals the
 * normalized distance to the ROI of the solution closest to the ROI (additional_penalty is set to 0 as
 * soon as a solution reaches the ROI). The final value to be logged (overall_value) is therefore computed
 * in the following way:
 *
 * overall_value = best_value - current_value + additional_penalty
 *
 * @note Other indicators are yet to be implemented.
 */
typedef struct {

  char *name;                /**< @brief Name of the indicator used for identification and the output. */

  FILE *dat_file;            /**< @brief File for logging indicator values at predefined values. */
  FILE *tdat_file;           /**< @brief File for logging indicator values at predefined evaluations. */
  FILE *info_file;           /**< @brief File for logging summary information on algorithm performance. */

  int target_hit;            /**< @brief Whether the target was hit in the latest evaluation. */
  coco_observer_targets_t *targets;
                             /**< @brief Triggers based on target values. */
  int evaluation_logged;     /**< @brief Whether the whether the latest evaluation was logged. */
  coco_observer_evaluations_t *evaluations;
                             /**< @brief Triggers based on numbers of evaluations. */

  double best_value;         /**< @brief The best known indicator value for this problem. */
  double current_value;      /**< @brief The current indicator value. */
  double additional_penalty; /**< @brief Additional penalty for solutions outside the ROI. */
  double overall_value;      /**< @brief The overall value of the indicator tested for target hits. */
  double previous_value;     /**< @brief The previous overall value of the indicator. */

} logger_biobj_indicator_t;

/**
 * @brief The bi-objective logger data type.
 *
 * @note Some fields from the observers (coco_observer as well as observer_biobj) need to be copied here
 * because the observers can be deleted before the logger is finalized and we need these fields for
 * finalization.
 */
typedef struct {
  observer_biobj_log_nondom_e log_nondom_mode;
                                 /**< @brief Mode for archiving nondominated solutions. */
  FILE *adat_file;               /**< @brief File for archiving nondominated solutions (all or final). */

  int log_vars;                  /**< @brief Whether to log the decision values. */

  int precision_x;               /**< @brief Precision for outputting decision values. */
  int precision_f;               /**< @brief Precision for outputting objective values. */

  size_t number_of_evaluations;  /**< @brief The number of evaluations performed so far. */
  size_t number_of_variables;    /**< @brief Dimension of the problem. */
  size_t number_of_objectives;   /**< @brief Number of objectives (clearly equal to 2). */
  size_t suite_dep_instance;     /**< @brief Suite-dependent instance number of the observed problem. */

  size_t previous_evaluations;   /**< @brief The number of evaluations from the previous call to the logger. */

  avl_tree_t *archive_tree;      /**< @brief The tree keeping currently non-dominated solutions. */
  avl_tree_t *buffer_tree;       /**< @brief The tree with pointers to nondominated solutions that haven't
                                      been logged yet. */

  /* Indicators (TODO: Implement others!) */
  int compute_indicators;        /**< @brief Whether to compute the indicators. */
  logger_biobj_indicator_t *indicators[LOGGER_BIOBJ_NUMBER_OF_INDICATORS];
                                 /**< @brief The implemented indicators. */
} logger_biobj_data_t;

/**
 * @brief The type for the node's item in the AVL tree as used by the bi-objective logger.
 *
 * Contains information on the exact objective values (y) and their rounded normalized values (normalized_y).
 * The exact values are used for output, while archive update and indicator computation use the normalized
 * values.
 */
typedef struct {
  double *x;                 /**< @brief The decision values of this solution. */
  double *y;                 /**< @brief The values of objectives of this solution. */
  double *normalized_y;      /**< @brief The values of normalized objectives of this solution. */
  size_t evaluation_number;  /**< @brief The evaluation number of when the solution was created. */

  double indicator_contribution[LOGGER_BIOBJ_NUMBER_OF_INDICATORS];
                      /**< @brief The contribution of this solution to the overall indicator values. */
  int within_ROI;     /**< @brief Whether the solution is within the region of interest (ROI). */

} logger_biobj_avl_item_t;

/**
 * @brief Creates and returns the information on the solution in the form of a node's item in the AVL tree.
 */
static logger_biobj_avl_item_t* logger_biobj_node_create(const coco_problem_t *problem,
                                                         const double *x,
                                                         const double *y,
                                                         const size_t evaluation_number,
                                                         const size_t dim,
                                                         const size_t num_obj) {

  size_t i;

  /* Allocate memory to hold the data structure logger_biobj_node_t */
  logger_biobj_avl_item_t *item = (logger_biobj_avl_item_t*) coco_allocate_memory(sizeof(*item));

  /* Allocate memory to store the (copied) data of the new node */
  item->x = coco_allocate_vector(dim);
  item->y = coco_allocate_vector(num_obj);

  /* Copy the data */
  for (i = 0; i < dim; i++)
    item->x[i] = x[i];
  for (i = 0; i < num_obj; i++)
    item->y[i] = y[i];

  /* Compute the normalized y */
  item->normalized_y = mo_normalize(item->y, problem->best_value, problem->nadir_value, num_obj);
  item->within_ROI = mo_is_within_ROI(item->normalized_y, num_obj);

  item->evaluation_number = evaluation_number;
  for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++)
    item->indicator_contribution[i] = 0;

  return item;
}

/**
 * @brief Frees the data of the given logger_biobj_avl_item_t.
 */
static void logger_biobj_node_free(logger_biobj_avl_item_t *item, void *userdata) {

  coco_free_memory(item->x);
  coco_free_memory(item->y);
  coco_free_memory(item->normalized_y);
  coco_free_memory(item);
  (void) userdata; /* To silence the compiler */
}

/**
 * @brief Defines the ordering of AVL tree nodes based on the value of the last objective.
 *
 * @note This ordering is used by the archive_tree.
 */
static int avl_tree_compare_by_last_objective(const logger_biobj_avl_item_t *item1,
                                              const logger_biobj_avl_item_t *item2,
                                              void *userdata) {
  if (coco_double_almost_equal(item1->normalized_y[1], item2->normalized_y[1], mo_precision))
    return 0;
  else if (item1->normalized_y[1] < item2->normalized_y[1])
    return -1;
  else
    return 1;

  (void) userdata; /* To silence the compiler */
}

/**
 * @brief Defines the ordering of AVL tree nodes based on the evaluation number (the time when the nodes were
 * created).
 *
 * @note This ordering is used by the buffer_tree.
 */
static int avl_tree_compare_by_eval_number(const logger_biobj_avl_item_t *item1,
                                           const logger_biobj_avl_item_t *item2,
                                           void *userdata) {
  if (item1->evaluation_number < item2->evaluation_number)
    return -1;
  else if (item1->evaluation_number > item2->evaluation_number)
    return 1;
  else
    return 0;

  (void) userdata; /* To silence the compiler */
}

/**
 * @brief Outputs the AVL tree to the given file. Returns the number of nodes in the tree.
 */
static size_t logger_biobj_tree_output(FILE *file,
                                       const avl_tree_t *tree,
                                       const size_t dim,
                                       const size_t num_obj,
                                       const int log_vars,
                                       const int precision_x,
                                       const int precision_f) {

  avl_node_t *solution;
  size_t i;
  size_t j;
  size_t number_of_nodes = 0;

  if (tree->tail) {
    /* There is at least a solution in the tree to output */
    solution = tree->head;
    while (solution != NULL) {
      fprintf(file, "%lu\t", (unsigned long) ((logger_biobj_avl_item_t*) solution->item)->evaluation_number);
      for (j = 0; j < num_obj; j++)
        fprintf(file, "%.*e\t", precision_f, ((logger_biobj_avl_item_t*) solution->item)->y[j]);
      if (log_vars) {
        for (i = 0; i < dim; i++)
          fprintf(file, "%.*e\t", precision_x, ((logger_biobj_avl_item_t*) solution->item)->x[i]);
      }
      fprintf(file, "\n");
      solution = solution->next;
      number_of_nodes++;
    }
  }

  return number_of_nodes;
}

/**
 * @brief Updates the archive and buffer trees with the given node.
 *
 * Checks for domination and updates the archive tree and the values of the indicators if the given node is
 * not weakly dominated by existing nodes in the archive tree. This is where the main computation of
 * indicator values takes place.
 *
 * @return 1 if the update was performed and 0 otherwise.
 */
static int logger_biobj_tree_update(logger_biobj_data_t *logger,
                                    logger_biobj_avl_item_t *node_item) {

  avl_node_t *node, *next_node, *new_node;
  int trigger_update = 0;
  int dominance;
  size_t i;
  int previous_unavailable = 0;

  /* Find the first point that is not worse than the new point (NULL if such point does not exist) */
  node = avl_item_search_right(logger->archive_tree, node_item, NULL);

  if (node == NULL) {
    /* The new point is an extreme point */
    trigger_update = 1;
    next_node = logger->archive_tree->head;
  } else {
    dominance = mo_get_dominance(node_item->normalized_y,
        ((logger_biobj_avl_item_t*) node->item)->normalized_y, logger->number_of_objectives);
    if (dominance > -1) {
      trigger_update = 1;
      next_node = node->next;
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        if (logger->compute_indicators) {
          for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            logger->indicators[i]->current_value -= ((logger_biobj_avl_item_t*) node->item)->indicator_contribution[i];
          }
        }
        avl_item_delete(logger->buffer_tree, node->item);
        avl_node_delete(logger->archive_tree, node);
      }
    } else {
      /* The new point is dominated or equal to an existing one, nothing more to do */
      trigger_update = 0;
    }
  }

  if (!trigger_update) {
    logger_biobj_node_free(node_item, NULL);
  } else {
    /* Perform tree update */
    while (next_node != NULL) {
      /* Check the dominance relation between the new node and the next node. There are only two possibilities:
       * dominance = 0: the new node and the next node are nondominated
       * dominance = 1: the new node dominates the next node */
      node = next_node;
      dominance = mo_get_dominance(node_item->normalized_y,
          ((logger_biobj_avl_item_t*) node->item)->normalized_y, logger->number_of_objectives);
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        if (logger->compute_indicators) {
          for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            logger->indicators[i]->current_value -= ((logger_biobj_avl_item_t*) node->item)->indicator_contribution[i];
          }
        }
        next_node = node->next;
        avl_item_delete(logger->buffer_tree, node->item);
        avl_node_delete(logger->archive_tree, node);
      } else {
        break;
      }
    }

    new_node = avl_item_insert(logger->archive_tree, node_item);
    assert(new_node != NULL);
    avl_item_insert(logger->buffer_tree, node_item);

    if (logger->compute_indicators) {
      if (node_item->within_ROI) {
        /* Compute indicator value for new node and update the indicator value of the affected nodes */
        logger_biobj_avl_item_t *next_item, *previous_item;

        if (new_node->next != NULL) {
          next_item = (logger_biobj_avl_item_t*) new_node->next->item;
          if (next_item->within_ROI) {
            for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
              logger->indicators[i]->current_value -= next_item->indicator_contribution[i];
              if (strcmp(logger->indicators[i]->name, "hyp") == 0) {
                next_item->indicator_contribution[i] = (node_item->normalized_y[0] - next_item->normalized_y[0])
                    * (1 - next_item->normalized_y[1]);
                assert(next_item->indicator_contribution[i] > 0);
              } else {
                coco_error(
                    "logger_biobj_tree_update(): Indicator computation not implemented yet for indicator %s",
                    logger->indicators[i]->name);
              }
              logger->indicators[i]->current_value += next_item->indicator_contribution[i];
            }
          }
        }

        previous_unavailable = 0;
        if (new_node->prev != NULL) {
          previous_item = (logger_biobj_avl_item_t*) new_node->prev->item;
          if (previous_item->within_ROI) {
            for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
              if (strcmp(logger->indicators[i]->name, "hyp") == 0) {
                node_item->indicator_contribution[i] = (previous_item->normalized_y[0] - node_item->normalized_y[0])
                    * (1 - node_item->normalized_y[1]);
                assert(node_item->indicator_contribution[i] > 0);
              } else {
                coco_error(
                    "logger_biobj_tree_update(): Indicator computation not implemented yet for indicator %s",
                    logger->indicators[i]->name);
              }
            }
          } else {
            previous_unavailable = 1;
          }
        } else {
          previous_unavailable = 1;
        }

        if (previous_unavailable) {
          /* Previous item does not exist or is out of ROI, use reference point instead */
          for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            if (strcmp(logger->indicators[i]->name, "hyp") == 0) {
              node_item->indicator_contribution[i] = (1 - node_item->normalized_y[0])
                  * (1 - node_item->normalized_y[1]);
              assert(node_item->indicator_contribution[i] > 0);
            } else {
              coco_error(
                  "logger_biobj_tree_update(): Indicator computation not implemented yet for indicator %s",
                  logger->indicators[i]->name);
            }
          }
        }

        for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
          if (strcmp(logger->indicators[i]->name, "hyp") == 0) {
            assert(node_item->indicator_contribution[i] >= 0);
            logger->indicators[i]->current_value += node_item->indicator_contribution[i];
          }
        }
      }
    }
  }

  return trigger_update;
}

/**
 * @brief Initializes the indicator with name indicator_name.
 *
 * Opens files for writing and resets counters.
 */
static logger_biobj_indicator_t *logger_biobj_indicator(const logger_biobj_data_t *logger,
                                                        const coco_observer_t *observer,
                                                        const coco_problem_t *problem,
                                                        const char *indicator_name) {

  observer_biobj_data_t *observer_biobj;
  logger_biobj_indicator_t *indicator;
  char *prefix, *file_name, *path_name;
  int info_file_exists = 0;

  indicator = (logger_biobj_indicator_t *) coco_allocate_memory(sizeof(*indicator));
  assert(observer);
  assert(observer->data);
  observer_biobj = (observer_biobj_data_t *) observer->data;

  indicator->name = coco_strdup(indicator_name);

  indicator->best_value = suite_biobj_get_best_value(indicator->name, problem->problem_id);
  indicator->target_hit = 0;
  indicator->evaluation_logged = 0;
  indicator->current_value = 0;
  indicator->additional_penalty = DBL_MAX;
  indicator->overall_value = 0;
  indicator->previous_value = 0;

  indicator->targets = coco_observer_targets(observer->number_target_triggers, observer->target_precision);
  indicator->evaluations = coco_observer_evaluations(observer->base_evaluation_triggers, problem->number_of_variables);

  /* Prepare the info file */
  path_name = coco_allocate_string(COCO_PATH_MAX);
  memcpy(path_name, observer->result_folder, strlen(observer->result_folder) + 1);
  coco_create_directory(path_name);
  file_name = coco_strdupf("%s_%s.info", problem->problem_type, indicator_name);
  coco_join_path(path_name, COCO_PATH_MAX, file_name, NULL);
  info_file_exists = coco_file_exists(path_name);
  indicator->info_file = fopen(path_name, "a");
  if (indicator->info_file == NULL) {
    coco_error("logger_biobj_indicator() failed to open file '%s'.", path_name);
    return NULL; /* Never reached */
  }
  coco_free_memory(file_name);
  coco_free_memory(path_name);

  /* Prepare the tdat file */
  path_name = coco_allocate_string(COCO_PATH_MAX);
  memcpy(path_name, observer->result_folder, strlen(observer->result_folder) + 1);
  coco_join_path(path_name, COCO_PATH_MAX, problem->problem_type, NULL);
  coco_create_directory(path_name);
  prefix = coco_remove_from_string(problem->problem_id, "_i", "_d");
  file_name = coco_strdupf("%s_%s.tdat", prefix, indicator_name);
  coco_join_path(path_name, COCO_PATH_MAX, file_name, NULL);
  indicator->tdat_file = fopen(path_name, "a");
  if (indicator->tdat_file == NULL) {
    coco_error("logger_biobj_indicator() failed to open file '%s'.", path_name);
    return NULL; /* Never reached */
  }
  coco_free_memory(file_name);
  coco_free_memory(path_name);

  /* Prepare the dat file */
  path_name = coco_allocate_string(COCO_PATH_MAX);
  memcpy(path_name, observer->result_folder, strlen(observer->result_folder) + 1);
  coco_join_path(path_name, COCO_PATH_MAX, problem->problem_type, NULL);
  coco_create_directory(path_name);
  file_name = coco_strdupf("%s_%s.dat", prefix, indicator_name);
  coco_join_path(path_name, COCO_PATH_MAX, file_name, NULL);
  indicator->dat_file = fopen(path_name, "a");
  if (indicator->dat_file == NULL) {
    coco_error("logger_biobj_indicator() failed to open file '%s'.", path_name);
    return NULL; /* Never reached */
  }

  /* Output header information to the info file */
  if (!info_file_exists) {
    /* Output algorithm name */
    assert(problem->suite);
    /* TODO: Use this once suite can be read by the postprocessing
    fprintf(indicator->info_file,
        "suite = '%s', algorithm = '%s', indicator = '%s', folder = '%s', coco_version = '%s'\n%% %s",
        problem->suite->suite_name, observer->algorithm_name, indicator_name, problem->problem_type,
        coco_version, observer->algorithm_info);*/
    fprintf(indicator->info_file,
        "algorithm = '%s', indicator = '%s', folder = '%s', coco_version = '%s'\n%% %s",
        observer->algorithm_name, indicator_name, problem->problem_type,
        coco_version, observer->algorithm_info);
    if (logger->log_nondom_mode == LOG_NONDOM_READ)
      fprintf(indicator->info_file, " (reconstructed)");
  }
  if ((observer_biobj->previous_function != problem->suite_dep_function)
    || (observer_biobj->previous_dimension != problem->number_of_variables)) {
    fprintf(indicator->info_file, "\nfunction = %2lu, ", (unsigned long) problem->suite_dep_function);
    fprintf(indicator->info_file, "dim = %2lu, ", (unsigned long) problem->number_of_variables);
    fprintf(indicator->info_file, "%s", file_name);
  }

  coco_free_memory(prefix);
  coco_free_memory(file_name);
  coco_free_memory(path_name);

  /* Output header information to the dat file */
  fprintf(indicator->dat_file, "%%\n%% index = %lu, name = %s\n", (unsigned long) problem->suite_dep_index,
      problem->problem_name);
  fprintf(indicator->dat_file, "%% instance = %lu, reference value = %.*e\n",
      (unsigned long) problem->suite_dep_instance, logger->precision_f, indicator->best_value);
  fprintf(indicator->dat_file, "%% function evaluation | indicator value | target hit\n");

  /* Output header information to the tdat file */
  fprintf(indicator->tdat_file, "%%\n%% index = %lu, name = %s\n", (unsigned long) problem->suite_dep_index,
      problem->problem_name);
  fprintf(indicator->tdat_file, "%% instance = %lu, reference value = %.*e\n",
      (unsigned long) problem->suite_dep_instance, logger->precision_f, indicator->best_value);
  fprintf(indicator->tdat_file, "%% function evaluation | indicator value\n");

  return indicator;
}

/**
 * @brief Outputs the final information about this indicator.
 */
static void logger_biobj_indicator_finalize(logger_biobj_indicator_t *indicator, const logger_biobj_data_t *logger) {

  /* Log the last eval_number in the dat file if wasn't already logged */
  if (!indicator->target_hit) {
    fprintf(indicator->dat_file, "%lu\t%.*e\t%.*e\n", (unsigned long) logger->number_of_evaluations,
        logger->precision_f, indicator->overall_value, logger->precision_f,
        ((coco_observer_targets_t *) indicator->targets)->value);
  }

  /* Log the last eval_number in the tdat file if wasn't already logged */
  if (!indicator->evaluation_logged) {
    fprintf(indicator->tdat_file, "%lu\t%.*e\n", (unsigned long) logger->number_of_evaluations,
        logger->precision_f, indicator->overall_value);
  }

  /* Log the information in the info file */
  fprintf(indicator->info_file, ", %lu:%lu|%.1e", (unsigned long) logger->suite_dep_instance,
      (unsigned long) logger->number_of_evaluations, indicator->overall_value);
  fflush(indicator->info_file);
}

/**
 * @brief Frees the memory of the given indicator.
 */
static void logger_biobj_indicator_free(void *stuff) {

  logger_biobj_indicator_t *indicator;

  assert(stuff != NULL);
  indicator = (logger_biobj_indicator_t *) stuff;

  if (indicator->name != NULL) {
    coco_free_memory(indicator->name);
    indicator->name = NULL;
  }

  if (indicator->dat_file != NULL) {
    fclose(indicator->dat_file);
    indicator->dat_file = NULL;
  }

  if (indicator->tdat_file != NULL) {
    fclose(indicator->tdat_file);
    indicator->tdat_file = NULL;
  }

  if (indicator->info_file != NULL) {
    fclose(indicator->info_file);
    indicator->info_file = NULL;
  }

  if (indicator->targets != NULL){
    coco_free_memory(indicator->targets);
    indicator->targets = NULL;
  }

  if (indicator->evaluations != NULL){
    coco_observer_evaluations_free(indicator->evaluations);
    indicator->evaluations = NULL;
  }

  coco_free_memory(stuff);

}

/*
 * @brief Outputs the information according to the observer options.
 *
 * Outputs to the:
 * - dat file, if the archive was updated and a new target was reached for an indicator;
 * - tdat file, if the number of evaluations matches one of the predefined numbers.
 *
 * Note that a target is reached when
 * best_value - current_value + additional_penalty <= relative_target_value
 *
 * The relative_target_value is a target for indicator difference, not the actual indicator value!
 */
static void logger_biobj_output(logger_biobj_data_t *logger,
                                const int update_performed,
                                const logger_biobj_avl_item_t *node_item) {

  size_t i, j;
  logger_biobj_indicator_t *indicator;

  if (logger->compute_indicators) {
    for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {

      indicator = logger->indicators[i];
      indicator->target_hit = 0;
      indicator->previous_value = indicator->overall_value;

      /* If the update was performed, update the overall indicator value */
      if (update_performed) {
        /* Compute the overall_value of an indicator */
        if (strcmp(indicator->name, "hyp") == 0) {
          if (coco_double_almost_equal(indicator->current_value, 0, mo_precision)) {
            /* Update the additional penalty for hypervolume (the minimal distance from the nondominated set
             * to the ROI) */
            double new_distance = mo_get_distance_to_ROI(node_item->normalized_y, logger->number_of_objectives);
            indicator->additional_penalty = coco_double_min(indicator->additional_penalty, new_distance);
            assert(indicator->additional_penalty >= 0);
          } else {
            indicator->additional_penalty = 0;
          }
          indicator->overall_value = indicator->best_value - indicator->current_value
              + indicator->additional_penalty;
        } else {
          coco_error("logger_biobj_evaluate(): Indicator computation not implemented yet for indicator %s",
              indicator->name);
        }

        /* Check whether a target was hit */
        indicator->target_hit = coco_observer_targets_trigger(indicator->targets, indicator->overall_value);
      }

      /* Log to the dat file if a target was hit */
      if (indicator->target_hit) {
        fprintf(indicator->dat_file, "%lu\t%.*e\t%.*e\n", (unsigned long) logger->number_of_evaluations,
            logger->precision_f, indicator->overall_value, logger->precision_f,
            ((coco_observer_targets_t *) indicator->targets)->value);
      }

      if (logger->log_nondom_mode == LOG_NONDOM_READ) {
        /* Log to the tdat file the previous indicator value if any evaluation number between the previous and
         * this one matches one of the predefined evaluation numbers. */
        for (j = logger->previous_evaluations + 1; j < logger->number_of_evaluations; j++) {
          indicator->evaluation_logged = coco_observer_evaluations_trigger(indicator->evaluations, j);
          if (indicator->evaluation_logged) {
            fprintf(indicator->tdat_file, "%lu\t%.*e\n", (unsigned long) j, logger->precision_f,
                indicator->previous_value);
          }
        }
      }

      /* Log to the tdat file if the number of evaluations matches one of the predefined numbers */
      indicator->evaluation_logged = coco_observer_evaluations_trigger(indicator->evaluations,
          logger->number_of_evaluations);
      if (indicator->evaluation_logged) {
        fprintf(indicator->tdat_file, "%lu\t%.*e\n", (unsigned long) logger->number_of_evaluations,
            logger->precision_f, indicator->overall_value);
      }

    }
  }
}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information according to
 * observer options.
 */
static void logger_biobj_evaluate(coco_problem_t *problem, const double *x, double *y) {

  logger_biobj_data_t *logger;
  logger_biobj_avl_item_t *node_item;
  int update_performed;
  coco_problem_t *inner_problem;

  logger = (logger_biobj_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  /* Evaluate function */
  coco_evaluate_function(inner_problem, x, y);
  logger->number_of_evaluations++;

  node_item = logger_biobj_node_create(inner_problem, x, y, logger->number_of_evaluations, logger->number_of_variables,
      logger->number_of_objectives);

  /* Update the archive with the new solution, if it is not dominated by or equal to existing solutions in
   * the archive */
  update_performed = logger_biobj_tree_update(logger, node_item);

  /* If the archive was updated and you need to log all nondominated solutions, output the new solution to
   * nondom_file */
  if (update_performed && (logger->log_nondom_mode == LOG_NONDOM_ALL)) {
    logger_biobj_tree_output(logger->adat_file, logger->buffer_tree, logger->number_of_variables,
        logger->number_of_objectives, logger->log_vars, logger->precision_x, logger->precision_f);
    avl_tree_purge(logger->buffer_tree);

    /* Flush output so that impatient users can see progress. */
    fflush(logger->adat_file);
  }

  /* Output according to observer options */
  logger_biobj_output(logger, update_performed, node_item);
}

/**
 * Sets the number of evaluations, adds the objective vector to the archive and outputs information according
 * to observer options (but does not output the archive).
 *
 * @note Vector y must point to a correctly sized allocated memory region and the given evaluation number must
 * be larger than the existing one.
 *
 * @param problem The given COCO problem.
 * @param evaluation The number of evaluations.
 * @param y The objective vector.
 * @return 1 if archive was updated was done and 0 otherwise.
 */
int coco_logger_biobj_feed_solution(coco_problem_t *problem, const size_t evaluation, const double *y) {

  logger_biobj_data_t *logger;
  logger_biobj_avl_item_t *node_item;
  int update_performed;
  coco_problem_t *inner_problem;
  double *x;
  size_t i;

  assert(problem != NULL);
  logger = (logger_biobj_data_t *) coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);
  assert(logger->log_nondom_mode == LOG_NONDOM_READ);

  /* Set the number of evaluations */
  logger->previous_evaluations = logger->number_of_evaluations;
  if (logger->previous_evaluations >= evaluation)
    coco_error("coco_logger_biobj_reconstruct(): Evaluation %lu came before evaluation %lu. Note that "
        "the evaluations need to be always increasing.", logger->previous_evaluations, evaluation);
  logger->number_of_evaluations = evaluation;

  /* Update the archive with the new solution */
  x = coco_allocate_vector(problem->number_of_variables);
  for (i = 0; i < problem->number_of_variables; i++)
    x[i] = 0;
  node_item = logger_biobj_node_create(inner_problem, x, y, logger->number_of_evaluations,
      logger->number_of_variables, logger->number_of_objectives);
  coco_free_memory(x);

  /* Update the archive */
  update_performed = logger_biobj_tree_update(logger, node_item);

  /* Output according to observer options */
  logger_biobj_output(logger, update_performed, node_item);

  return update_performed;
}

/**
 * @brief Outputs the final nondominated solutions to the archive file.
 */
static void logger_biobj_finalize(logger_biobj_data_t *logger) {

  avl_tree_t *resorted_tree;
  avl_node_t *solution;

  /* Re-sort archive_tree according to time stamp and then output it */
  resorted_tree = avl_tree_construct((avl_compare_t) avl_tree_compare_by_eval_number, NULL);

  if (logger->archive_tree->tail) {
    /* There is at least a solution in the tree to output */
    solution = logger->archive_tree->head;
    while (solution != NULL) {
      avl_item_insert(resorted_tree, solution->item);
      solution = solution->next;
    }
  }

  logger_biobj_tree_output(logger->adat_file, resorted_tree, logger->number_of_variables,
      logger->number_of_objectives, logger->log_vars, logger->precision_x, logger->precision_f);

  avl_tree_destruct(resorted_tree);
}

/**
 * @brief Frees the memory of the given biobjective logger.
 */
static void logger_biobj_free(void *stuff) {

  logger_biobj_data_t *logger;
  size_t i;

  assert(stuff != NULL);
  logger = (logger_biobj_data_t *) stuff;

  if (logger->log_nondom_mode == LOG_NONDOM_FINAL) {
     logger_biobj_finalize(logger);
  }

  if (logger->compute_indicators) {
    for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
      logger_biobj_indicator_finalize(logger->indicators[i], logger);
      logger_biobj_indicator_free(logger->indicators[i]);
    }
  }

  if (((logger->log_nondom_mode == LOG_NONDOM_ALL) || (logger->log_nondom_mode == LOG_NONDOM_FINAL)) &&
      (logger->adat_file != NULL)) {
    fprintf(logger->adat_file, "%% evaluations = %lu\n", (unsigned long) logger->number_of_evaluations);
    fclose(logger->adat_file);
    logger->adat_file = NULL;
  }

  avl_tree_destruct(logger->archive_tree);
  avl_tree_destruct(logger->buffer_tree);

}

/**
 * @brief Initializes the biobjective logger.
 *
 * Copies all observer field values that are needed after initialization into logger field values for two
 * reasons:
 * - If the observer is deleted before the suite, the observer is not available anymore when the logger
 * is finalized.
 * - This reduces function calls.
 */
static coco_problem_t *logger_biobj(coco_observer_t *observer, coco_problem_t *inner_problem) {

  coco_problem_t *problem;
  logger_biobj_data_t *logger_biobj;
  observer_biobj_data_t *observer_biobj;
  const char nondom_folder_name[] = "archive";
  char *path_name, *file_name = NULL;
  size_t i;

  if (inner_problem->number_of_objectives != 2) {
    coco_error("logger_biobj(): The bi-objective logger cannot log a problem with %d objective(s)",
        inner_problem->number_of_objectives);
    return NULL; /* Never reached. */
  }

  logger_biobj = (logger_biobj_data_t *) coco_allocate_memory(sizeof(*logger_biobj));

  logger_biobj->number_of_evaluations = 0;
  logger_biobj->previous_evaluations = 0;
  logger_biobj->number_of_variables = inner_problem->number_of_variables;
  logger_biobj->number_of_objectives = inner_problem->number_of_objectives;
  logger_biobj->suite_dep_instance = inner_problem->suite_dep_instance;

  observer_biobj = (observer_biobj_data_t *) observer->data;
  /* Copy values from the observes that you might need even if they do not exist any more */
  logger_biobj->log_nondom_mode = observer_biobj->log_nondom_mode;
  logger_biobj->compute_indicators = observer_biobj->compute_indicators;
  logger_biobj->precision_x = observer->precision_x;
  logger_biobj->precision_f = observer->precision_f;

  if (((observer_biobj->log_vars_mode == LOG_VARS_LOW_DIM) && (inner_problem->number_of_variables > 5))
      || (observer_biobj->log_vars_mode == LOG_VARS_NEVER))
    logger_biobj->log_vars = 0;
  else
    logger_biobj->log_vars = 1;

  /* Initialize logging of nondominated solutions into the archive file */
  if ((logger_biobj->log_nondom_mode == LOG_NONDOM_ALL) ||
      (logger_biobj->log_nondom_mode == LOG_NONDOM_FINAL)) {

    /* Create the path to the file */
    path_name = coco_allocate_string(COCO_PATH_MAX);
    memcpy(path_name, observer->result_folder, strlen(observer->result_folder) + 1);
    coco_join_path(path_name, COCO_PATH_MAX, nondom_folder_name, NULL);
    coco_create_directory(path_name);

    /* Construct file name */
    if (logger_biobj->log_nondom_mode == LOG_NONDOM_ALL)
      file_name = coco_strdupf("%s_nondom_all.adat", inner_problem->problem_id);
    else if (logger_biobj->log_nondom_mode == LOG_NONDOM_FINAL)
      file_name = coco_strdupf("%s_nondom_final.adat", inner_problem->problem_id);
    coco_join_path(path_name, COCO_PATH_MAX, file_name, NULL);
    coco_free_memory(file_name);

    /* Open and initialize the archive file */
    logger_biobj->adat_file = fopen(path_name, "a");
    if (logger_biobj->adat_file == NULL) {
      coco_error("logger_biobj() failed to open file '%s'.", path_name);
      return NULL; /* Never reached */
    }
    coco_free_memory(path_name);

    /* Output header information */
    fprintf(logger_biobj->adat_file, "%% instance = %lu, name = %s\n",
        (unsigned long) inner_problem->suite_dep_instance, inner_problem->problem_name);
    if (logger_biobj->log_vars) {
      fprintf(logger_biobj->adat_file, "%% function evaluation | %lu objectives | %lu variables\n",
          (unsigned long) inner_problem->number_of_objectives,
          (unsigned long) inner_problem->number_of_variables);
    } else {
      fprintf(logger_biobj->adat_file, "%% function evaluation | %lu objectives \n",
          (unsigned long) inner_problem->number_of_objectives);
    }
  }

  /* Initialize the AVL trees */
  logger_biobj->archive_tree = avl_tree_construct((avl_compare_t) avl_tree_compare_by_last_objective,
      (avl_free_t) logger_biobj_node_free);
  logger_biobj->buffer_tree = avl_tree_construct((avl_compare_t) avl_tree_compare_by_eval_number, NULL);

  /* Initialize the indicators */
  if (logger_biobj->compute_indicators) {
    for (i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++)
      logger_biobj->indicators[i] = logger_biobj_indicator(logger_biobj, observer, inner_problem, logger_biobj_indicators[i]);

    observer_biobj->previous_function = (long) inner_problem->suite_dep_function;
    observer_biobj->previous_dimension = (long) inner_problem->number_of_variables;
  }

  problem = coco_problem_transformed_allocate(inner_problem, logger_biobj, logger_biobj_free, observer->observer_name);
  problem->evaluate_function = logger_biobj_evaluate;

  return problem;
}
#line 371 "code-experiments/src/coco_observer.c"
#line 1 "code-experiments/src/logger_toy.c"
/**
 * @file logger_toy.c
 * @brief Implementation of the toy logger.
 *
 * Logs the evaluation number, function value the target hit and all the variables each time a target has
 * been hit.
 */

#include <stdio.h>
#include <assert.h>

#line 13 "code-experiments/src/logger_toy.c"

#line 15 "code-experiments/src/logger_toy.c"
#line 16 "code-experiments/src/logger_toy.c"
#line 17 "code-experiments/src/logger_toy.c"
#line 1 "code-experiments/src/observer_toy.c"
/**
 * @file observer_toy.c
 * @brief Implementation of the toy observer.
 */

#line 7 "code-experiments/src/observer_toy.c"
#line 8 "code-experiments/src/observer_toy.c"

static coco_problem_t *logger_toy(coco_observer_t *observer, coco_problem_t *problem);
static void logger_toy_free(void *logger);

/**
 * @brief The toy observer data type.
 */
typedef struct {
  FILE *log_file;            /**< @brief File used for logging. */
} observer_toy_data_t;

/**
 * @brief Frees memory of the toy observer data structure.
 */
static void observer_toy_free(void *stuff) {

  observer_toy_data_t *data;

  assert(stuff != NULL);
  data = (observer_toy_data_t *) stuff;

  if (data->log_file != NULL) {
    fclose(data->log_file);
    data->log_file = NULL;
  }

}

/**
 * @brief Initializes the toy observer.
 *
 * Possible options:
 * - file_name: string (name of the output file; default value is "first_hitting_times.dat")
 */
static void observer_toy(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys) {

  observer_toy_data_t *observer_toy;
  char *string_value;
  char *file_name;

  /* Sets the valid keys for toy observer options
   * IMPORTANT: This list should be up-to-date with the code and the documentation */
  const char *known_keys[] = { "file_name" };
  *option_keys = coco_option_keys_allocate(sizeof(known_keys) / sizeof(char *), known_keys);

  observer_toy = (observer_toy_data_t *) coco_allocate_memory(sizeof(*observer_toy));

  /* Read file_name and number_of_targets from the options and use them to initialize the observer */
  string_value = coco_allocate_string(COCO_PATH_MAX);
  if (coco_options_read_string(options, "file_name", string_value) == 0) {
    strcpy(string_value, "first_hitting_times.dat");
  }

  /* Open log_file */
  file_name = coco_allocate_string(COCO_PATH_MAX);
  memcpy(file_name, observer->result_folder, strlen(observer->result_folder) + 1);
  coco_create_directory(file_name);
  coco_join_path(file_name, COCO_PATH_MAX, string_value, NULL);

  observer_toy->log_file = fopen(file_name, "a");
  if (observer_toy->log_file == NULL) {
    coco_error("observer_toy(): failed to open file %s.", file_name);
    return; /* Never reached */
  }

  coco_free_memory(string_value);
  coco_free_memory(file_name);

  observer->logger_allocate_function = logger_toy;
  observer->logger_free_function = logger_toy_free;
  observer->data_free_function = observer_toy_free;
  observer->data = observer_toy;
}
#line 18 "code-experiments/src/logger_toy.c"

/**
 * @brief The toy logger data type.
 */
typedef struct {
  FILE *log_file;                    /**< @brief Pointer to the file already prepared for logging. */
  coco_observer_targets_t *targets;  /**< @brief Triggers based on target values. */
  size_t number_of_evaluations;      /**< @brief The number of evaluations performed so far. */
  int precision_x;                   /**< @brief Precision for outputting decision values. */
  int precision_f;                   /**< @brief Precision for outputting objective values. */
} logger_toy_data_t;

/**
 * @brief Frees the memory of the given toy logger.
 */
static void logger_toy_free(void *stuff) {

  logger_toy_data_t *logger;

  assert(stuff != NULL);
  logger = (logger_toy_data_t *) stuff;

  if (logger->targets != NULL){
    coco_free_memory(logger->targets);
    logger->targets = NULL;
  }

}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information based on the
 * targets that have been hit.
 */
static void logger_toy_evaluate(coco_problem_t *problem, const double *x, double *y) {

  logger_toy_data_t *logger = (logger_toy_data_t *) coco_problem_transformed_get_data(problem);
  size_t i;

  coco_evaluate_function(coco_problem_transformed_get_inner_problem(problem), x, y);
  logger->number_of_evaluations++;

  /* Output the solution when a new target that has been hit */
  if (coco_observer_targets_trigger(logger->targets, y[0])) {
    fprintf(logger->log_file, "%lu\t%.*e\t%.*e", (unsigned long) logger->number_of_evaluations,
    		logger->precision_f, y[0], logger->precision_f, logger->targets->value);
    for (i = 0; i < problem->number_of_variables; i++) {
      fprintf(logger->log_file, "\t%.*e", logger->precision_x, x[i]);
    }
    fprintf(logger->log_file, "\n");
  }

  /* Flush output so that impatient users can see the progress */
  fflush(logger->log_file);
}

/**
 * @brief Initializes the toy logger.
 */
static coco_problem_t *logger_toy(coco_observer_t *observer, coco_problem_t *inner_problem) {

  logger_toy_data_t *logger_toy;
  coco_problem_t *problem;

  if (inner_problem->number_of_objectives != 1) {
    coco_warning("logger_toy(): The toy logger shouldn't be used to log a problem with %d objectives",
        inner_problem->number_of_objectives);
  }

  /* Initialize the logger_toy_data_t object instance */
  logger_toy = (logger_toy_data_t *) coco_allocate_memory(sizeof(*logger_toy));
  logger_toy->number_of_evaluations = 0;
  logger_toy->targets = coco_observer_targets(observer->number_target_triggers, observer->target_precision);
  logger_toy->log_file = ((observer_toy_data_t *) observer->data)->log_file;
  logger_toy->precision_x = observer->precision_x;
  logger_toy->precision_f = observer->precision_f;

  problem = coco_problem_transformed_allocate(inner_problem, logger_toy, logger_toy_free, observer->observer_name);
  problem->evaluate_function = logger_toy_evaluate;

  /* Output initial information */
  assert(coco_problem_get_suite(inner_problem));
  fprintf(logger_toy->log_file, "\n");
  fprintf(logger_toy->log_file, "suite = '%s', problem_id = '%s', problem_name = '%s', coco_version = '%s'\n",
          coco_problem_get_suite(inner_problem)->suite_name, coco_problem_get_id(inner_problem),
          coco_problem_get_name(inner_problem), coco_version);
  fprintf(logger_toy->log_file, "%% evaluation number | function value | target hit | %lu variables \n",
  		(unsigned long) inner_problem->number_of_variables);

  return problem;
}
#line 372 "code-experiments/src/coco_observer.c"

/**
 * Currently, three observers are supported:
 * - "bbob" is the observer for single-objective (both noisy and noiseless) problems with known optima, which
 * creates *.info, *.dat, *.tdat and *.rdat files and logs the distance to the optimum.
 * - "bbob-biobj" is the observer for bi-objective problems, which creates *.info, *.dat and *.tdat files for
 * the given indicators, as well as an archive folder with *.adat files containing nondominated solutions.
 * - "toy" is a simple observer that logs when a target has been hit.
 *
 * @param observer_name A string containing the name of the observer. Currently supported observer names are
 * "bbob", "bbob-biobj", "toy". Strings "no_observer", "" or NULL return NULL.
 * @param observer_options A string of pairs "key: value" used to pass the options to the observer. Some
 * observer options are general, while others are specific to some observers. Here we list only the general
 * options, see observer_bbob, observer_biobj and observer_toy for options of the specific observers.
 * - "result_folder: NAME" determines the folder within the "exdata" folder into which the results will be
 * output. If the folder with the given name already exists, first NAME_001 will be tried, then NAME_002 and
 * so on. The default value is "default".
 * - "algorithm_name: NAME", where NAME is a short name of the algorithm that will be used in plots (no
 * spaces are allowed). The default value is "ALG".
 * - "algorithm_info: STRING" stores the description of the algorithm. If it contains spaces, it must be
 * surrounded by double quotes. The default value is "" (no description).
 * - "number_target_triggers: VALUE" defines the number of targets between each 10**i and 10**(i+1)
 * (equally spaced in the logarithmic scale) that trigger logging. The default value is 100.
 * - "target_precision: VALUE" defines the precision used for targets (there are no targets for
 * abs(values) < target_precision). The default value is 1e-8.
 * - "number_evaluation_triggers: VALUE" defines the number of evaluations to be logged between each 10**i
 * and 10**(i+1). The default value is 20.
 * - "base_evaluation_triggers: VALUES" defines the base evaluations used to produce an additional
 * evaluation-based logging. The numbers of evaluations that trigger logging are every
 * base_evaluation * dimension * (10**i). For example, if base_evaluation_triggers = "1,2,5", the logger will
 * be triggered by evaluations dim*1, dim*2, dim*5, 10*dim*1, 10*dim*2, 10*dim*5, 100*dim*1, 100*dim*2,
 * 100*dim*5, ... The default value is "1,2,5".
 * - "precision_x: VALUE" defines the precision used when outputting variables and corresponds to the number
 * of digits to be printed after the decimal point. The default value is 8.
 * - "precision_f: VALUE" defines the precision used when outputting f values and corresponds to the number of
 * digits to be printed after the decimal point. The default value is 15.
 *
 * @return The constructed observer object or NULL if observer_name equals NULL, "" or "no_observer".
 */
coco_observer_t *coco_observer(const char *observer_name, const char *observer_options) {

  coco_observer_t *observer;
  char *path, *result_folder, *algorithm_name, *algorithm_info;
  const char *outer_folder_name = "exdata";
  int precision_x, precision_f;

  size_t number_target_triggers;
  size_t number_evaluation_triggers;
  double target_precision;
  char *base_evaluation_triggers;

  coco_option_keys_t *known_option_keys, *given_option_keys, *additional_option_keys, *redundant_option_keys;

  /* Sets the valid keys for observer options
   * IMPORTANT: This list should be up-to-date with the code and the documentation */
  const char *known_keys[] = { "result_folder", "algorithm_name", "algorithm_info",
      "number_target_triggers", "target_precision", "number_evaluation_triggers", "base_evaluation_triggers",
      "precision_x", "precision_f" };
  additional_option_keys = NULL; /* To be set by the chosen observer */

  if (0 == strcmp(observer_name, "no_observer")) {
    return NULL;
  } else if (strlen(observer_name) == 0) {
    coco_warning("Empty observer_name has no effect. To prevent this warning use 'no_observer' instead");
    return NULL;
  }

  result_folder = coco_allocate_string(COCO_PATH_MAX);
  algorithm_name = coco_allocate_string(COCO_PATH_MAX);
  algorithm_info = coco_allocate_string(5 * COCO_PATH_MAX);
  /* Read result_folder, algorithm_name and algorithm_info from the observer_options and use
   * them to initialize the observer */
  if (coco_options_read_string(observer_options, "result_folder", result_folder) == 0) {
    strcpy(result_folder, "default");
  }
  /* Create the result_folder inside the "exdata" folder */
  path = coco_allocate_string(COCO_PATH_MAX);
  memcpy(path, outer_folder_name, strlen(outer_folder_name) + 1);
  coco_join_path(path, COCO_PATH_MAX, result_folder, NULL);
  coco_create_unique_directory(&path);
  coco_info("Results will be output to folder %s", path);

  if (coco_options_read_string(observer_options, "algorithm_name", algorithm_name) == 0) {
    strcpy(algorithm_name, "ALG");
  }

  if (coco_options_read_string(observer_options, "algorithm_info", algorithm_info) == 0) {
    strcpy(algorithm_info, "");
  }

  number_target_triggers = 100;
  if (coco_options_read_size_t(observer_options, "number_target_triggers", &number_target_triggers) != 0) {
    if (number_target_triggers == 0)
      number_target_triggers = 100;
  }

  target_precision = 1e-8;
  if (coco_options_read_double(observer_options, "target_precision", &target_precision) != 0) {
    if ((target_precision > 1) || (target_precision <= 0))
      target_precision = 1e-8;
  }

  number_evaluation_triggers = 20;
  if (coco_options_read_size_t(observer_options, "number_evaluation_triggers", &number_evaluation_triggers) != 0) {
    if (number_evaluation_triggers < 4)
      number_evaluation_triggers = 20;
  }

  base_evaluation_triggers = coco_allocate_string(COCO_PATH_MAX);
  if (coco_options_read_string(observer_options, "base_evaluation_triggers", base_evaluation_triggers) == 0) {
    strcpy(base_evaluation_triggers, "1,2,5");
  }

  precision_x = 8;
  if (coco_options_read_int(observer_options, "precision_x", &precision_x) != 0) {
    if ((precision_x < 1) || (precision_x > 32))
      precision_x = 8;
  }

  precision_f = 15;
  if (coco_options_read_int(observer_options, "precision_f", &precision_f) != 0) {
    if ((precision_f < 1) || (precision_f > 32))
      precision_f = 15;
  }

  observer = coco_observer_allocate(path, observer_name, algorithm_name, algorithm_info,
      number_target_triggers, target_precision, number_evaluation_triggers, base_evaluation_triggers,
      precision_x, precision_f);

  coco_free_memory(path);
  coco_free_memory(result_folder);
  coco_free_memory(algorithm_name);
  coco_free_memory(algorithm_info);
  coco_free_memory(base_evaluation_triggers);

  /* Here each observer must have an entry - a call to a specific function that sets the additional_option_keys
   * and the following observer fields:
   * - logger_allocate_function
   * - logger_free_function
   * - data_free_function
   * - data */
  if (0 == strcmp(observer_name, "toy")) {
    observer_toy(observer, observer_options, &additional_option_keys);
  } else if (0 == strcmp(observer_name, "bbob")) {
    observer_bbob(observer, observer_options, &additional_option_keys);
  } else if (0 == strcmp(observer_name, "bbob-biobj")) {
    observer_biobj(observer, observer_options, &additional_option_keys);
  } else {
    coco_warning("Unknown observer!");
    return NULL;
  }

  /* Check for redundant option keys */
  known_option_keys = coco_option_keys_allocate(sizeof(known_keys) / sizeof(char *), known_keys);
  coco_option_keys_add(&known_option_keys, additional_option_keys);
  given_option_keys = coco_option_keys(observer_options);

  if (given_option_keys) {
    redundant_option_keys = coco_option_keys_get_redundant(known_option_keys, given_option_keys);

    if ((redundant_option_keys != NULL) && (redundant_option_keys->count > 0)) {
      /* Warn the user that some of given options are being ignored and output the valid options */
      char *output_redundant = coco_option_keys_get_output_string(redundant_option_keys,
          "coco_observer(): Some keys in observer options were ignored:\n");
      char *output_valid = coco_option_keys_get_output_string(known_option_keys,
          "Valid keys for observer options are:\n");
      coco_warning("%s%s", output_redundant, output_valid);
      coco_free_memory(output_redundant);
      coco_free_memory(output_valid);
    }

    coco_option_keys_free(given_option_keys);
    coco_option_keys_free(redundant_option_keys);
  }
  coco_option_keys_free(known_option_keys);
  coco_option_keys_free(additional_option_keys);

  return observer;
}

/**
 * Wraps the observer's logger around the problem if the observer is not NULL and invokes the initialization
 * of this logger.
 *
 * @param problem The given COCO problem.
 * @param observer The COCO observer, whose logger will wrap the problem.
 *
 * @returns The observed problem in the form of a new COCO problem instance or the same problem if the
 * observer is NULL.
 */
coco_problem_t *coco_problem_add_observer(coco_problem_t *problem, coco_observer_t *observer) {

  if (problem == NULL)
	  return NULL;

  if ((observer == NULL) || (observer->is_active == 0)) {
    coco_warning("The problem will not be observed. %s",
        observer == NULL ? "(observer == NULL)" : "(observer not active)");
    return problem;
  }

  assert(observer->logger_allocate_function);
  return observer->logger_allocate_function(observer, problem);
}

/**
 * Frees the observer's logger and returns the inner problem.
 *
 * @param problem The observed COCO problem.
 * @param observer The COCO observer, whose logger was wrapping the problem.
 *
 * @returns The unobserved problem as a pointer to the inner problem or the same problem if the problem
 * was not observed.
 */
coco_problem_t *coco_problem_remove_observer(coco_problem_t *problem, coco_observer_t *observer) {

  coco_problem_t *problem_unobserved;
  char *prefix;

  if ((observer == NULL) || (observer->is_active == 0)) {
    coco_warning("The problem was not observed. %s",
        observer == NULL ? "(observer == NULL)" : "(observer not active)");
    return problem;
  }

  /* Check that we are removing the observer that is actually wrapping the problem.
   *
   * This is a hack - it assumes that the name of the problem is formatted as "observer_name(problem_name)".
   * While not elegant, it does the job and is better than nothing. */
  prefix = coco_remove_from_string(problem->problem_name, "(", "");
  if (strcmp(prefix, observer->observer_name) != 0) {
    coco_error("coco_problem_remove_observer(): trying to remove observer %s instead of %s",
        observer->observer_name, prefix);
  }
  coco_free_memory(prefix);

  /* Keep the inner problem and remove the logger data */
  problem_unobserved = coco_problem_transformed_get_inner_problem(problem);
  coco_problem_transformed_free_data(problem);
  problem = NULL;

  return problem_unobserved;
}
#line 1 "code-experiments/src/coco_archive.c"
/**
 * @file coco_archive.c
 * @brief Definitions of functions regarding COCO archives.
 *
 * COCO archives are used to do some pre-processing on the bi-objective archive files. Namely, through a
 * wrapper written in Python, these functions are used to merge archives and compute their hypervolumes.
 */

#line 10 "code-experiments/src/coco_archive.c"
#line 11 "code-experiments/src/coco_archive.c"
#line 12 "code-experiments/src/coco_archive.c"
#line 13 "code-experiments/src/coco_archive.c"

/**
 * @brief The COCO archive structure.
 *
 * The archive structure is used for pre-processing archives of non-dominated solutions.
 */
struct coco_archive_s {

  avl_tree_t *tree;              /**< @brief The AVL tree with non-dominated solutions. */
  double *ideal;                 /**< @brief The ideal point. */
  double *nadir;                 /**< @brief The nadir point. */

  size_t number_of_objectives;   /**< @brief Number of objectives (clearly equal to 2). */

  int is_up_to_date;             /**< @brief Whether archive fields have been updated since last addition. */
  size_t number_of_solutions;    /**< @brief Number of solutions in the archive. */
  double hypervolume;            /**< @brief Hypervolume of the solutions in the archive. */

  avl_node_t *current_solution;  /**< @brief Current solution (to return). */
  avl_node_t *extreme1;          /**< @brief Pointer to the first extreme solution. */
  avl_node_t *extreme2;          /**< @brief Pointer to the second extreme solution. */
  int extremes_already_returned; /**< @brief Whether the extreme solutions have already been returned. */
};

/**
 * @brief The type for the node's item in the AVL tree used by the archive.
 *
 * Contains information on the rounded normalized objective values (normalized_y), which are used for
 * computing the indicators and the text, which is used for output.
 */
typedef struct {
  double *normalized_y;      /**< @brief The values of normalized objectives of this solution. */
  char *text;                /**< @brief The text describing the solution (the whole line of the archive). */
} coco_archive_avl_item_t;

/**
 * @brief Creates and returns the information on the solution in the form of a node's item in the AVL tree.
 */
static coco_archive_avl_item_t* coco_archive_node_item_create(const double *y,
                                                              const double *ideal,
                                                              const double *nadir,
                                                              const size_t num_obj,
                                                              const char *text) {

  /* Allocate memory to hold the data structure coco_archive_avl_item_t */
  coco_archive_avl_item_t *item = (coco_archive_avl_item_t*) coco_allocate_memory(sizeof(*item));

  /* Compute the normalized y */
  item->normalized_y = mo_normalize(y, ideal, nadir, num_obj);

  item->text = coco_strdup(text);
  return item;
}

/**
 * @brief Frees the data of the given coco_archive_avl_item_t.
 */
static void coco_archive_node_item_free(coco_archive_avl_item_t *item, void *userdata) {
  coco_free_memory(item->normalized_y);
  coco_free_memory(item->text);
  coco_free_memory(item);
  (void) userdata; /* To silence the compiler */
}

/**
 * @brief Defines the ordering of AVL tree nodes based on the value of the last objective.
 */
static int coco_archive_compare_by_last_objective(const coco_archive_avl_item_t *item1,
                                                  const coco_archive_avl_item_t *item2,
                                                  void *userdata) {
  if (coco_double_almost_equal(item1->normalized_y[1], item2->normalized_y[1], mo_precision))
    return 0;
  else if (item1->normalized_y[1] < item2->normalized_y[1])
    return -1;
  else
    return 1;

  (void) userdata; /* To silence the compiler */
}

/**
 * @brief Allocates memory for the archive and initializes its fields.
 */
static coco_archive_t *coco_archive_allocate(void) {

  /* Allocate memory to hold the data structure coco_archive_t */
  coco_archive_t *archive = (coco_archive_t*) coco_allocate_memory(sizeof(*archive));

  /* Initialize the AVL tree */
  archive->tree = avl_tree_construct((avl_compare_t) coco_archive_compare_by_last_objective,
      (avl_free_t) coco_archive_node_item_free);

  archive->ideal = NULL;                /* To be allocated in coco_archive() */
  archive->nadir = NULL;                /* To be allocated in coco_archive() */
  archive->number_of_objectives = 2;
  archive->is_up_to_date = 0;
  archive->number_of_solutions = 0;
  archive->hypervolume = 0.0;

  archive->current_solution = NULL;
  archive->extreme1 = NULL;             /* To be set in coco_archive() */
  archive->extreme2 = NULL;             /* To be set in coco_archive() */
  archive->extremes_already_returned = 0;

  return archive;
}

/**
 * The archive always contains the two extreme solutions
 */
coco_archive_t *coco_archive(const char *suite_name,
                             const size_t function,
                             const size_t dimension,
                             const size_t instance) {

  coco_archive_t *archive = coco_archive_allocate();
  int output_precision = 15;
  coco_suite_t *suite;
  char *suite_instance = coco_strdupf("instances: %lu", (unsigned long) instance);
  char *suite_options = coco_strdupf("dimensions: %lu function_indices: %lu",
  		(unsigned long) dimension, (unsigned long) function);
  coco_problem_t *problem;
  char *text;
  int update;

  suite = coco_suite(suite_name, suite_instance, suite_options);
  if (suite == NULL) {
    coco_error("coco_archive(): cannot create suite '%s'", suite_name);
    return NULL; /* Never reached */
  }
  problem = coco_suite_get_next_problem(suite, NULL);
  if (problem == NULL) {
    coco_error("coco_archive(): cannot create problem f%02lu_i%02lu_d%02lu in suite '%s'",
    		(unsigned long) function, (unsigned long) instance, (unsigned long) dimension, suite_name);
    return NULL; /* Never reached */
  }

  /* Store the ideal and nadir points */
  archive->ideal = coco_duplicate_vector(problem->best_value, 2);
  archive->nadir = coco_duplicate_vector(problem->nadir_value, 2);

  /* Add the extreme points to the archive */
  text = coco_strdupf("0\t%.*e\t%.*e\n", output_precision, archive->nadir[0], output_precision, archive->ideal[1]);
  update = coco_archive_add_solution(archive, archive->nadir[0], archive->ideal[1], text);
  coco_free_memory(text);
  assert(update == 1);

  text = coco_strdupf("0\t%.*e\t%.*e\n", output_precision, archive->ideal[0], output_precision, archive->nadir[1]);
  update = coco_archive_add_solution(archive, archive->ideal[0], archive->nadir[1], text);
  coco_free_memory(text);
  assert(update == 1);

  archive->extreme1 = archive->tree->head;
  archive->extreme2 = archive->tree->tail;
  assert(archive->extreme1 != archive->extreme2);

  coco_free_memory(suite_instance);
  coco_free_memory(suite_options);
  coco_suite_free(suite);

  return archive;
}

int coco_archive_add_solution(coco_archive_t *archive, const double y1, const double y2, const char *text) {

  coco_archive_avl_item_t* insert_item;
  avl_node_t *node, *next_node;
  int update = 0;
  int dominance;

  double *y = coco_allocate_vector(2);
  y[0] = y1;
  y[1] = y2;
  insert_item = coco_archive_node_item_create(y, archive->ideal, archive->nadir,
      archive->number_of_objectives, text);
  coco_free_memory(y);

  /* Find the first point that is not worse than the new point (NULL if such point does not exist) */
  node = avl_item_search_right(archive->tree, insert_item, NULL);

  if (node == NULL) {
    /* The new point is an extreme point */
    update = 1;
    next_node = archive->tree->head;
  } else {
    dominance = mo_get_dominance(insert_item->normalized_y, ((coco_archive_avl_item_t*) node->item)->normalized_y,
        archive->number_of_objectives);
    if (dominance > -1) {
      update = 1;
      next_node = node->next;
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
      	assert((node != archive->extreme1) && (node != archive->extreme2));
      	avl_node_delete(archive->tree, node);
      }
    } else {
      /* The new point is dominated or equal to an existing one, ignore */
      update = 0;
    }
  }

  if (!update) {
    coco_archive_node_item_free(insert_item, NULL);
  } else {
    /* Perform tree update */
    while (next_node != NULL) {
      /* Check the dominance relation between the new node and the next node. There are only two possibilities:
       * dominance = 0: the new node and the next node are nondominated
       * dominance = 1: the new node dominates the next node */
      node = next_node;
      dominance = mo_get_dominance(insert_item->normalized_y, ((coco_archive_avl_item_t*) node->item)->normalized_y,
          archive->number_of_objectives);
      if (dominance == 1) {
        next_node = node->next;
        /* The new point dominates the next point, remove the next point */
        assert((node != archive->extreme1) && (node != archive->extreme2));
      	avl_node_delete(archive->tree, node);
      } else {
        break;
      }
    }

    if(avl_item_insert(archive->tree, insert_item) == NULL) {
      coco_warning("Solution %s did not update the archive", text);
      update = 0;
    }

    archive->is_up_to_date = 0;
  }

  return update;
}

/**
 * @brief Updates the archive fields returned by the getters.
 */
static void coco_archive_update(coco_archive_t *archive) {

  double hyp;

  if (!archive->is_up_to_date) {

    avl_node_t *node, *left_node;
    coco_archive_avl_item_t *node_item, *left_node_item;

    /* Updates number_of_solutions */

    archive->number_of_solutions = avl_count(archive->tree);

    /* Updates hypervolume */

    node = archive->tree->head;
    archive->hypervolume = 0; /* Hypervolume of the extreme point equals 0 */
    while (node->next) {
      /* Add hypervolume contributions of the other points that are within ROI */
      left_node = node->next;
      node_item = (coco_archive_avl_item_t *) node->item;
      left_node_item = (coco_archive_avl_item_t *) left_node->item;
      if (mo_is_within_ROI(left_node_item->normalized_y, archive->number_of_objectives)) {
        hyp = 0;
        if (mo_is_within_ROI(node_item->normalized_y, archive->number_of_objectives))
          hyp = (node_item->normalized_y[0] - left_node_item->normalized_y[0]) * (1 - left_node_item->normalized_y[1]);
        else
          hyp = (1 - left_node_item->normalized_y[0]) * (1 - left_node_item->normalized_y[1]);
        assert(hyp >= 0);
         archive->hypervolume += hyp;
      }
      node = left_node;
    }

    archive->is_up_to_date = 1;
    archive->current_solution = NULL;
    archive->extremes_already_returned = 0;
  }

}

const char *coco_archive_get_next_solution_text(coco_archive_t *archive) {

  char *text;

  coco_archive_update(archive);

  if (!archive->extremes_already_returned) {

    if (archive->current_solution == NULL) {
      /* Return the first extreme */
      text = ((coco_archive_avl_item_t *) archive->extreme1->item)->text;
      archive->current_solution = archive->extreme2;
      return text;
    }

    if (archive->current_solution == archive->extreme2) {
      /* Return the second extreme */
      text = ((coco_archive_avl_item_t *) archive->extreme2->item)->text;
      archive->extremes_already_returned = 1;
      archive->current_solution = archive->tree->head;
      return text;
    }

  } else {

    if (archive->current_solution == NULL)
      return "";

    if ((archive->current_solution == archive->extreme1) || (archive->current_solution == archive->extreme2)) {
      /* Skip this one */
      archive->current_solution = archive->current_solution->next;
      return coco_archive_get_next_solution_text(archive);
    }

    /* Return the current solution and move to the next */
    text = ((coco_archive_avl_item_t *) archive->current_solution->item)->text;
    archive->current_solution = archive->current_solution->next;
    return text;
  }

  return NULL; /* This point should never be reached. */
}

size_t coco_archive_get_number_of_solutions(coco_archive_t *archive) {
  coco_archive_update(archive);
  return archive->number_of_solutions;
}

double coco_archive_get_hypervolume(coco_archive_t *archive) {
  coco_archive_update(archive);
  return archive->hypervolume;
}

void coco_archive_free(coco_archive_t *archive) {

  assert(archive != NULL);

  avl_tree_destruct(archive->tree);
  coco_free_memory(archive->ideal);
  coco_free_memory(archive->nadir);
  coco_free_memory(archive);

}
#line 1 "code-experiments/src/coco_runtime_c.c"
/**
 * @file coco_runtime_c.c
 * @brief Generic COCO runtime implementation for the C language.
 *
 * Other language interfaces might want to replace this so that memory allocation and error handling go
 * through the respective language runtime.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#line 14 "code-experiments/src/coco_runtime_c.c"
#line 15 "code-experiments/src/coco_runtime_c.c"

void coco_error(const char *message, ...) {
  va_list args;

  fprintf(stderr, "COCO FATAL ERROR: ");
  va_start(args, message);
  vfprintf(stderr, message, args);
  va_end(args);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

void coco_warning(const char *message, ...) {
  va_list args;

  if (coco_log_level >= COCO_WARNING) {
    fprintf(stderr, "COCO WARNING: ");
    va_start(args, message);
    vfprintf(stderr, message, args);
    va_end(args);
    fprintf(stderr, "\n");
  }
}

void coco_info(const char *message, ...) {
  va_list args;

  if (coco_log_level >= COCO_INFO) {
    fprintf(stdout, "COCO INFO: ");
    va_start(args, message);
    vfprintf(stdout, message, args);
    va_end(args);
    fprintf(stdout, "\n");
    fflush(stdout);
  }
}

/**
 * A function similar to coco_info that prints only the given message without any prefix and without
 * adding a new line.
 */
void coco_info_partial(const char *message, ...) {
  va_list args;

  if (coco_log_level >= COCO_INFO) {
    va_start(args, message);
    vfprintf(stdout, message, args);
    va_end(args);
    fflush(stdout);
  }
}

void coco_debug(const char *message, ...) {
  va_list args;

  if (coco_log_level >= COCO_DEBUG) {
    fprintf(stdout, "COCO DEBUG: ");
    va_start(args, message);
    vfprintf(stdout, message, args);
    va_end(args);
    fprintf(stdout, "\n");
    fflush(stdout);
  }
}

void *coco_allocate_memory(const size_t size) {
  void *data;
  if (size == 0) {
    coco_error("coco_allocate_memory() called with 0 size.");
    return NULL; /* never reached */
  }
  data = malloc(size);
  if (data == NULL)
    coco_error("coco_allocate_memory() failed.");
  return data;
}

void coco_free_memory(void *data) {
  free(data);
}
