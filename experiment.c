/**
 * Benchmarking ES-With-Search-Path algorithm on a COCO suite.
 *
 * Set the global parameter BUDGET_MULTIPLIER to suit your needs.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "algo4.hpp"

/**
 * The maximal budget for evaluations done by an optimization algorithm equals dimension * BUDGET_MULTIPLIER.
 * Increase the budget multiplier value gradually to see how it affects the runtime.
 */
static const size_t BUDGET_MULTIPLIER = 1e4;

/**
 * The maximal number of independent restarts allowed for an algorithm that restarts itself.
 */
static const size_t INDEPENDENT_RESTARTS = 1e5;

/**
 * The random seed. Change if needed.
 */
static const uint32_t RANDOM_SEED = 0xdeadbeef;

/**
 * A pointer to the problem to be optimized (needed in order to simplify the interface between the optimization
 * algorithm and the COCO platform).
 */
static coco_problem_t *PROBLEM;

/**
 * The function that calls the evaluation of the first vector on the problem to be optimized and stores the
 * evaluation result in the second vector.
 */
static void evaluate_function(const double *x, double *y) {
  coco_evaluate_function(PROBLEM, x, y);
}

/* Declarations of all functions implemented in this file (so that their order is not important): */
void experiment(const char *suite_name,
                        const char *observer_name,
                        coco_random_state_t *random_generator);

void my_random_search(evaluate_function_t evaluate,
                      const size_t dimension,
                      const size_t number_of_objectives,
                      const double *lower_bounds,
                      const double *upper_bounds,
                      const size_t max_budget,
                      coco_random_state_t *random_generator);

/* Structure and functions needed for timing the experiment */
typedef struct {
	size_t number_of_dimensions;
	size_t current_idx;
	char **output;
	size_t previous_dimension;
	size_t cumulative_evaluations;
	time_t start_time;
	time_t overall_start_time;
} timing_data_t;
static timing_data_t *timing_data_initialize(coco_suite_t *suite);
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem);
static void timing_data_finalize(timing_data_t *timing_data);

/**
 * The main method initializes the random number generator and calls the example experiment on the
 * bi-objective suite.
 */
int main(void) {

  coco_random_state_t *random_generator = coco_random_new(RANDOM_SEED);

  /* Change the log level to "warning" to get less output */
  coco_set_log_level("info");

  printf("Running the experiment... (might take time, be patient)\n");
  fflush(stdout);

  // Run the experiment on the bbob suite
  experiment("bbob", "bbob", random_generator);

  printf("Done!\n");
  fflush(stdout);

  coco_random_free(random_generator);

  return 0;
}

/**
 * A simple example of benchmarking random search on a suite with instances from 2016 that can serve also as
 * a timing experiment.
 *
 * @param suite_name Name of the suite (use "bbob" for the single-objective and "bbob-biobj" for the
 * bi-objective suite).
 * @param observer_name Name of the observer (use "bbob" for the single-objective and "bbob-biobj" for the
 * bi-objective observer).
 * @param random_generator The random number generator.
 */
void experiment(const char *suite_name,
                        const char *observer_name,
                        coco_random_state_t *random_generator) {

  size_t run;
  coco_suite_t *suite;
  coco_observer_t *observer;
  timing_data_t *timing_data;

  /* Set some options for the observer. See documentation for other options. */
  char *observer_options =
      coco_strdupf("result_folder: ES-With-Search-Path_on_%s_%d "
                   "algorithm_name: ES-With-Search-Path "
                   "algorithm_info: \"ES-With-Search-Path algorithm\"", suite_name, BUDGET_MULTIPLIER);

  /* Initialize the suite and observer */
  suite = coco_suite(suite_name, "year: 2016", "dimensions: 2,3,5,10,20,40");
  observer = coco_observer(observer_name, observer_options);
  coco_free_memory(observer_options);

  /* Initialize timing */
  timing_data = timing_data_initialize(suite);

  /* Iterate over all problems in the suite */
  while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {

    size_t dimension = coco_problem_get_dimension(PROBLEM);

    /* Run the algorithm at least once */
    for (run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

      size_t evaluations_done = coco_problem_get_evaluations(PROBLEM);
      long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - (long) evaluations_done;

      /* Break the loop if the target was hit or there are no more remaining evaluations */
      if (coco_problem_final_target_hit(PROBLEM) || (evaluations_remaining <= 0))
      break;

      /* Call the optimization algorithm for the remaining number of evaluations */
      algo4(evaluate_function,
        dimension,
        coco_problem_get_number_of_objectives(PROBLEM),
        coco_problem_get_smallest_values_of_interest(PROBLEM),
        coco_problem_get_largest_values_of_interest(PROBLEM),
        (size_t) evaluations_remaining);

        /* Break the loop if the algorithm performed no evaluations or an unexpected thing happened */
        if (coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
          printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
          (unsigned long) evaluations_done, (unsigned long) dimension * BUDGET_MULTIPLIER);
          break;
        }
        else if (coco_problem_get_evaluations(PROBLEM) < evaluations_done)
        coco_error("Something unexpected happened - function evaluations were decreased!");
      }
      /* Keep track of time */
      timing_data_time_problem(timing_data, PROBLEM);
  }

  /* Output and finalize the timing data */
  timing_data_finalize(timing_data);

  coco_observer_free(observer);
  coco_suite_free(suite);

}

/**
 * A random search algorithm that can be used for single- as well as multi-objective optimization.
 *
 * @param evaluate The evaluation function used to evaluate the solutions.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param max_budget The maximal number of evaluations.
 * @param random_generator Pointer to a random number generator able to produce uniformly and normally
 * distributed random numbers.
 */
void my_random_search(evaluate_function_t evaluate,
                      const size_t dimension,
                      const size_t number_of_objectives,
                      const double *lower_bounds,
                      const double *upper_bounds,
                      const size_t max_budget,
                      coco_random_state_t *random_generator) {

  double *x = coco_allocate_vector(dimension);
  double *y = coco_allocate_vector(number_of_objectives);
  double range;
  size_t i, j;

  for (i = 0; i < max_budget; ++i) {

    /* Construct x as a random point between the lower and upper bounds */
    for (j = 0; j < dimension; ++j) {
      range = upper_bounds[j] - lower_bounds[j];
      x[j] = lower_bounds[j] + coco_random_uniform(random_generator) * range;
    }

    /* Call the evaluate function to evaluate x on the current problem (this is where all the COCO logging
     * is performed) */
    evaluate(x, y);

  }

  coco_free_memory(x);
  coco_free_memory(y);
}


/**
 * Allocates memory for the timing_data_t object and initializes it.
 */
static timing_data_t *timing_data_initialize(coco_suite_t *suite) {

	timing_data_t *timing_data = (timing_data_t *) coco_allocate_memory(sizeof(*timing_data));
	size_t function_idx, dimension_idx, instance_idx, i;

	/* Find out the number of all dimensions */
	coco_suite_decode_problem_index(suite, coco_suite_get_number_of_problems(suite) - 1, &function_idx,
			&dimension_idx, &instance_idx);
	timing_data->number_of_dimensions = dimension_idx + 1;
	timing_data->current_idx = 0;
	timing_data->output = (char **) coco_allocate_memory(timing_data->number_of_dimensions * sizeof(char *));
	for (i = 0; i < timing_data->number_of_dimensions; i++) {
		timing_data->output[i] = NULL;
	}
	timing_data->previous_dimension = 0;
	timing_data->cumulative_evaluations = 0;
	time(&timing_data->start_time);
	time(&timing_data->overall_start_time);

	return timing_data;
}

/**
 * Keeps track of the total number of evaluations and elapsed time. Produces an output string when the
 * current problem is of a different dimension than the previous one or when NULL.
 */
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem) {

	double elapsed_seconds = 0;

	if ((problem == NULL) || (timing_data->previous_dimension != coco_problem_get_dimension(problem))) {

		/* Output existing timing information */
		if (timing_data->cumulative_evaluations > 0) {
			time_t now;
			time(&now);
			elapsed_seconds = difftime(now, timing_data->start_time) / (double) timing_data->cumulative_evaluations;
			timing_data->output[timing_data->current_idx++] = coco_strdupf("d=%lu done in %.2e seconds/evaluation\n",
					timing_data->previous_dimension, elapsed_seconds);
		}

		if (problem != NULL) {
			/* Re-initialize the timing_data */
			timing_data->previous_dimension = coco_problem_get_dimension(problem);
			timing_data->cumulative_evaluations = coco_problem_get_evaluations(problem);
			time(&timing_data->start_time);
		}

	} else {
		timing_data->cumulative_evaluations += coco_problem_get_evaluations(problem);
	}
}

/**
 * Outputs and finalizes the given timing data.
 */
static void timing_data_finalize(timing_data_t *timing_data) {

	/* Record the last problem */
	timing_data_time_problem(timing_data, NULL);

  if (timing_data) {
  	size_t i;
  	double elapsed_seconds;
		time_t now;
		int hours, minutes, seconds;

		time(&now);
		elapsed_seconds = difftime(now, timing_data->overall_start_time);

  	printf("\n");
  	for (i = 0; i < timing_data->number_of_dimensions; i++) {
    	if (timing_data->output[i]) {
				printf("%s", timing_data->output[i]);
				coco_free_memory(timing_data->output[i]);
    	}
    }
  	hours = (int) elapsed_seconds / 3600;
  	minutes = ((int) elapsed_seconds % 3600) / 60;
  	seconds = (int)elapsed_seconds - (hours * 3600) - (minutes * 60);
  	printf("Total elapsed time: %dh%02dm%02ds\n", hours, minutes, seconds);

    coco_free_memory(timing_data->output);
    coco_free_memory(timing_data);
  }
}
