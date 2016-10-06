// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
#include <iostream>
#include <math.h>
#include <random> // Normal distribution

#include "coco.h"

using namespace std;


/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);


/**
 * A grid search optimizer that can be used for single- as well as multi-objective optimization.
 *
 * @param evaluate The evaluation function used to evaluate the solutions.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param max_budget The maximal number of evaluations.
 *
 * If max_budget is not enough to cover even the smallest possible grid, only the first max_budget
 * nodes of the grid are evaluated.
 */

void algo4(evaluate_function_t evaluate,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t max_budget) {

  // ======================================================
  // INITIALIZATION
  // ======================================================
  /*
  float X[n]; // To be randomely initialized
  float X_best[n];
  float Sigma[n]; // idem
  float s_sigma = 0; // search path ? not a vector ?
  bool happy = false;
  int counter = 0;
  */



}

#if 0
void my_grid_search(evaluate_function_t evaluate,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t max_budget){

  double *x = coco_allocate_vector(dimension);
  double *y = coco_allocate_vector(number_of_objectives);
  long *nodes = (long *) coco_allocate_memory(sizeof(long) * dimension);
  double *grid_step = coco_allocate_vector(dimension);
  size_t i, j;
  size_t evaluations = 0;
  long max_nodes = (long) floor(pow((double) max_budget, 1.0 / (double) dimension)) - 1;

  /* Take care of the borderline case */
  if (max_nodes < 1) max_nodes = 1;

  /* Initialization */
  for (j = 0; j < dimension; j++) {
    nodes[j] = 0;
    grid_step[j] = (upper_bounds[j] - lower_bounds[j]) / (double) max_nodes;
  }

  while (evaluations < max_budget) {

    /* Construct x and evaluate it */
    for (j = 0; j < dimension; j++) {
      x[j] = lower_bounds[j] + grid_step[j] * (double) nodes[j];
    }

    /* Call the evaluate function to evaluate x on the current problem (this is where all the COCO logging
     * is performed) */
    evaluate(x, y);
    evaluations++;

    /* Inside the grid, move to the next node */
    if (nodes[0] < max_nodes) {
      nodes[0]++;
    }

    /* At an outside node of the grid, move to the next level */
    else if (max_nodes > 0) {
      for (j = 1; j < dimension; j++) {
        if (nodes[j] < max_nodes) {
          nodes[j]++;
          for (i = 0; i < j; i++)
            nodes[i] = 0;
          break;
        }
      }

      /* At the end of the grid, exit */
      if ((j == dimension) && (nodes[j - 1] == max_nodes))
        break;
    }
  }

  coco_free_memory(x);
  coco_free_memory(y);
  coco_free_memory(nodes);
  coco_free_memory(grid_step);
}

#endif




#if 0
void printArray(float* array, size_t n)
{
  for (size_t i = 0; i < n; i++)
    cout << array[i] << " ";
  
  cout << endl;
}

void printMatrix(float** array, size_t n, size_t m)
{
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < m; j++) {
      cout << array[i][j] << " ";
    }
    cout << endl;
  }
}

void elementProduct(float* a, float* b, float* result, size_t n)
{
  for (size_t i = 0; i < n; i++)
    result[i] = a[i]*b[i];
}

void arraySum(float* a, float* b, float* result, size_t n)
{
  for (size_t i = 0; i < n; i++)
    result[i] = a[i] + b[i];
}

/**
 * Sum component wise an array of k vectors of dimension n
 * In fact, it sums the rows of a matrix of size k*n
 * @param mat    a matrix
 * @param result an array containing the result of the operation
 * @param k  number of vectors (= nb of rows)
 * @param n the dimension (= nb of columns)
 */
void sumVectors(float** mat, float* result, size_t k, size_t n)
{
  for (size_t j = 0; j < n; j++) {
    result[j] = 0;
    for (size_t i = 0; i < k; i++) {
      result[j] += mat[i][j];
    }
  }
}

/**
 * Fill a matrix with random numbers drawn from a normal distribution
 * @param mat     a matrix of dimension nb_rows x nb_cols
 * @param N       the normal distribution
 * @param gen     a random generator
 * @param nb_rows 
 * @param nb_cols 
 */
void normalMatrix(float** mat, normal_distribution<> N, mt19937 gen, 
                  size_t nb_rows, size_t nb_cols)
{
  for (size_t i = 0; i < nb_rows; i++)
  {
    for (size_t j = 0; j < nb_cols; j++)
    {
      mat[i][j] = N(gen);
    }
  }
}

void fmin(void* fitnessFunction, int n, int lambda, int maxIterations)
{
  int mu = (int) lambda/4;
  float cs = sqrt(mu/(n+mu)); // normalization factor?
  float d = 1 + sqrt(mu/n); // damping factor for cs ?
  float d_i = 3*n; // or vector ?
  
  // Random generator (normal distribution)
  random_device rd;
  mt19937 gen(rd()); //Mersenne Twister 19937 generator
  normal_distribution<> N(0,1);
  
  // ======================================================
  // INITIALIZATION
  // ======================================================
  float X[n]; // To be randomely initialized
  float X_best[n];
  float Sigma[n]; // idem
  float s_sigma = 0; // search path ? not a vector ?
  bool happy = false;
  int counter = 0;
  bool population[mu];
  while (!happy && counter < maxIterations)
  {
    float Z[lambda][n];
    float X_k[lambda][n];
    
    for (size_t k = 0; k < lambda; k++)
    {
      normalVector(Z, N, gen, lambda, n);
      X_k[k] = arraySum(X + elementProduct(Sigma, Z[k]), n);
    }
    // Select the mu best solution + update X_best
    select_mu_best(mu, lambda, X_k, Z, fitnessFunction, population);
    
    // Better to do that inline (because we need several output)
    s_sigma = (1-cs)*s_sigma + sqrt(cs*(2-cs))*(sqrt(mu)/mu) * sumVectors(Z, population, lambda, n);
    Sigma = elementProduct(Sigma, \
      exp(abs(s_sigma)/(2*d_i) - 1) \
      * exp((cs/d)*abs(s_sigma)/(esperanceNormalDistri) - 1)
      ,n);
    // Warning here, sum of vectors
    X = (1/mu)*sumVectors(X_k, population, lambda, n);
    
    counter++;
  }
  std::cout << "/* Best Solution :  */" << printArray(X_best, n) << std::endl;  
}
#endif