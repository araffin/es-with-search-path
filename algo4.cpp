// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
#include "algo4.hpp"
#include <assert.h>

using namespace std;

/**
 * A grid search optimizer that can be used for single- as well as multi-objective optimization.
 *
 * @param evaluate The evaluation function used to evaluate the solutions.
 * @param dimension The number of variables (= n in our case)
 * @param number_of_objectives The number of objectives. = 1 in our case
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
                    const size_t max_budget)
{

  // ======================================================
  // INITIALIZATION
  // ======================================================

  assert(number_of_objectives == 1);

  double X[dimension] = {0}; // To be randomely initialized - double (not float) to be used with the function evaluate
  double X_best[dimension];
  double Sigma[dimension]; // idem
  double s_sigma[dimension] = {0}; // search path --> vector !
  bool happy = false;
  int counter = 0;
  int lambda; // lambda is given, how do we have it ? 
  int mu = (int) lambda/4;
  // d and d_i uninitilized ? --> done !
  double d = 1 + sqrt(mu/dimension); 
  double di = 3*dimension; 
  double c_sigma = sqrt(mu/(dimension + mu)); 
  // offspring population
  double** X_k;
  // Mutation vectors
  double** Z;
  // Matrix Initialization
  X_k = new double*[lambda];
  Z = new double*[lambda];
  for (size_t i = 0; i < lambda; i++ )
  {
    X_k[i] = new double[dimension];
    Z[i] = new double[dimension];
  }
  // initialize matrix elements
  for (size_t i = 0; i < lambda; i++ )
  {
     for (size_t j = 0; j < dimension; j++ )
     {
        X_k[i][j] = 0;
     }
  }

  // Random generator (normal distribution)
  random_device rd;
  mt19937 gen(rd()); //Mersenne Twister 19937 generator
  normal_distribution<> N(0,1);

  // a t-on besoin d'un critère d'arrêt ? Sinon on fait juste le nb d'itérations comme dans les exemples
  double stop_criterion = 0.0002; // juste pour tester les 3 dernières lignes de la boucle

  // printArray(X, dimension);
  // double *x = coco_allocate_vector(dimension);
  double *y = coco_allocate_vector(number_of_objectives);

  for (int i = 0; i < dimension ; i ++)
  {
    Sigma[i] = 1;
  }

  while (!happy && counter < max_budget)
  {

    for (size_t k = 0; k < lambda; k++)
    {
      // z_k = N(0,I)
      normalMatrix(Z, N, gen, lambda, dimension);
      double product_result[dimension];
      elementProduct(Sigma, Z[k], product_result, dimension);
      // x_k = x + sigma o z_k
      arraySum(X, product_result, X_k[k], dimension);
    }
    /*
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
    */

    // utilise t-on vraiment un critère d'arrêt ?
    evaluate(X,y);
    happy = (y[0] < stop_criterion);

    counter++;
  }
  // Free memory
  freeMatrix(X_k, lambda);
  freeMatrix(Z, lambda);
  coco_free_memory(y);
}


void printArray(double* array, size_t n)
{
  for (size_t i = 0; i < n; i++)
    cout << array[i] << " ";

  cout << endl;
}

void printMatrix(double** array, size_t n, size_t m)
{
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < m; j++) {
      cout << array[i][j] << " ";
    }
    cout << endl;
  }
}

void elementProduct(double* a, double* b, double* result, size_t n)
{
  for (size_t i = 0; i < n; i++)
    result[i] = a[i]*b[i];
}

void arraySum(double* a, double* b, double* result, size_t n)
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
void sumVectors(double** mat, double* result, size_t k, size_t n)
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
void normalMatrix(double** mat, normal_distribution<> N, mt19937 gen,
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

void freeMatrix(double** mat, size_t nb_rows)
{
  for (size_t i = 0; i < nb_rows; i++ )
  {
     delete mat[i];
  }
  delete mat;
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
