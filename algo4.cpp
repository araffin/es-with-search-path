// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
// Post Process : python -m bbob_pproc [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
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
                    const size_t max_budget, 
                    coco_random_state_t *random_generator)
{

  // ======================================================
  // INITIALIZATION
  // ======================================================

  assert(number_of_objectives == 1);

  double X[dimension] = {0}; // To be randomely initialized - double (not float) to be used with the function evaluate
  double X_tmp[dimension] = {0}; // Used when sorting the solutions
  double Z_tmp[dimension] = {0}; // Used when sorting the solutions
  double Sigma[dimension]; // idem
  double s_sigma[dimension] = {0}; // search path --> vector !
  bool happy = false;
  int counter = 0;
  size_t lambda = 50;
  size_t mu = (size_t) lambda/4;
  // Damping factors
  double d = 1 + sqrt((double)mu/(double)dimension); 
  size_t di = 3*dimension; 
  double c_sigma = sqrt((double)mu/(double)(dimension + mu)); 
  // offspring population
  double** X_k;
  // Mutation vectors
  double** Z;
  double fitness[lambda];
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
        X_k[i][j] = lower_bounds[j] + (upper_bounds[j] - lower_bounds[j])*coco_random_uniform(random_generator);
     }
  }

  // Random generator (normal distribution)
  random_device rd;
  mt19937 gen(rd()); //Mersenne Twister 19937 generator
  normal_distribution<> N(0,1);

  // a t-on besoin d'un critère d'arrêt ? Sinon on fait juste le nb d'itérations comme dans les exemples
  double stop_criterion = 0.0002; // juste pour tester les 3 dernières lignes de la boucle

  // double *x = coco_allocate_vector(dimension);
  double *y = coco_allocate_vector(number_of_objectives);

  for (size_t j = 0; j < dimension ; j++)
  {
    Sigma[j] = 1;
  }

  while (!happy && counter < max_budget)
  {

    for (size_t k = 0; k < lambda; k++)
    {
      // z_k = N(0,I)
      normalMatrix(Z, N, gen, lambda, dimension);

      // temp variable to store intermediate result
      double product_result[dimension];
      elementProduct(Sigma, Z[k], product_result, dimension);
      // x_k = x + sigma o z_k
      arraySum(X, product_result, X_k[k], dimension);

      //check boundaries
      for(size_t j = 0; j < dimension; j++)
      {
        X_k[k][j] = fmax(lower_bounds[j], fmin(X_k[k][j], upper_bounds[j]));
      }
    }


    // Select the mu best
    for (size_t k = 0; k < lambda; k++)
    {
      // X_tmp = X_k[k]
      vectorCopy(X_k[k], X_tmp, dimension);
      // Z_tmp = Z[k]
      vectorCopy(Z[k], Z_tmp, dimension);

      //Evaluate f(x_k) and store the result in y
      evaluate(X_tmp, y);

      //Variable used to sort the X_k array
      size_t position = k;
      if(k > 0)
      {
        //Find the position of x_k given f(x_k)
        while(*y <= fitness[position-1] && position > 0) position--;
        for (size_t l = k; l > position; l--)
        {
          //Move x_k and f(x_k) backwards till position
          fitness[l] = fitness[l - 1];
          vectorCopy(X_k[l - 1], X_k[l], dimension);
          vectorCopy(Z[l - 1], Z[l], dimension);
        }
      }
      // Repostion f(x_k)
      fitness[position] = *y;
      // Reposition x_k : X_k[pos] = X_tmp
      vectorCopy(X_tmp, X_k[position], dimension);
      // Reposition z_k : Z[pos] = Z_tmp
      vectorCopy(Z_tmp, Z[position], dimension);

    }

    // uptdate s_sigma // TO BE CHECKED => CHECKED by toni
    for(size_t j = 0; j < dimension; j++)
    {
      double sum = 0;  // sum = sum(zk); zk in P
      for(size_t k = 0; k < mu; k++) // we take the mu best in Z
      {
        sum = sum + Z[k][j];
      }
      s_sigma[j] = (1-c_sigma)*s_sigma[j] \
                  + sqrt(c_sigma*(2 - c_sigma)*(double)mu)/double(mu) * sum; 
    }

    // update Sigma
    double exp1=0;
    double exp2 = exp((normE(s_sigma, dimension) / (sqrt(dimension)*(1 - 1/4*dimension + 1/(21*pow(dimension, 2))))) -1);
    exp2 = pow(exp2, c_sigma/d);
    for(size_t j = 0; j < dimension; j++)
    {
      exp1 = exp((abs(Sigma[j]) / (sqrt(2)/sqrt(coco_pi))) -1);
      exp1 = pow(exp1, 1/di);
      Sigma[j] = Sigma[j]*exp1*exp2;
    }


    // update X
    for(size_t j = 0; j < dimension; j++)
    {
      double sumXk = 0;
      for (size_t k = 0; k < mu; k++)
      {
        sumXk = sumXk + X_k[k][j];
      }
      X[j] = (double)sumXk/(double)mu; 

    }


    /*
    // Better to do that inline (because we need several output)
    s_sigma = (1-cs)*s_sigma + sqrt(cs*(2-cs))*(sqrt(mu)/mu) * sumVectors(Z, population, lambda, n);
    Sigma = elementProduct(Sigma, \
      exp(abs(s_sigma)/(2*d_i) - 1) \
      * exp((cs/d)*abs(s_sigma)/(esperanceNormalDistri) - 1)
      ,n);
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


/* Selection function
So far I have considered that this function will classify the best member of X_k and Z at the fist positions
ie P will be the mu first columns of X_k and Z (to avoid having new matrix) but I don't know if it is the best way

*/

void select_mu_best(double mu, double lambda, double** X_k, double** Z, evaluate_function_t evaluate)
{

}

double normE(double *m, size_t dimension)
{
  double out=0;
  for (int i=0; i<dimension; ++i)
  {
    out += pow(m[i], 2);
  }
  return sqrt(out);
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

void vectorCopy(double* a, double* b, size_t n)
{
  for (size_t i = 0; i < n; i++)
    b[i] = a[i];
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
