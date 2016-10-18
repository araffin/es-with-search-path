// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
// Post Process : python2 -m bbob_pproc [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
#include "algo4.hpp"
#include <assert.h>
#define PI 3.14159265359
#define EVAL_CRIT 3 // The number of step we look at to decide if we are stuck
using namespace std;

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
                    const size_t max_budget)
{
  // ======================================================
  // INITIALIZATION
  // ======================================================

  // Check that we are in mono-objective
  assert(number_of_objectives == 1);

  double X[dimension]; // To be randomely initialized - double (not float) to be used with the function evaluate
  double X_tmp[dimension]; // Used when sorting the solutions
  double Z_tmp[dimension]; // Used when sorting the solutions
  double Sigma[dimension];
  double s_sigma[dimension];
  bool happy = false;

  // Stop if there is no change in the last EVAL_CRIT values greater that stop_criterion
  double stop_criterion = 1e-9;
  // Counter for the number of function evaluation budget
  int counter = 1;
  // adaptative lambda, as suggested in https://arxiv.org/pdf/1604.00772.pdf
  size_t lambda = 8 + int(3*log10(dimension));
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
  double E_HALF_NORMAL = sqrt(2./PI);
  double _dim = double(dimension);
  double E_MULTIDIM_NORMAL = sqrt(_dim)*(1-1./(4*_dim)+1./(21*_dim*_dim));

  uint criterion_counter = 0;
  double last_evaluations[EVAL_CRIT];

  // Matrix Initialization
  X_k = new double*[lambda];
  Z = new double*[lambda];
  for (size_t i = 0; i < lambda; i++ )
  {
    X_k[i] = new double[dimension];
    Z[i] = new double[dimension];
  }

  // Random generator (normal distribution)
  random_device rd;
  mt19937 gen(rd()); //Mersenne Twister 19937 generator
  normal_distribution<> N(0,1);

  // Init Solution and Sigma (step size)
  for (size_t j = 0; j < dimension; j++)
  {
    X[j] = lower_bounds[j] + (upper_bounds[j] - lower_bounds[j]) / 2.0;
    // Based on the +/-3sigma rule to obtain a 99.7% confidence interval (cf cma-es.c code)
    Sigma[j] = (upper_bounds[j] - lower_bounds[j]) / 6.0;
  }

  // With this algorithm, single value
  double *y = coco_allocate_vector(number_of_objectives);

  // ======================================================
  // RUN
  // ======================================================

  // happy depends on the stop_criterion ; counter counts the calls to evaluate
  while (!happy && counter < max_budget)
  {
    for (size_t k = 0; k < lambda; k++)
    {
      for (size_t j = 0; j < dimension; j++) {
        // Create z following a distribution N(0,I)
        Z[k][j] = N(gen);
        X_k[k][j] = X[j] + Sigma[j]*Z[k][j];
        //check boundaries
        X_k[k][j] = fmax(lower_bounds[j], fmin(X_k[k][j], upper_bounds[j]));
      }
    }

    // Evaluate the solutions
    for (size_t k = 0; k < lambda; k++)
    {
      evaluate(X_k[k], y);
      fitness[k] = y[0];
    }

    // Select the mu best (shell sort)
    // Shell sort : better when using greater lambda
    int gaps[8] = {701, 301, 132, 57, 23, 10, 4, 1};

    // // Start with the largest gap and work down to a gap of 1
    for (size_t g = 0; g < 8; g++)
    {
      size_t gap = gaps[g];
      // Do a gapped insertion sort for this gap size.
      // The first gap elements fitness[0..gap-1] are already in gapped order
      // keep adding one more element until the entire array is gap sorted
      for (size_t k = gap; k < lambda; k += 1)
      {
        double y_tmp = fitness[k];
        for (size_t j = 0; j < dimension; j++)
        {
          X_tmp[j] = X_k[k][j];
          Z_tmp[j] = Z[k][j];
        }
        size_t pos = k;
        // shift earlier gap-sorted elements up
        // until the correct location for fitness[k] is found
        while(pos >= gap && fitness[pos-gap] > y_tmp)
        {
          fitness[pos] = fitness[pos-gap];
          for (size_t j = 0; j < dimension; j++)
          {
            X_k[pos][j] = X_k[pos-gap][j];
            Z[pos][j] = Z[pos-gap][j];
          }
          pos -= gap;
        }
        // Reposition f(x_k)
        fitness[pos] = y_tmp;
        for (size_t j = 0; j < dimension; j++)
        {
          // Reposition x_k : X_k[pos] = X_tmp
          X_k[pos][j] = X_tmp[j];
          // Reposition z_k : Z[pos] = Z_tmp
          Z[pos][j] = Z_tmp[j];
        }
      }
    }

    // update s_sigma
    for(size_t j = 0; j < dimension; j++)
    {
      double sum = 0;  // sum = sum(zk); zk in P
      for(size_t k = 0; k < mu; k++) // we take the mu best in Z
      {
        sum = sum + Z[k][j];
      }
      s_sigma[j] = (1-c_sigma)*s_sigma[j] \
                  + sqrt(c_sigma*(2 - c_sigma))*sqrt((double)mu)/double(mu) * sum;
    }

    // update Sigma
    double exp1 = 0;
    double exp2 = exp((c_sigma/d)*((normE(s_sigma, dimension)/E_MULTIDIM_NORMAL)-1));
    for(size_t j = 0; j < dimension; j++)
    {
      exp1 = exp((1./double(di))*((abs(s_sigma[j])/E_HALF_NORMAL)-1));
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
      X[j] = fmax(lower_bounds[j], fmin(X[j], upper_bounds[j]));
    }

    evaluate(X,y);
    // Check if we are still optimising the solution ; stop if we are stuck
    last_evaluations[criterion_counter] = y[0];
    criterion_counter = (criterion_counter + 1) % EVAL_CRIT;
    if (counter >= EVAL_CRIT)
    {
      double min_eval = last_evaluations[0];
      double max_eval = last_evaluations[0];
      // Calcul of the minimum and maximum during the last EVAL_CRIT evaluations
      for (size_t i = 1; i < EVAL_CRIT; i++)
      {
        min_eval = fmin(min_eval, last_evaluations[i]);
        max_eval = fmax(max_eval, last_evaluations[i]);
      }
      happy = ((max_eval - min_eval) < stop_criterion);
    }
    // Update of the counter of evaluation calls. The number of calls per loop is fixed,
    // hence the use of a constant for the update.
    counter += int(lambda) + 1;
  }

  // Free memory
  freeMatrix(X_k, lambda);
  freeMatrix(Z, lambda);
  coco_free_memory(y);
}


/**
 * Euclidian norm
 *
 * @param m pointer to a vector
 * @param dimension size of the vector
 * @return Euclidian norm of the vector
 */
double normE(double *m, size_t dimension)
{
  double out = 0;
  for (size_t i = 0; i < dimension; i++)
  {
    out += m[i]*m[i];
  }
  return sqrt(out);
}

/**
 * Destructor for handcrafted matrices
 *
 * @param mat pointer to a matrix
 * @param nb_rows Number of rows of the matrix
 */
void freeMatrix(double** mat, size_t nb_rows)
{
  for (size_t i = 0; i < nb_rows; i++ )
  {
     delete mat[i];
  }
  delete mat;
}
