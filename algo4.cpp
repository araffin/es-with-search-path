// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
// Post Process : python2 -m bbob_pproc [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
#include "algo4.hpp"
#include <assert.h>
#define PI 3.14159265359
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

  double X[dimension]; // To be randomely initialized - double (not float) to be used with the function evaluate
  double X_tmp[dimension]; // Used when sorting the solutions
  double Z_tmp[dimension]; // Used when sorting the solutions
  double Sigma[dimension]; // idem
  double s_sigma[dimension]; // search path --> vector !
  bool happy = false;
  int counter = 0;
  size_t lambda = 20;
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
  uniform_real_distribution<double> random_uniform(0.0,1.0);
  
  for (size_t j = 0; j < dimension; j++)
  {
    X[j] = lower_bounds[j] + (upper_bounds[j] - lower_bounds[j])*random_uniform(gen);
  }
  // a t-on besoin d'un critère d'arrêt ? Sinon on fait juste le nb d'itérations comme dans les exemples
  // double stop_criterion = 0.0002; // juste pour tester les 3 dernières lignes de la boucle

  // double *x = coco_allocate_vector(dimension);
  double *y = coco_allocate_vector(number_of_objectives);

  for (size_t j = 0; j < dimension ; j++)
  {
    Sigma[j] = 0.1;
  }

  while (!happy && counter < max_budget)
  {
    for (size_t k = 0; k < lambda; k++)
    {
      // z_k = N(0,I)
      // x_k = x + sigma o z_k
      for (size_t j = 0; j < dimension; j++) {
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
    // Select the mu best (insertion sort)
    for (size_t k = 0; k < lambda-1; k++)
    {
      double y_tmp = fitness[k];
      // X_tmp = X_k[k]
      vectorCopy(X_k[k], X_tmp, dimension);
      // Z_tmp = Z[k]
      vectorCopy(Z[k], Z_tmp, dimension);
      
      size_t pos = k;
      while (pos > 0 && fitness[pos-1] > y_tmp)
      {
        fitness[pos] = fitness[pos-1];
        vectorCopy(X_k[pos-1], X_k[pos], dimension);
        vectorCopy(Z[pos-1], Z[pos], dimension);
        pos--;
      }
      // Repostion f(x_k)
      fitness[pos] = y_tmp;
      // Reposition x_k : X_k[pos] = X_tmp
      vectorCopy(X_tmp, X_k[pos], dimension);
      // Reposition z_k : Z[pos] = Z_tmp
      vectorCopy(Z_tmp, Z[pos], dimension);
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
    // happy = (y[0] < stop_criterion);
    counter++;
  }
  // Free memory
  freeMatrix(X_k, lambda);
  freeMatrix(Z, lambda);
  coco_free_memory(y);
}

double normE(double *m, size_t dimension)
{
  double out = 0;
  for (size_t i = 0; i < dimension; i++)
  {
    out += m[i]*m[i];
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

void freeMatrix(double** mat, size_t nb_rows)
{
  for (size_t i = 0; i < nb_rows; i++ )
  {
     delete mat[i];
  }
  delete mat;
}