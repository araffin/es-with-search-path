// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
// Post Process : python -m bbob_pproc [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
#include <iostream>
#include <math.h>
#include <random> // Normal distribution
#include <assert.h>
#include "test_algo.h"

using namespace std;

/* Selection function
So far I have considered that this function will classify the best member of X_k and Z at the fist positions
ie P will be the mu first columns of X_k and Z (to avoid having new matrix) but I don't know if it is the best way

*/

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


#define PI 3.14159265359

void algo4(void (*evaluate)(double*, double*),
                    const size_t dimension,
                    const size_t number_of_objectives,
                    double *lower_bounds,
                    double *upper_bounds,
                    const size_t max_budget)
{

  // ======================================================
  // INITIALIZATION
  // ======================================================

  assert(number_of_objectives == 1);

  // double stop_criterion = 0.00001; 
  double X[dimension] = {0}; // To be randomely initialized - double (not float) to be used with the function evaluate
  double X_tmp[dimension] = {0}; // Used when sorting the solutions
  double Z_tmp[dimension] = {0}; // Used when sorting the solutions
  double Sigma[dimension]; // idem
  double s_sigma[dimension] = {0}; // search path --> vector !
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
  double _d = double(dimension);
  double E_MULTIDIM_NORMAL = sqrt(_d)*(1-1./(4*_d)+1./(21*_d*_d));
  // Matrix Initialization
  X_k = new double*[lambda];
  Z = new double*[lambda];
  for (size_t k = 0; k < lambda; k++)
  {
    X_k[k] = new double[dimension];
    Z[k] = new double[dimension];
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
  
  // initialize matrix elements
  for (size_t k = 0; k < lambda; k++)
  {
    for (size_t j = 0; j < dimension; j++)
    {
      X_k[k][j] = lower_bounds[j] + (upper_bounds[j] - lower_bounds[j])*random_uniform(gen);
    }
  }
  // double *x = coco_allocate_vector(dimension);
  double y[number_of_objectives];

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
      // printArray(X_k[k],dimension);

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
      // printArray(X_k[k],dimension);
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
    double exp1 = 0;
    double exp2 = exp((normE(s_sigma, dimension)/E_MULTIDIM_NORMAL) -1);
    exp2 = pow(exp2, c_sigma/d);
    for(size_t j = 0; j < dimension; j++)
    {
      exp1 = exp((abs(Sigma[j])/E_HALF_NORMAL -1));
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
      X[j] = fmax(lower_bounds[j], fmin(X[j], upper_bounds[j]));
    }

    // utilise t-on vraiment un critère d'arrêt ?
    evaluate(X,y);
    // happy = (y[0] < stop_criterion);
    // printArray(y,number_of_objectives);
    // printArray(X,dimension);

    counter++;
  }
  evaluate(X,y);
  printArray(X,dimension);
  printArray(y,number_of_objectives);
  // Free memory
  freeMatrix(X_k, lambda);
  freeMatrix(Z, lambda);
  // coco_free_memory(y);
}

void poly(double* x, double* y)
{
  y[0] = 4*x[0]*x[0] + 2*x[0] + 5;
}

void poly2(double* x, double* y)
{
  y[0] = 18*x[0]*x[0] + 2*x[0] -20;
}

void poly3(double* x, double* y)
{
  y[0] = 4*x[1]*x[1] + 2*x[1] + 5;
}

void twoDim(double* x, double* y)
{
  double a = 4*x[1]*x[1] + 2*x[1] + 5;
  double b = 18*x[0]*x[0] + 2*x[0] - 20;
  y[0] = a + b;
}

int main(int argc, char const *argv[]) {
  size_t dim = 2;
  size_t number_of_objectives = 1;
  size_t budget = 100;
  double lower[dim];
  double upper[dim];
  for (size_t j = 0; j < dim; j++) {
    lower[j] = -10;
    upper[j] = 10;
  }
  algo4(twoDim, dim, number_of_objectives, lower, upper, budget);
  return 0;
}