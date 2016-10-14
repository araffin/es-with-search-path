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
  double X[dimension]; // To be randomely initialized - double (not float) to be used with the function evaluate
  double X_tmp[dimension]; // Used when sorting the solutions
  double Z_tmp[dimension]; // Used when sorting the solutions
  double Sigma[dimension]; // idem
  double s_sigma[dimension]; // search path --> vector !
  bool happy = false;
  int counter = 0;
  size_t lambda = 40;
  size_t mu = (size_t) lambda/4;
  // Damping factors
  double d = 1. + sqrt((double)mu/(double)dimension); 
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
  
  // double *x = coco_allocate_vector(dimension);
  double y[number_of_objectives];

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
    
    // Select the mu best (shell sort)
    // Shell sort : better when using greater lambda
    int gaps[8] = {701, 301, 132, 57, 23, 10, 4, 1};

    // Start with the largest gap and work down to a gap of 1
    for (size_t g = 0; g < 8; g++)
    {
      size_t gap = gaps[g];
      // Do a gapped insertion sort for this gap size.
      // The first gap elements fitness[0..gap-1] are already in gapped order
      // keep adding one more element until the entire array is gap sorted
      for (size_t k = gap; k < lambda; k += 1)
      {
        double y_tmp = fitness[k];
        // X_tmp = X_k[k]
        vectorCopy(X_k[k], X_tmp, dimension);
        // Z_tmp = Z[k]
        vectorCopy(Z[k], Z_tmp, dimension);
        size_t pos = k;
        // shift earlier gap-sorted elements up 
        // until the correct location for fitness[k] is found
        while(pos >= gap && fitness[pos-gap] > y_tmp)
        {
          fitness[pos] = fitness[pos-gap];
          vectorCopy(X_k[pos-gap], X_k[pos], dimension);
          vectorCopy(Z[pos-gap], Z[pos], dimension);
          pos -= gap;
        }
        // Repostion f(x_k)
        fitness[pos] = y_tmp;
        // Reposition x_k : X_k[pos] = X_tmp
        vectorCopy(X_tmp, X_k[pos], dimension);
        // Reposition z_k : Z[pos] = Z_tmp
        vectorCopy(Z_tmp, Z[pos], dimension);
      }
    }

    // Select the mu best (insertion sort)
    // Better when lambda small
    // for (size_t k = 0; k < lambda-1; k++)
    // {
    //   double y_tmp = fitness[k];
    //   // X_tmp = X_k[k]
    //   vectorCopy(X_k[k], X_tmp, dimension);
    //   // Z_tmp = Z[k]
    //   vectorCopy(Z[k], Z_tmp, dimension);
    //   
    //   size_t pos = k;
    //   while (pos > 0 && fitness[pos-1] > y_tmp)
    //   {
    //     fitness[pos] = fitness[pos-1];
    //     vectorCopy(X_k[pos-1], X_k[pos], dimension);
    //     vectorCopy(Z[pos-1], Z[pos], dimension);
    //     pos--;
    //   }
    //   // Repostion f(x_k)
    //   fitness[pos] = y_tmp;
    //   // Reposition x_k : X_k[pos] = X_tmp
    //   vectorCopy(X_tmp, X_k[pos], dimension);
    //   // Reposition z_k : Z[pos] = Z_tmp
    //   vectorCopy(Z_tmp, Z[pos], dimension);
    // }

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
  evaluate(X,y);
  printArray(X,dimension);
  printArray(y,number_of_objectives);
  // Free memory
  freeMatrix(X_k, lambda);
  freeMatrix(Z, lambda);
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

void euclidianNorm(double* x, double* y)
{
  y[0] = normE(x, 40);
}

int main(int argc, char const *argv[]) {
  size_t dim = 40;
  size_t number_of_objectives = 1;
  size_t budget = 1000;
  double lower[dim];
  double upper[dim];
  for (size_t j = 0; j < dim; j++) {
    lower[j] = -10;
    upper[j] = 10;
  }
  algo4(euclidianNorm, dim, number_of_objectives, lower, upper, budget);
  return 0;
}