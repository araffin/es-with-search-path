// Evolution Strategy with search path
// as describe in https://www.lri.fr/~hansen/es-overview-2015.pdf [Algo 4]
#include <iostream>
#include <math.h>
#include <random> // Normal distribution

using namespace std;

void printArray(float* array, int n)
{
  for (size_t i = 0; i < n; i++)
    cout << array[i] << " ";
  
  cout << endl;
}

void printMatrix(float** array, int n, int m)
{
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < m; j++) {
      cout << array[i][j] << " ";
    }
    cout << endl;
  }
}

void elementProduct(float* a, float* b, float* result, int n)
{
  for (size_t i = 0; i < n; i++)
    result[i] = a[i]*b[i];
}
/**
 * Sum of two arrays
 * @param a      [description]
 * @param b      [description]
 * @param result [description]
 * @param n      [description]
 */
void arraySum(float* a, float* b, float* result, int n)
{
  for (size_t i = 0; i < n; i++)
    result[i] = a[i] + b[i];
}

/**
 * Sum of k vectors (of dimension n)
 * For each column sum the k rows of the matrix a 
 * @param a      matrix (array of array)
 * @param result [description]
 * @param k      [description]
 * @param n      [description]
 */
void sumVectors(float** a, float* result, int k, int n)
{
 for (size_t j = 0; j < n; j++) {
   result[j] = 0;
   for (size_t i = 0; i < k; i++) {
     result[j] += a[i][j];
   }
 }
}

// Or something like :
// template <size_t rows, size_t cols>
// void process_2d_array_template(int (&array)[rows][cols])
// for 2D arrays ?
void normalVector(float** Z, normal_distribution<> N, mt19937 gen, int lambda, int n)
{
  for (size_t k = 0; k < lambda; k++) {
    for (size_t j = 0; j < n; j++) {
      Z[k][j] = N(gen());
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
      
    );
    
    counter++;
  }
  std::cout << "/* Best Solution :  */" << printArray(X_best, n) << std::endl;  
}

int main(int argc, char const *argv[])
{
  float **mat;
  int rows = 2, cols = 4, i, j;
  // dynamically allocate memory for a matrix of size m x n
  mat = new float*[rows];
  for ( i = 0; i < rows; i++ ) {
     mat[i] = new float[cols];
  }
  // initialize matrix elements
  for ( i = 0; i < rows; i++ ) {
     for ( j = 0; j < cols; j++ ) {
        mat[i][j] = i + j;         
     }
  }
  // display matrix 
  printMatrix(mat, rows, cols);

  // deallocate (free) the memory
  for ( i = 0; i < rows; i++ ) {
     delete mat[i];
  }
  delete mat;
  return 0;
}