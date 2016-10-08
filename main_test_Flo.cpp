
#include "stdio.h"
#include <math.h>
#include <iostream>

using namespace std;



int main()
{

size_t dimension = 3;

int lambda = 4; // lambda is given, how do we have it ? 
int mu = (int) lambda/4;

cout << "lambda,mu = " << lambda << " , " << mu << endl; 
double c_sigma = sqrt((double)mu/(dimension + mu)); 
cout << "c_sigma = " << c_sigma << endl; 
double s_sigma[dimension] = {0}; // search path --> vector !
double X[dimension] = {0};
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
        Z[i][j] = 0.1;
     }
  }

  cout << "X_k = " << X_k << endl; 

  cout << "Z = " << Z << endl;

  // uptdate s_sigma // TO BE CHECKED
  for(int i = 0; i< dimension; i ++)
  {
    double sum = 0;  // sum = sum(zk); zk in P
    for(int j = 0; j < mu ;  j++) // we take the mu best in Z
    {
      sum = sum + Z[j][i]; 
      cout << Z[j][i] << endl;
      cout << "sum = " << sum << endl; 
    }

    s_sigma[i] = (1-c_sigma)*s_sigma[i] + sqrt(c_sigma*(2 - c_sigma)*mu)/mu * sum;
    cout << "s_sigma[" << i << "] = " << s_sigma[i] << endl; 
  }


    // update X
    for(int i = 0 ; i < dimension ; i ++)
    {
      double sumXk = 0; 
      for (int j = 0 ; j < mu ; j++)
      {
        sumXk = sumXk + X_k[j][i];
      }
      X[i] = (double)sumXk/mu;
      cout << X[i] << endl;
    }


  cout << "s_sigma = " << s_sigma << endl;



}