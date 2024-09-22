#include <iostream>
#include <vector>
#include <cmath>
#include "fMatrix.h"
using namespace std;

int main()
{
    Float A[3] = {1.1, 2.2, 3.3};
    Float B[3] = {4.4, 5.5, 6.6};
    Float Z[9] = {1,0,2,-1,5,0,0,3,-9};


    Float M[9] = {1, 1, 3,
                    4 ,0, 4, 
                    5, 0, 5};

    Float R[9] = {5.0990 ,6.2757, 0.5883 ,
                    0.0000, 1.2710, 4.9629 ,
                    0.0000 ,0.0000, 0.1543}; 

    Float Q[9] = {0.1961, 0.6052,-0.7715, 
                    0.0000, 0.7868, 0.6172, 
                    0.9806, -0.1210, 0.1543};
	Float X[15] = {0.7922, 0.9340, 0.6555, 
				   0.9595, 0.6787, 0.1712, 
				   0.6557, 0.7577, 0.7060, 
				   0.0357, 0.7431, 0.0318, 
				   0.8491, 0.3922, 0.2769};
    fVector VecA(A, 3);
    fVector VecB(B, 3);
    fVector VecC(VecA);
    fMatrix MatA(Q, 3, 3);
    fMatrix MatB(R, 3, 3);
    fMatrix MatC(MatA);
    fMatrix MatZ(Z, 3, 3);
    fMatrix MatM(M, 3, 3);
    // Transp(MatZ).Show();
    // AATransp(MatZ).Show();
    // ATranspA(MatZ).Show();
    Inverse(MatZ).Show();
}