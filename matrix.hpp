// Tarek Frahi

#include <cmath>
#include <vector>
#include <numeric>
using namespace std;

typedef struct{
    vector<vector<double> > L;
    vector<vector<double> > U;
}LUstruct;


class CMatrixOperations {
public:
    CMatrixOperations(void);
    ~CMatrixOperations(void);
    
    // LU decomposition;
    LUstruct LU(vector<vector<double> > A);
    
    // Inverse of an upper triangular matrix
    vector<vector<double> >  MatUpTriangleInv(vector<vector<double> > U);
    
    // Inverse of a lower triangular matrix
    vector<vector<double> >  MatLowTriangleInv(vector<vector<double> > L);
    
    // Multiply two matrices together
    // First matrix is (n x k), Second is (k x m)
    // Resulting matrix is (n x m)
    vector<vector<double> > MMult(vector<vector<double> >A,vector<vector<double> >B, int n, int k, int m);
    
    // Inverse of a matrix through LU decomposition
    vector<vector<double> >  MInvLU(vector<vector<double> > A);
    
    // Transpose of a matrix
    vector<vector<double> > MTrans(vector<vector<double> >A);
    
    // Multiply a matrix by a vector
    vector<double> MVecMult(vector<vector<double> > X, vector<double> Y);
    
    // Regression beta coefficients
    vector<double> betavec(vector<vector<double> > X, vector<double> Y);
    
    // Mean of a vector
    double VecMean(vector<double> X);
    
    // Sum of a vector
    double VecSum(vector<double> X);
};