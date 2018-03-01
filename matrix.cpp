// Tarek Frahi

#include <cmath>
#include <vector>
#include <numeric>
#include "matrix.hpp"
using namespace std;

CMatrixOperations::CMatrixOperations(void) { }
CMatrixOperations::~CMatrixOperations(void) { }

// LU decomposition;
LUstruct CMatrixOperations::LU(vector<vector<double> > A) {
    int N = A.size();
    vector<vector<double> >  B(N,vector<double> (N));
    for(int i=0;i<=N-1;i++)
        for(int j=0;j<=N-1;j++)
            B[i][j] = A[i][j];
    
    for(int k=0;k<=N-2;k++) {
        for(int i=k+1;i<=N-1;i++)
            B[i][k] = B[i][k] / B[k][k];
        for(int j=k+1;j<=N-1;j++)
            for(int i=k+1;i<=N-1;i++)
                B[i][j] = B[i][j] - B[i][k]*B[k][j];
    }
    vector<vector<double> >  L(N,vector<double> (N));
    vector<vector<double> >  U(N,vector<double> (N));
    for(int i=0;i<=N-1;i++)	{
        L[i][i] = 1.0;
        for(int j=0;j<=N-1;j++) {
            if(i>j)
                L[i][j] = B[i][j];
            else
                U[i][j] = B[i][j];
        }
    }
    LUstruct Mats;
    Mats.L = L;
    Mats.U = U;
    return Mats;
}

// Inverse of an upper triangular matrix
vector<vector<double> >  CMatrixOperations::MatUpTriangleInv(vector<vector<double> > U) {
    int N = U.size();
    vector<vector<double> >  V(N,vector<double> (N));
    for(int j=N-1;j>=0;j--) {
        V[j][j] = 1.0 / U[j][j];
        for(int i=j-1;i>=0;i--)
            for(int k=i+1;k<=j;k++)
                V[i][j] -= 1.0 / U[i][i] * U[i][k] * V[k][j];
    }
    return V;
}

// Inverse of a lower triangular matrix
vector<vector<double> >  CMatrixOperations::MatLowTriangleInv(vector<vector<double> > L) {
    int N = L.size();
    vector<vector<double> >  V(N,vector<double> (N));
    for(int i=0;i<=N-1;i++)	{
        V[i][i] = 1.0 / L[i][i];
        for(int j=i-1;j>=0;j--)
            for(int k=i-1;k>=j;k--)
                V[i][j] -= 1.0 / L[i][i] * L[i][k] * V[k][j];
    }
    return V;
}

// Multiply two matrices together
// First matrix is (n x k), Second is (k x m)
// Resulting matrix is (n x m)
vector<vector<double> > CMatrixOperations::MMult(vector<vector<double> >A,vector<vector<double> >B, int n, int k, int m) {
    vector<vector<double> >  C(n, vector<double> (m));
    for (int j=0; j<=m-1; j++)
        for (int i=0; i<=n-1; i++) {
            C[i][j] = 0;
            for (int r=0; r<=k-1; r++)
                C[i][j] += A[i][r] * B[r][j];
        }
    return C;
}

// Inverse of a matrix through LU decomposition
vector<vector<double> >  CMatrixOperations::MInvLU(vector<vector<double> > A) {
    LUstruct Mats;
    Mats = LU(A);
    vector<vector<double> >  L = Mats.L;
    vector<vector<double> >  U = Mats.U;
    vector<vector<double> >  Linv = MatLowTriangleInv(L);
    vector<vector<double> >  Uinv = MatUpTriangleInv(U);
    int n = Linv.size();
    int k = Linv[0].size();
    int m = Uinv[0].size();
    vector<vector<double> >  Ainv = MMult(Uinv,Linv,n,k,m);
    return Ainv;
}
// Transpose of a matrix
vector<vector<double> > CMatrixOperations::MTrans(vector<vector<double> >A) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<double> >  C(cols,vector<double> (rows));
    for (int i=0; i<cols; i++)
        for (int j=0; j<rows; j++)
            C[i][j] = A[j][i];
    return C;
}

// Multiply a matrix by a vector
vector<double> CMatrixOperations::MVecMult(vector<vector<double> > X, vector<double> Y) {
    int rows = X.size();
    int cols = X[0].size();
    vector<double> XY(rows);
    for (int r=0; r<rows; r++) {
        XY[r] = 0.0;
        for (int c=0; c<cols; c++)
            XY[r] += X[r][c] * Y[c];
    }
    return XY;
}

// Regression beta coefficients
vector<double> CMatrixOperations::betavec(vector<vector<double> > X, vector<double> Y) {
    int rows = X.size();
    int cols = X[0].size();
    vector<vector<double> > Xt = MTrans(X);
    vector<vector<double> > XtX = MMult(Xt,X,cols,rows,cols);
    vector<vector<double> > XtXinv = MInvLU(XtX);
    vector<double> XtY = MVecMult(Xt,Y);
    return MVecMult(XtXinv,XtY);
}

// Mean of a vector
double CMatrixOperations::VecMean(vector<double> X) {
    double vectorsum = VecSum(X);
    int n = X.size();
    return vectorsum/n;
}

// Sum of a vector
double CMatrixOperations::VecSum(vector<double> X) {
    double initialize = 0.0;
    double vectorsum = accumulate(X.begin(),X.end(),initialize);
    return vectorsum;
}