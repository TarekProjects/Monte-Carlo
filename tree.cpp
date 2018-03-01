// Tarek Frahi

#include <algorithm>
#include <set>
#include <iterator>
#include <vector>
#include <numeric>
#include <cmath>
#include "matrix.hpp"
#include "tree.hpp"

class Tree {
    
    public:
    
    // Function for binomial tree
		double Binomial(double S0,
						double K, 
						double r,
						double q,
						double v, 
						double T,
						int N,
						string PutCall,
						char EuroAmer) 
		{
			int i,j;
			
			// Matrix for the stock price evolution and option price
			vector<vector<double> > S(N+1,vector<double> (N+1));
			vector<vector<double> > Op(N+1,vector<double> (N+1));

			double dt,u,d,p;
			
			// Quantities for the tree
			dt = T/N;
			u = exp(v*sqrt(dt));
			d = 1.0/u;
			p = (exp(r*dt)-d) / (u-d);
			
			// Build the binomial tree
			for (j=0; j<=N; j++)
				for (i=0; i<=j; i++)
					S[i][j] = S0*pow(u,j-i)*pow(d,i);
			
			// Compute terminal payoffs
			for (i=0; i<=N; i++) {
				if (PutCall=="Call")
					Op[i][N] = max(S[i][N] - K, 0.0);
				else
					Op[i][N] = max(K - S[i][N], 0.0);
			}
			
			// Backward recursion through the tree
			for (j=N-1; j>=0; j--)
				for (i=0; i<=j; i++) {
					if (EuroAmer=='E')
						Op[i][j] = exp(-r*dt)*(p*(Op[i][j+1]) + (1.0-p)*(Op[i+1][j+1]));
					else {
						if (PutCall=="Call")
							Op[i][j] = max(S[i][j] - K, exp(-r*dt)*(p*(Op[i][j+1]) + (1.0-p)*(Op[i+1][j+1])));
						else
							Op[i][j] = max(K - S[i][j], exp(-r*dt)*(p*(Op[i][j+1]) + (1.0-p)*(Op[i+1][j+1])));
					}
				}
			// Return the option price
			return Op[0][0];
		}
	
	
		// Function for trinomial tree
        double Trinomial(	double S0,
							double K, 
							double r,
							double q,
							double v, 
							double T,
							int N,
							string PutCall,
							char EuroAmer)
        {
            int i,j;
            double b = r - q;
            double dt = T/N;
            double u = exp(v*sqrt(2.0*dt));
            double pu = pow((exp(0.5*b*dt) - exp(-v*sqrt(0.5*dt))) / (exp(v*sqrt(0.5*dt)) - exp(-v*sqrt(0.5*dt))), 2);
            double pd = pow((exp(v*sqrt(0.5*dt)) - exp(0.5*b*dt))  / (exp(v*sqrt(0.5*dt)) - exp(-v*sqrt(0.5*dt))), 2);
            double pm = 1.0 - pu - pd;
            
            // Build the trinomial tree for the stock price evolution
            vector<vector<double> > S(2*N+1,vector<double> (N+1));
            for (j=0; j<=N; j++)
                for (i=0; i<=2*j; i++)
                    S[i][j] = S0*pow(u,double(j-i));
            
            // Boost Matrices for the Euro and American calls and puts
            vector<vector<double> > V(2*N+1,vector<double> (N+1));
            
            // Compute terminal payoffs
            for (i=0; i<=2*N; i++) {
                if (PutCall == "Call")
                    V[i][N] = max(S[i][N] - K, 0.0);
                else if (PutCall == "Put")
                    V[i][N] = max(K - S[i][N], 0.0);
            }
            
            // Backward recursion through the tree
            for (j=N-1; j>=0; j--) {
                for (i=0; i<=2*j; i++) {
                    if (EuroAmer == 'E') {
                        V[i][j] = exp(-r*dt)*(pu*V[i][j+1] + pm*V[i+1][j+1] + pd*V[i+2][j+1]);
                    }
                    else if (EuroAmer == 'A') {
                        if (PutCall == "Call")
                            V[i][j] = max(S[i][j] - K, exp(-r*dt)*(pu*V[i][j+1] + pm*V[i+1][j+1] + pd*V[i+2][j+1]));
                        else if (PutCall  == "Put")
                            V[i][j] = max(K - S[i][j], exp(-r*dt)*(pu*V[i][j+1] + pm*V[i+1][j+1] + pd*V[i+2][j+1]));
                    }
                }
            }
            // Output the results
            return V[0][0];
        }	
};
