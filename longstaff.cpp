// Tarek Frahi

#include <algorithm>
#include <set>
#include <iterator>
#include <vector>
#include <numeric>
#include <time.h>
#include <cmath>
#include "matrix.hpp"
#include "normale.hpp"
#include "longstaff.hpp"
using namespace std;

class Longstaff {
  
    public:
		
		Longstaff(void){}
		~Longstaff(void){}
		
        vector<double> AmericanLongstaff(	double S0,
						double K, 
						double r,
						double q,
						double v, 
						double T,
						int NT, 
						int NS,
						string PutCall) {
            
            CMatrixOperations Mat;
			Normal normal;
			
            // Time increment
            double dt = T/double(NT);

			// Generate Black Scholes stock price paths
			srand (time(0));
            double Z;
            vector<vector<double> > S(NT,vector<double> (NS));
            for (int s=0; s<=NS-1; s++) {
                for (int t=0; t<=NT-1; t++) {
                    Z = normal.normalRandom();
                    if (t==0)
                        S[t][s] = S0;
                    else {
                        S[t][s] = S[t-1][s]*exp((r-q-0.5*v*v)*dt + v*sqrt(dt)*Z);
                    }
                }
            }

            // Initialize the Cash Flows.
            vector<vector<double> > CF(NS,vector<double> (NT));
            
            // Set the last cash flows to the intrinsic value.
            for(int s=0;s<=NS-1;s++)
                if(PutCall == "Put")
                    CF[s][NT-1] = max(K - S[s][NT-1],0.0);
                else if(PutCall == "Call")
                    CF[s][NT-1] = max(S[s][NT-1] - K,0.0);
            
            // European price
            double EuroPrice = 0.0;
            for(int s=0;s<=NS-1;s++)
                EuroPrice += exp(-r*T)*CF[s][NT-1]/double(NS);
            
            // Work backwards through the stock prices until time t=2.
            // We could work through to time t=1 but the regression will not be
            // of full rank at time 1, so this is cleaner.
            for(int t=NT-2;t>=1;t--)
            {
                // Indices for stock paths in-the-money at time t
                vector<int> I(NS);
                for(int s=0;s<=NS-1;s++)
                {
                    I[s] = 0;
                    if(((PutCall == "Put") & (S[s][t] < K)) | ((PutCall == "Call") & (S[s][t] > K)))
                        I[s] = 1;
                }
                
                // Stock paths in-the-money at time t
                int NI = 0;
                vector<double> X;
                vector<int> Xi;
                for(int s=0;s<=NS-1;s++)
                    if(I[s] == 1) {
                        X.push_back(S[s][t]);
                        Xi.push_back(s);
                        NI += 1;
                    }
                
                // Cash flows at time t+1, discounted one period
                vector<double> Y(NI);
                for(int s=0;s<=NI-1;s++)
                    Y[s] = CF[Xi[s]][t+1]*exp(-r*dt);
                
                // Design matrix for regression to predict cash flows
                vector<vector<double> > Z(NI,vector<double> (3));
                for(int s=0;s<=NI-1;s++) {
                    Z[s][0] = 1.0;
                    Z[s][1] = (1.0 - X[s]);
                    Z[s][2] = (2.0 - 4.0*X[s] - X[s]*X[s])/2.0;
                }
                
                // Regression parameters and predicted cash flows
                vector<double> beta = Mat.betavec(Z,Y);
                vector<double> PredCF = Mat.MVecMult(Z,beta);
                
                // Indices for stock paths where immediate exercise is optimal
                // J[s] contains the path number
                vector<int> E;
                vector<int> J;
                int NE = 0;
                for(int s=0;s<=NI-1;s++)
                    if(((PutCall == "Put") & (K - X[s]>0)  &  (K - X[s]> PredCF[s])) |
                       ((PutCall == "Call") & (X[s] - K>0)  &  (X[s] - K> PredCF[s])))
                    {
                        J.push_back(s);
                        E.push_back(Xi[s]);
                        NE += 1;
                    }
                
                // All other paths --> Continuation is optimal
                vector<int> All(NS);
                for(int k=0;k<=NS-1;k++)
                    All[k] = k;
                
                // C contains indices for continuation paths
                // C = All - E so that All = (C union E);
                set<int> CoSet;
                set_difference(All.begin(),All.end(),E.begin(),E.end(),inserter(CoSet,CoSet.end()));
                int NC = CoSet.size();
                // Copy the set into a vector
                vector<int> Co(CoSet.begin(),CoSet.end());
                
                // Replace cash flows with exercise value where exercise is optimal
                for(int s=0;s<=NE-1;s++)
                    if(PutCall == "Put")
                        CF[E[s]][t] = max(K - X[J[s]],0.0);
                    else if(PutCall == "Call")
                        CF[E[s]][t] = max(X[J[s]] - K,0.0);
                for(int s=0;s<=NC-1;s++)
                    CF[Co[s]][t] = exp(-r*dt)*CF[Co[s]][t+1];
            }
            
            // Calculate the cash flow at time 2
            vector<double> Op(NS);
            for(int s=0;s<=NS-1;s++)
                Op[s] = CF[s][1];
            
            // Calculate the American price
            double AmerPrice = exp(-r*dt)*Mat.VecMean(Op);
            
            // Return the European and American prices
            vector<double> output(2);
            output[0] = EuroPrice;
            output[1] = AmerPrice;
            return output;
	    return AmerPrice;
        }
};
