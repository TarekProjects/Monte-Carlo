// Tarek Frahi

#include <cmath>
#include <numeric>
#include "normal.hpp"
#include "blackscholes.hpp"
using namespace std;

class BlackScholes {
    
    public:
		
		BlackScholes(void){}
		~BlackScholes(void){}

		BS BSPrice(	double S0,
						double K, 
						double r,
						double q,
						double v, 
						double T,
						string PutCall)
        {

			double d1 = (log(S0/K) + (r+pow(v,2)/2)*T)/(v*sqrt(T));
            double d2 = d1 - v*sqrt(T);
            double D1 = normal::N(d1);
            double D2 = normal::N(d2);
            double phi = exp(-pow(d1,2)/2)/sqrt(2*M_PI);
			
			if (PutCall=="Call") {
				double value = D1*S0 - D2*K*exp(-r*T);
				double delta = D1;         
				double gamma = phi/(S0*v*sqrt(T));
				double vega = S0*phi*sqrt(T);
				double theta = -S0*phi*v/(2*sqrt(T)) - r*K*D2*exp(-r*T);
				double rho = K*T*D2*exp(-r*T);   
				double zomma =  (phi*(d1*d2 - 1))/(S0*v*v*sqrt(T));
				double speed =  (-phi * ( (d1/(v*sqrt(T))) + 1 )) / (pow(S0, 2)*v*sqrt(T));
				double vanna =  phi*d2/v;
				double vomma =  S0*phi*sqrt(T)*d1*d2/v;
                return BS{value, delta, gamma, vega, theta, rho, vanna, vomma, speed, zomma};
				//return value;
			}
			
			if (PutCall=="Put") {
				double value =  D2*K*exp(-r*T) - D1*S0;
				double delta = - D1;
				double gamma = phi/(S0*v*sqrt(T));
				double vega = S0*phi*sqrt(T);
				double theta = -S0*phi*v/(2*sqrt(T)) + r*K*D2*exp(-r*T);
				double rho = -K*T*D2*exp(-r*T);
				double zomma =  (phi*(d1*d2 - 1))/(S0*v*v*sqrt(T));
				double speed =  (-phi * ( (d1/(v*sqrt(T))) + 1 )) / (pow(S0, 2)*v*sqrt(T));
				double vanna =  phi*d2/v;
				double vomma =  S0*phi*sqrt(T)*d1*d2/v;
                return BS{value, delta, gamma, vega, theta, rho, vanna, vomma, speed, zomma};
				//return value;
			}
        }
};