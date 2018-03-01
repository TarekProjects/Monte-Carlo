// Tarek Frahi

#include <cmath>
#include <numeric>
#include "normal.hpp"
using namespace std;

 typedef struct{
            double value;
            double delta;
            double gamma;
            double vega;
            double theta;
            double rho;
            double vanna;
            double vomma;
            double speed;
            double zomma;
        }BS;

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
};