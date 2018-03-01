// Tarek Frahi

#include <cmath>
#include "normal.hpp"
using namespace std;

class Normal {

	public:
		
		Normal(void){}
		~Normal(void){}
		
		double N(double x)
		{
			double a1 =  0.254829592;
			double a2 = -0.284496736;
			double a3 =  1.421413741;
			double a4 = -1.453152027;
			double a5 =  1.061405429;
			double p  =  0.3275911;

			int sign = 1;
			if (x < 0)
				sign = -1;
			x = fabs(x)/sqrt(2.0);
			double t = 1.0/(1.0 + p*x);
			double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
			return 0.5*(1.0 + sign*y);
		}

		double uniformRandom()
		{
			return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
		}
		
		double normalRandom()
		{
			double u1=uniformRandom();
			double u2=uniformRandom();
			return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
		}
};