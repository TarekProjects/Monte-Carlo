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
using namespace std;

class Longstaff {
  
    public:
		
		Longstaff(void);
		~Longstaff(void);
		
		vector<double> AmericanLongstaff(	double S0,
							double K, 
							double r,
							double q,
							double v, 
							double T,
							int NT, 
							int NS,
							string PutCall);
};
