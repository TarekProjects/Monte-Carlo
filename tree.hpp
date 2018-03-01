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
						char EuroAmer); 
	
		// Function for trinomial tree
        double Trinomial(	double S0,
							double K, 
							double r,
							double q,
							double v, 
							double T,
							int N,
							string PutCall,
							char EuroAmer);
	
};


