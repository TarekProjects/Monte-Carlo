// Tarek Frahi

#include <cmath>
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include "normal.hpp"
#include "blackscholes.hpp"

using namespace std;

double getInverseCDFValue(double p) {
    
    double a1 = -39.69683028665376;
    double a2 = 220.9460984245205;
    double a3 = -275.9285104469687;
    double a4 = 138.3577518672690;
    double a5 =-30.66479806614716;
    double a6 = 2.506628277459239;
    
    double b1 = -54.47609879822406;
    double b2 = 161.5858368580409;
    double b3 = -155.6989798598866;
    double b4 = 66.80131188771972;
    double b5 = -13.28068155288572;
    
    double c1 = -0.007784894002430293;
    double c2 = -0.3223964580411365;
    double c3 = -2.400758277161838;
    double c4 = -2.549732539343734;
    double c5 = 4.374664141464968;
    double c6 = 2.938163982698783;
    
    double d1 = 0.007784695709041462;
    double d2 = 0.3224671290700398;
    double d3 = 2.445134137142996;
    double d4 = 3.754408661907416;
    
    double p_low =  0.02425;
    double p_high = 1 - p_low;
    long double  q, r, e, u;
    long double x = 0.0;

    if (0 < p && p < p_low) {
        q = sqrt(-2*log(p));
        x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
    }

    if (p_low <= p && p <= p_high) {
        q = p - 0.5;
        r = q*q;
        x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
    }
    
    if (p_high < p && p < 1) {
        q = sqrt(-2*log(1-p));
        x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
    }
    
    if(( 0 < p)&&(p < 1)){
        e = 0.5 * erfc(-x/sqrt(2)) - p;
        u = e * sqrt(2*M_PI) * exp(x*x/2);
        x = x - u/(1 + x*u/2);
    }
    
    return x;
}


double *Stock_Geometric_Wiener(double S0,double mu,double Vol,double T,double u,double Dt){

    double *Result = new double[int(u)];
    Result[0] = S0;
    
    double St = S0;
    double EPS; //EPS as Wiener process
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,6);
    auto rd = std::bind ( distribution, generator );

    for (int i=1; i<u; i++) {
        double r = rd() / 6.0;
        EPS = getInverseCDFValue(r);
        St = St + mu * St * Dt + Vol * St * sqrt(Dt) * EPS;
        Result[i] = St;
    }
    return Result;
};

double *Stock_Arithmetic(double S0,double mu,double Vol,double T,double u,double Dt){
    
    double *Result = new double[int(u)];
    Result[0] = S0;
    
    double St = S0;
    
    for (int i=1; i<u; i++) {
        St = St + 0.1;
        Result[i] = St;
    }
    return Result;
};


double *Stock_Market_Crash(double S0,double mu,double Vol,double T,double u,double Dt){
    
    double *Result = new double[int(u)];
    Result[0] = S0;
    
    double St = S0;
    double EPS;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,6);
    auto rd = std::bind ( distribution, generator );
    
    for (int i=1; i<u; i++) {
        double r = rd() / 6.0;
        EPS = getInverseCDFValue(r);
        St = St - mu * St * Dt - Vol * St * sqrt(Dt) * EPS;
        Result[i] = St;
    }
    return Result;
};


double *Stock_Geometric_Wiener_Volatility(double S0 ,double mu ,
                                  double Vol_0 , double T ,
                                  double u ,double Dt ){
    double *Result = new double[int(u)];
    double St = S0;
    double Vol_t = Vol_0;
    Result[0] = St;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,6);
    auto rd = std::bind ( distribution, generator );

    for (int i=1; i<u; i++) {
        double r = rd() / 6.0;
        St = St + mu * St * Dt + Vol_t * St * getInverseCDFValue(r) * sqrt(Dt);
        Result[i] = St;
        Vol_t = Vol_t + 0.05;
    }
    return Result;
}

double *Stock_Arithmetic_Volatility(double S0 ,double mu ,
                                  double Vol_0 , double T ,
                                  double u ,double Dt ){
    double *Result = new double[int(u)];
    double St = S0;
    double Vol_t = Vol_0;
    Result[0] = St;
    for (int i=1; i<u; i++) {
        St = St + 0.1;
        Result[i] = St;
        Vol_t = Vol_t + 0.05;
    }
    return Result;
}

double *Stock_Market_Crash_Volatility(double S0 ,double mu ,
                                  double Vol_0 , double T ,
                                  double u ,double Dt ){
    double *Result = new double[int(u)];
    double St = S0;
    double Vol_t = Vol_0;
    Result[0] = St;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,6);
    auto rd = std::bind ( distribution, generator );
    
    for (int i=1; i<u; i++) {
        double r = rd() / 6.0;
        St = St - mu * St * Dt - Vol_t * St * getInverseCDFValue(r) * sqrt(Dt);
        Result[i] = St;
        Vol_t = Vol_t + 0.05;
    }
    return Result;
}


int main() {
    
    double S0, mu, Vol, T, Dt, u,r;
    double K_call = 100.0;  // Strike price
    
    S0 	= 96.01; // Underlying price
    T 	= 73.0/365.0; // Time to maturity
    mu 	= 40.0/100.0; // Stock drift
    Vol = 10.0/100.0; // Volatility of the underlying
    Dt 	= 0.001; // Time step
    u 	= 117; // Number of iterations
    r 	=	0.0187; // Risk-free rate (LIBOR 1.87%)
    
    // Deterministic Volatility
    {
        
        double *S_Geo 	= Stock_Geometric_Wiener_Volatility(S0, mu, Vol, T, u, Dt);
        double *S_Ari 	= Stock_Arithmetic_Volatility(S0, mu, Vol, T, u, Dt);
        double *S_Crash = Stock_Market_Crash_Volatility(S0, mu, Vol, T, u, Dt);
        
		Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_Deterministic_Volatility.csv");
        
        file << "Time,Volatility, type, Spot,Price,Delta,Gamma,Vega,Theta, type, Spot,Price,Delta,Gamma,Vega,Theta, type, Spot,Price,Delta,Gamma,Vega,Theta" << endl;
        
        for (int d = 0; d < u; ++d) {
            
            const double s_geo 		= S_Geo[d];
            const double s_ari 		= S_Ari[d];
            const double s_crash 	= S_Crash[d];
            
            call_geo 	= BS.BSPrice(s_geo, K_call, r, 0, Vol, T-Dt*d, "Call");	
            call_ari 	= BS.BSPrice(s_ari, K_call, r, 0, Vol, T-Dt*d, "Call");
            call_crash 	= BS.BSPrice(s_crash, K_call, r, 0, Vol, T-Dt*d, "Call");

            file << d*Dt
            << "," << Vol
            << "," << "Geometric Vol."
            << "," << s_geo
            << "," << call_geo.value
            << "," << call_geo.delta
            << "," << call_geo.gamma
            << "," << call_geo.vega
            << "," << call_geo.theta
            << "," << "Arithmetic Vol."
            << "," << s_ari
            << "," << call_ari.value
            << "," << call_ari.delta
            << "," << call_ari.gamma
            << "," << call_ari.vega
            << "," << call_ari.theta
            << "," << "Market Crash Vol."
            << "," << s_crash
            << "," << call_crash.value
            << "," << call_crash.delta
            << "," << call_crash.gamma
            << "," << call_crash.vega
            << "," << call_crash.theta
            << endl;
            
            Vol += 0.05;
        }
        
        file.close();
    }
    
    // Hedging Frequency: Dt = 0.001
    {
        
        double *S_Geo 	= Stock_Geometric_Wiener_Volatility(S0, mu, Vol, T, u, Dt);
        double *S_Ari 	= Stock_Arithmetic_Volatility(S0, mu, Vol, T, u, Dt);
        double *S_Crash = Stock_Market_Crash_Volatility(S0, mu, Vol, T, u, Dt);
        
        Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_Hedging_Frequency_001.csv");
        
        file << "Time, type, Spot,Price,Delta,Gamma,Vega,Theta, type, Spot,Price,Delta,Gamma,Vega,Theta, type, Spot,Price,Delta,Gamma,Vega,Theta" << endl;
        
        for (int d = 0; d < u; ++d) {
            
            const double s_geo 		= S_Geo[d];
            const double s_ari 		= S_Ari[d];
            const double s_crash 	= S_Crash[d];
            
            call_geo 	= BS.BSPrice(s_geo, K_call, r, 0, Vol, T-Dt*d, "Call");	
            call_ari 	= BS.BSPrice(s_ari, K_call, r, 0, Vol, T-Dt*d, "Call");
            call_crash 	= BS.BSPrice(s_crash, K_call, r, 0, Vol, T-Dt*d, "Call");
            
            file << d*Dt
            << "," << "Geometric Vol."
            << "," << s_geo
            << "," << call_geo.value
            << "," << call_geo.delta
            << "," << call_geo.gamma
            << "," << call_geo.vega
            << "," << call_geo.theta
            << "," << "Arithmetic Vol."
            << "," << s_ari
            << "," << call_ari.value
            << "," << call_ari.delta
            << "," << call_ari.gamma
            << "," << call_ari.vega
            << "," << call_ari.theta
            << "," << "Market Crash Vol."
            << "," << s_crash
            << "," << call_crash.value
            << "," << call_crash.delta
            << "," << call_crash.gamma
            << "," << call_crash.vega
            << "," << call_crash.theta
            << endl;
            
        }
        
        file.close();
    }

    
    // Hedging Frequency: Dt = 0.003
    {
        

        Vol = 20.0/100.0; // Volatility of the underlying
        Dt = 0.003; // Time step
        u = 39; // Number of iterations
        
        double *S_Geo 	= Stock_Geometric_Wiener(S0, mu, Vol, T, u, Dt);
        double *S_Ari 	= Stock_Arithmetic(S0, mu, Vol, T, u, Dt);
        double *S_Crash = Stock_Market_Crash(S0, mu, Vol, T, u, Dt);
        
        Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_Hedging_Frequency_003.csv");
        
        file << "Time, type, Spot,Price,Delta,Gamma,Vega,Theta, type, Spot,Price,Delta,Gamma,Vega,Theta, type, Spot,Price,Delta,Gamma,Vega,Theta" << endl;
        
        for (int d = 0; d < u; ++d) {
            
            const double s_geo 		= S_Geo[d];
            const double s_ari 		= S_Ari[d];
            const double s_crash 	= S_Crash[d];
            
            call_geo 	= BS.BSPrice(s_geo, K_call, r, 0, Vol, T-Dt*d, "Call");	
            call_ari 	= BS.BSPrice(s_ari, K_call, r, 0, Vol, T-Dt*d, "Call");
            call_crash 	= BS.BSPrice(s_crash, K_call, r, 0, Vol, T-Dt*d, "Call");
            
            file << d*Dt
            << "," << "Geometric Vol."
            << "," << s_geo
            << "," << call_geo.value
            << "," << call_geo.delta
            << "," << call_geo.gamma
            << "," << call_geo.vega
            << "," << call_geo.theta
            << "," << "Arithmetic Vol."
            << "," << s_ari
            << "," << call_ari.value
            << "," << call_ari.delta
            << "," << call_ari.gamma
            << "," << call_ari.vega
            << "," << call_ari.theta
            << "," << "Market Crash Vol."
            << "," << s_crash
            << "," << call_crash.value
            << "," << call_crash.delta
            << "," << call_crash.gamma
            << "," << call_crash.vega
            << "," << call_crash.theta
            << endl;
            
        }

        file.close();
    }
    
    // BASF Greeks
    {
        
        Vol = 20.0/100.0; // Volatility of the underlying
        Dt = 0.002; // Time Step
        u = 60; // Number of iterations
        
        double *S_Geo 	= Stock_Geometric_Wiener(S0, mu, Vol, T, u, Dt);
        double *S_Ari 	= Stock_Arithmetic(S0, mu, Vol, T, u, Dt);
        double *S_Crash = Stock_Market_Crash(S0, mu, Vol, T, u, Dt);
        
        Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_Greeks.csv");
        
        file << "Time, type, Spot,Price,Delta,Gamma,Vega,Theta,Rho, type, Spot,Price,Delta,Gamma,Vega,Theta,Rho, type, Spot,Price,Delta,Gamma,Vega,Theta,Rho" << endl;
        
        for (int d = 0; d < u; ++d) {
            
            const double s_geo 		= S_Geo[d];
            const double s_ari 		= S_Ari[d];
            const double s_crash 	= S_Crash[d];
            
            call_geo 	= BS.BSPrice(s_geo, K_call, r, 0, Vol, T-Dt*d, "Call");	
            call_ari 	= BS.BSPrice(s_ari, K_call, r, 0, Vol, T-Dt*d, "Call");
            call_crash 	= BS.BSPrice(s_crash, K_call, r, 0, Vol, T-Dt*d, "Call");
            
            file << d*Dt
            << "," << "Geometric Vol."
            << "," << s_geo
            << "," << call_geo.value
            << "," << call_geo.delta
            << "," << call_geo.gamma
            << "," << call_geo.vega
            << "," << call_geo.theta
            << "," << call_geo.rho
            << "," << "Arithmetic Vol."
            << "," << s_ari
            << "," << call_ari.value
            << "," << call_ari.delta
            << "," << call_ari.gamma
            << "," << call_ari.vega
            << "," << call_ari.theta
            << "," << call_ari.rho
            << "," << "Market Crash Vol."
            << "," << s_crash
            << "," << call_crash.value
            << "," << call_crash.delta
            << "," << call_crash.gamma
            << "," << call_crash.vega
            << "," << call_crash.theta
            << "," << call_crash.rho
            << endl;
            
        }

        
        file.close();
    }
    
    
    // Market crash
    {
        
        mu = 80.0/100.0; // Stock drift
        Vol = 20.0/100.0; // Volatility of the underlying
        Dt = 0.002; // Time steps
        u = 60; // Numbers of iterations
        
        double *S_Geo 	= Stock_Geometric_Wiener(S0, mu, Vol, T, u, Dt);
        double *S_Ari 	= Stock_Arithmetic(S0, mu, Vol, T, u, Dt);
        double *S_Crash = Stock_Market_Crash(S0, mu, Vol, T, u, Dt);
        
        Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_Market_crash.csv");
        
        file << "Time, type, Spot,Price,Delta,Gamma,Vega,Theta,Rho, type, Spot,Price,Delta,Gamma,Vega,Theta,Rho, type, Spot,Price,Delta,Gamma,Vega,Theta,Rho" << endl;
        
        for (int d = 0; d < u; ++d) {
            
            const double s_geo 		= S_Geo[d];
            const double s_ari 		= S_Ari[d];
            const double s_crash 	= S_Crash[d];
            
            call_geo 	= BS.BSPrice(s_geo, K_call, r, 0, Vol, T-Dt*d, "Call");	
            call_ari 	= BS.BSPrice(s_ari, K_call, r, 0, Vol, T-Dt*d, "Call");
            call_crash 	= BS.BSPrice(s_crash, K_call, r, 0, Vol, T-Dt*d, "Call");
            
            file << d*Dt
            << "," << "Geometric Vol."
            << "," << s_geo
            << "," << call_geo.value
            << "," << call_geo.delta
            << "," << call_geo.gamma
            << "," << call_geo.vega
            << "," << call_geo.theta
            << "," << call_geo.rho
            << "," << "Arithmetic Vol."
            << "," << s_ari
            << "," << call_ari.value
            << "," << call_ari.delta
            << "," << call_ari.gamma
            << "," << call_ari.vega
            << "," << call_ari.theta
            << "," << call_ari.rho
            << "," << "Market Crash Vol."
            << "," << s_crash
            << "," << call_crash.value
            << "," << call_crash.delta
            << "," << call_crash.gamma
            << "," << call_crash.vega
            << "," << call_crash.theta
            << "," << call_crash.rho
            << endl;
            
        }
        
        file.close();
    }
    

    // Generation of 2 options traded:
    // Traded Option 1 (Spot 63, Strike 75, Maturity 1.00, Volatility 10%)
    // Traded Option 2 (Spot 43 ,Strike 46, Maturity 0.77, Volatility 20%)
    {
        double Vol_1,T_1;
        T 		= 1.00; 		// Time to maturity, option 1
        mu 		= 40.0/100.0; 	// Stock drift
        Vol 	= 10.0/100.0; 	// Volatility of the underlying, option 1
        Dt 		= 0.002; 		// Time Steps
        u 		= 150; 			// Number of iterations
        Vol_1 	= 0.02; 		// Volatility of the underlying, option 2
        T_1 	= 0.77; 		// Time to maturity, option 2

        double *S 	= Stock_Arithmetic(63.0,mu, Vol, T, u, Dt);
        double *S_1 = Stock_Geometric_Wiener(43.0, mu, Vol_1, T_1, u, Dt);
        
        double K_call = 75.0;  // Strike price option 1
        double K_call_1 = 46.0;  // Strike price option 2
        Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_options_traded.csv");
        
        file << "day, option,Spot,Price,Delta,Gamma,Vega,Theta,Rho, option ,Spot,Price,Delta,Gamma,Vega,Theta,Rho" << endl;
        for (int d = 0; d < u; ++d) {
            
            const double s 		= S[d];
            const double s_1 	= S_1[d];
			
            callBS 		= BS.BSPrice(s, K_call, r, 0, Vol, T-Dt*d, "Call");
            callBS_1 	= BS.BSPrice(s_1, K_call_1, r, 0, Vol_1, T_1-Dt*d, "Call");

            file << d*Dt
            << "," << "option 1"
            << "," << s
            << "," << callBS.value
            << "," << callBS.delta
            << "," << callBS.gamma
            << "," << callBS.vega
            << "," << callBS.theta
            << "," << callBS.rho
            << "," << "option 2"
            << "," << s_1
            << "," << callBS_1.value
            << "," << callBS_1.delta
            << "," << callBS_1.gamma
            << "," << callBS_1.vega
            << "," << callBS_1.theta
            << "," << callBS_1.rho
            << endl;
        }
        
        file.close();
    }

    // Trendy VS Rangy stock
    {

        Vol = 5.0/100.0; // Trendy Volatility of the underlying
        Dt = 0.002; // Time steps
        u = 60; // Number of iterations
        
        double mu_rangy = 10.0/100.0; // Rangy Stock drift
        double Vol_rangy = 80.0/100.0; // Rangy Volatility of the underlying
        
        double *Trendy_ar = Stock_Arithmetic(S0, mu, Vol, T, u, Dt);       
        double *Trendy_geo = Stock_Geometric_Wiener(S0, mu, Vol, T, u, Dt);
        double *Rangy = Stock_Geometric_Wiener(S0, mu, Vol, T, u, Dt);
        
        Blackscholes BS;
        
        ofstream file;
        file.open ("BASF_rangy_trendy.csv");
        
        file << "Time, option,Spot,Price,Delta,Gamma,Theta, option,Spot,Price,Delta,Gamma,Theta, option,Spot,Price,Delta,Gamma,Theta" << endl;
        for (int d = 0; d < u; ++d) {
            
            const double s_t_ar = Trendy_ar[d];
            const double s_t_geo = Trendy_geo[d];
            const double s_r = Rangy[d];

            //Call Trendy Arithmetic
            call_t_ar = BS.BSPrice(s_t_ar, K_call, r, 0, Vol, T-d*Dt, "Call");     
            //Call Trendy Geo
            call_t_geo = BS.BSPrice(s_t_geo, K_call, r, 0, Vol, T-d*Dt, "Call");          
            //Call Rangy Geo
            call_r = BS.BSPrice(s_r, K_call, r, 0, Vol_rangy, T-d*Dt, "Call");
            file << d*Dt
            << "," << "Call Trendy Arithmetic"
            << "," << s_t_ar
            << "," << call_t_ar.value
            << "," << call_t_ar.delta
            << "," << call_t_ar.gamma
            << "," << call_t_ar.theta
            << "," << "Call Trendy Geo"
            << "," << s_t_geo
            << "," << call_t_geo.value
            << "," << call_t_geo.delta
            << "," << call_t_geo.gamma
            << "," << call_t_geo.theta
            << "," << "Call Rangy Geo"
            << "," << s_r
            << "," << call_r.value
            << "," << call_r.delta
            << "," << call_r.gamma
            << "," << call_r.theta
            << endl;
        }
        
        file.close();

    }
    return 0;
}
