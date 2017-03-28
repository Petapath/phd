#define _USE_MATH_DEFINES
 
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <limits>

const double h = 6.626176e-34;
const double c = 2.997925e8;
const double k = 1.380662e-23;


double toJoule( double eV ){
	return eV*(1.60218e-19);
}

double T1(int n, double x){
	return ((x*x + 2*(1/n + x)/n)/n)*exp(-n*x);
}

double T2(int n, double x){
	return ((1/n + x)/n)*exp(-n*x);
}

double T3(int n, double x){
	return exp(-n*x);
}

double integral_e2inf( double E_eV, double T_K, double u_eV ){
	double x=toJoule(E_eV-u_eV)/(k*T_K);
	int iter=std::min(20.0+20.0/x, 512.0);
	//printf( "-------------- %f %d\n", x, iter);
	
	double t1=0.0,t2=0.0,t3=0.0;
	
	for( int i=1;i<iter;i++ ) {
		t1+=T1(i,x);
		t2+=T2(i,x);
		t3+=T3(i,x);
	}
	
	return ( 2*pow(k,3)*pow(T_K,3)*t1 + 
	         4*toJoule(u_eV)*pow(k,2)*pow(T_K,2)*t2 +
		 2*k*T_K*pow(toJoule(u_eV),2)*t3 ) / (pow(c,2) * pow(h,3));
}

// flux units: photons s-1 m-2 sr-1
double N( double E1_eV, double E2_eV, double T_K, double u_eV, double etendue ) {
	double r1 = integral_e2inf( E1_eV, T_K, u_eV );
	double r2 = 0.0;

	if( E2_eV != std::numeric_limits<double>::infinity() ) 
		r2 = integral_e2inf( E2_eV, T_K, u_eV );
	
	return etendue*(r1-r2);
}

// Power = radiance L * etendue, where radiance units: W m-2 sr-1 and etendue sr
double P( double T_K, double etendue ) {
	return (2*pow(M_PI*k*T_K,4))/(15*pow(h,3)*pow(c,2))*etendue;
}

int main() {
	printf( "Solar photon flux: %.1f\n", N( 0,std::numeric_limits<double>::infinity(),6000,0, 6.8e-5) );
	printf( "Incoming power   : %.1f\n", P( 6000, 6.8e-5) );
	printf( "Amb   photon flux: %.1f\n", N( 0,std::numeric_limits<double>::infinity(), 300,0, M_PI - 6.8e-5) ); 	// absorbed ambient

	// simple balance without the upconverter 

	double incoming = P( 6000, 6.8e-5);	       	

	// both loops in eV units
	double step=0.05;
	for( double Eg=0.05; Eg<4.0; Eg+=step ) {
		printf( "---- Eg=%f\n", Eg);
		for( double u=0; u<Eg; u+=step ) {
	//		printf( "\n\t\t %f\n",N( Eg,std::numeric_limits<double>::infinity(),6000,0, 6.8e-5         ) );	// absorbed solar incoming 
	//		printf( "\t\t %f\n", N( Eg,std::numeric_limits<double>::infinity(), 300,0, M_PI - 6.8e-5 ) ); 	// absorbed ambient
	//		printf( "\t\t%f\n",-N( Eg,std::numeric_limits<double>::infinity(), 300,u, M_PI         ) ); 	// radiated back 

			double absorbed =  N( Eg,std::numeric_limits<double>::infinity(),6000,0, 6.8e-5        ) +	// absorbed solar incoming 
					   N( Eg,std::numeric_limits<double>::infinity(), 300,0, M_PI - 6.8e-5 ) + 	// absorbed ambient
					  -N( Eg,std::numeric_limits<double>::infinity(), 300,u, M_PI          ) ; 	// radiated back 

			// absorbed photon flux -> power=I*V -> u=qV
			printf( "%f %f %f\n" , Eg, u, (absorbed*toJoule(u)*100.0)/incoming );

		}
	}

	return 0;
}

