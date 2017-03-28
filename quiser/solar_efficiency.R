options(digits=22)

h <- 6.62607004e-34
hbar <- 1.054571800e-34 
cvac <- 299792458
c <- cvac
k    <- 1.38064852e-23
evj  <- 1.60218e-19
es   <- 6.8e-5
u    <- 0

normalizationFactor <- (4*pi^3*hbar^3*cvac^2)^-1

N <- function(E1,E2,T,u,e) {
  planck <- function(E) {
    (E*evj)^2 * ( exp(((E-u)*evj)/(k*T))-1 )^-1
  }
  #e*normalizationFactor*integrate(planck,E1,E2)$val
}

N(1,10,6000,0,es)
N(1,100,6000,0,es)
N(1,1000,6000,0,es)
N(1,10000,6000,0,es)
N(1,Inf,6000,0,es)




AN <- function(E,T,u) {
  c1 = (2*k^3*T^3)/(c^2*h^3);
  c2 = (4*u*k^2*T^2)/(c^2*h^3);
  c3 = (2*k*T*u^2)/(c^2*h^3);
  
  res=0.0;
  
  x=(E*evj - u*evj)/(k*T);
  
  for( n in 1:10000 ) {
    res = res + c1*(2/n^3+2*x/n^2+E^2/n)*exp(-n*x);
    res = res + c2*(n^-2+x/n)*exp(-n*x);
    res = res + c3*exp(-n*x);
  }
  res;
}

