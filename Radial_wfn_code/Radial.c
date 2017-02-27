#include<stdio.h>
#include<math.h>


// Bohr radius in SI units
double a0_SI = 5.2917721067e-11;
// Bohr radius in atomic units
double a0 = 1.0;

double L_Z_A(double x);
double L_O_A(int alpha, double x);
double L_KP1_A(int k, int alpha, double Lk, double Lkm1, double x);
double Laguerre(int k, int alpha, double x);
double Radial_Orb(int n, int l, double r);
double factorial(int n);

int main() {

  double pi = 3.14159265359;
  double Z = 1;  // Nuclear charge of hydrogen in atomic units
  double a0 = 1;  // Bohr radius in atomic units
  double h = 1;  // Planck's constant in atomic units
  double m = 1;  // mass of electron in atomic units
  int i, n, l;

  double value, x, dvalue;
  for (i=0; i<100; i++) {
    x = -2 + 0.1*i; 
    value = Laguerre( 3, 2, x);
    printf("  %f  %f \n",x,value);
  }

  return 0;
 }

double L_Z_A(double x) {

  return 1.;

}

double L_O_A(int alpha, double x) {
  double term;
  term = 1. + alpha - x;
  return term;

}

double L_KP1_A(int kp1, int alpha, double Lk, double Lkm1, double x) {
  int k = kp1-1;
  double num, denom;
  num = (2*k+1+alpha-x)*Lk-(k+alpha)*Lkm1;
  denom = k+1;
  return num/denom;

}

// k is the subscript and alpha is the superscript
// in the associated Laguere polynomial
// L_k^alpha ( x ) 
double Laguerre(int k, int alpha, double x) {

  double Lkp1, Lk, Lkm1;
  int tk;

  if (k==0) return L_Z_A( x );
  else if (k==1) return L_O_A(alpha, x);
  else if (k==2) {
     Lk = L_O_A(alpha, x);
     Lkm1 = L_Z_A(x);
     L_KP1_A( k, alpha, Lk, Lkm1, x);
  }

  else {

    tk = 1;
    Lk = L_O_A(alpha, x);
    Lkm1 = L_Z_A( x );

    while (tk<k) {
 
      Lkp1 = L_KP1_A( tk+1, alpha, Lk, Lkm1, x);
      tk++;
      Lkm1 = Lk;
      Lk = Lkp1;

    }
    return Lkp1;  

  }
}

//R_n,l(r) = sqrt( (n-1-l)!/(2*n*((n+l)!)^3) ) * (2/(n*a0))^{l+3/2} * r^l * exp(-r/(n*a0)) * L_{n+1}^{2l+1} (2*r/(n*a0))
double Radial_Orb(int n, int l, double r) {

  // Prefactor 1:  sqrt( (n-1-l)!/(2*n*((n+l)!)^3) )
  double pre1, num, numfac, denom, denomfac;
  num = (n - 1 - l);
  numfac = factorial(num);
  denomfac = factorial((n+1));
  denom = 2*n*pow(denomfac,3);
  pre1 = sqrt(numfac/denom);

  // Prefactor 2: (2/(n*a0))^{l+3/2}
  double pre2, arg2, pow2;
  pow2 = l+3./2.;
  arg2 = (2./(n*a0));
  pre2 = pow(arg2,pow2);

  // Term 1: r^l
  double term1;
  term1 = pow(r,l);

  // Term 2:  exp(-r/(n*a0))
  double term2;
  term2 = exp(-r/(n*a0));

  // Term 3:  L_{n+1}^{2l+1} (2*r/(n*a0)) 
  double term3;
  term3 = Laguerre( (n+1), (2*l+1), (2*r/(n*a0)) );

  return pre1*pre2*term1*term2*term3;

}

double factorial(int n) {
  int i;
  double result = 1;
 
  for (i = 1; i <= n; i++) {
    result = result * i;
  }
  return result;
}
