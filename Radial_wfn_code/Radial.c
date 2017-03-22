#include"Radial.h"

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
  /*else if (k==2) {
     Lk = L_O_A(alpha, x);
     Lkm1 = L_Z_A(x);
     L_KP1_A( k, alpha, Lk, Lkm1, x);
  }
  */
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

//R_n,l(r) = sqrt( (2/(n*a0))^3  * (n - l - 1)!/(2*n((n+l)!) ) )  * exp(-r/(n*a0)) * (2r/(n*a0))^l  L_{n-l-1}^{2l+1} (2*r/(n*a0))
double Radial_Orb(int n, int l, double r) {
  double a0 = 1.;
  // Prefactor 1:  sqrt( (2/(n*a0))^3  * (n - l - 1)!/(2*n((n+l)!) ) )
  double term1 = pow((2./(n*a0)),3.);
  double term2a = factorial ( (n - l - 1) );
  double term2b = 2*n*factorial(n+l);
  double term2aob = term2a/term2b;

  double pre1 = sqrt(term1*term2aob);

  // Term 2:  exp(-r/(n*a0))
  double term2;
  term2 = exp(-r/(n*a0));

  // Term 3: (2r/(n*a0))^l 
  double term3;
  term3 = pow((2*r/n),l);
  //printf("  term1 is %f\n",term1);

  // Term 4:  L_{n+1}^{2l+1} (2*r/(n*a0)) 
  double term4;
  term4 = Laguerre( (n-l-1), (2*l+1), (2*r/(n*a0)) );
  //printf("  term3 is %f\n",term3);

  return pre1*term2*term3*term4;

}

double factorial(int n) {
  int i;
  double result = 1;
 
  for (i = 1; i <= n; i++) {
    result = result * i;
  }
  return result;
}

double AngularIntegral( int l, int m, int lp, int mp) {

double fl, flp, integral;  



 //fl = 0;
 //flp = 0;
  
 if (m==mp){


   if (l==(lp+1)){
     
     fl = sqrt(((l+1.)*(l+1.)-(m*m))/(4.*(l+1.)*(l+1.)-1.));
   }
   else { fl = 0;}

   
   if (l==(lp-1)) {
     flp = sqrt(((l*l-m*m)/(4.*l*l-1.)));
    }
    else { flp = 0;}
  
  }
 else { fl = 0;
        flp = 0; }


integral = fl+flp;
                
printf(" Angular Integral is %f\n",integral);
return integral;

}
