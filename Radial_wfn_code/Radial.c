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

void N_AngularIntegral(int li, int mi, int lf, int mf, double *mur, double *mui) {
  int i, j, max;

  double theta, phi;
  double dt, dp;

  double Yri, Yrf, Yii, Yif;
  double complex sum = 0. + 0*I;
  double pi=4.*atan(1.0);
  max = 5000;

  dt = pi/max;
  dp = 2*pi/max;

  if (mi==mf && ( (lf-li)==1 || (li-lf)==1) ) {
    for (i=0; i<=max; i++) {

      theta = i*dt;

      //for (j=0; j<=max; j++) {
        j=0;
        phi = j*dp;

        Spherical_Y(li, mi, theta, phi, &Yri, &Yii);
        Spherical_Y(lf, mf, theta, phi, &Yrf, &Yif);

        sum += (Yri - I*Yii)*(Yrf + I*Yif)*cos(theta)*sin(theta)*dt;

    //}
  }

  *mur = 2*pi*creal(sum);
  *mui = 2*pi*cimag(sum);

  }
  else {
  *mur = 0.;
  *mui = 0.;
  }

}

void Spherical_Y(int l, int m, double theta, double phi, double *YR, double *YI) {

  int mp;
  double ctheta, pfac, P;
  double complex y;

  mp = m;
  // Legendre Polynomial function will only take positive values of m
  if (m<0) {

    mp = abs(m);
  }

  // Prefactor for Y
  pfac = prefac( mp, l );

  // cosine of theta
  ctheta = cos(theta);

  // Legendre Polynomial P_l^m (cos(theta))
  P = plgndr( l, mp, ctheta);

  // Spherical Harmonic = prefac*P_l^m(cos(theta))*exp(i*m*phi)
  y = pfac*P*cexp(I*m*phi);

  *YR = creal(y);
  *YI = cimag(y);

}

double plgndr(int l, int m, double x) {
//Computes the associated Legendre polynomial P m
//l (x). Here m and l are integers satisfying
//0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1.
//void nrerror(char error_text[]);

  float fact,pll,pmm,pmmp1,somx2;

  int i,ll;

  if (m<0) m*=-1;

  if (m < 0 || m > l || fabs(x) > 1.0) {
  //if (m>l || fabs(x) > 1.0) {
    printf("Bad arguments in routine plgndr\n");
    printf("  l is %i and m is %i\n",l,m);
    exit(0);
  }

  pmm=1.0;   //Compute P^m_m .

  if (m > 0) {

    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;

    for (i=1;i<=m;i++) {

      pmm *= -fact*somx2;
      fact += 2.0;

    }
  }

  if (l == m)
    return pmm;

  else {    //Compute P^m_m+1

    pmmp1=x*(2*m+1)*pmm;

    if (l == (m+1))

      return pmmp1;

    else {   //Compute P^m_l, l>m+1

      for (ll=m+2;ll<=l;ll++) {

        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;

     }

     return pll;
   }
  }
}


double prefac(int m, int l) {


  double p, num1, num2, denom1, denom2;
  double pi=4.*atan(1.0);
  num1 = 2*l+1;
  num2 = factorial( (l-m) );

  denom1 = 4*pi;
  denom2 = factorial( (l+m) );


  p = sqrt((num1/denom1)*(num2/denom2));

  return p;

}

double Legendre(int l, int m, double theta) {


  int mp;
  double ctheta, pfac, P;
  double y;

  mp = m;
  // Legendre Polynomial function will only take positive values of m
  if (m<0) {

     mp = abs(m);
  }
  // Prefactor 
  pfac = prefac( mp, l );

  // cosine of theta
  ctheta = cos(theta);

  // Legendre Polynomial P_l^m (cos(theta))
  P = plgndr( l, mp, ctheta);

  // Spherical Harmonic = prefac*P_l^m(cos(theta))*exp(i*m*phi)
  y = pfac*P;

  return y;
}


