#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

double pi=4.*atan(1.0);
void Spherical_Y(int l, int m, double theta, double phi, double *YR, double *YI);
double Legendre(int l, int m, double theta);
double prefac(int m, int l);
double factorial(int n);
double plgndr(int l, int m, double x);
void N_AngularIntegral(int li, int lf, int mi, int mf, double *mur, double *mui);
//void TransitionDipole_Simp(int ni, int nf, int li, int lf, int mi, int mf, double R, double *mur, double *mui);

double Ang(int l, int m, int lp, int mp);

int main() {

  double int_r, int_i;
  int li, lf, mi, mf;

  for (li=1; li<=3; li++) {

    for (mi=-li; mi<=li; mi++) {

      lf=li-1;
      
      for (mf=-lf; mf<=lf; mf++) { 

        N_AngularIntegral(li,mi,lf,mf,&int_r, &int_i);

        printf("  (%i %i| %i %i) %12.10f + i*%12.10f\n",
            li,mi,lf,mf,int_r, int_i);
      }
      
      lf=li+1;

      for (mf=-lf; mf<=lf; mf++){
        N_AngularIntegral(li,mi,lf,mf,&int_r, &int_i);

        printf("  (%i %i| %i %i) %12.10f + i*%12.10f\n",
          li,mi,lf,mf,int_r, int_i);
      }

    }
  }

}

double Ang(int l, int m, int lp, int mp) {

  double num1, num2, denom1, denom2, t1, t2;
  double intval;
  num1 = (l+1.)*(l+1.)-m*m;
  denom1 = 4.*(l+1.)*(l+1.) - 1.;
  t1 = sqrt(num1/denom1);

  num2 = l*l-m*m;
  denom2 = 4.*l*l-1.;
  t2 = sqrt(num2/denom2);

  t2=0;
  if (m!=mp) { 

    intval = 0.;
    

  }
  else {

    if (l==(lp+1)) {

      intval = t1;      

    }
    else if (l==(lp-1)) {

      intval = t2;
    }
    else {

      intval = 0.;
    }

  }

  return intval;

}

void N_AngularIntegral(int li, int mi, int lf, int mf, double *mur, double *mui) {
  int i, j, max;

  double theta, phi;
  double dt, dp;

  double Yri, Yrf, Yii, Yif;
  double complex sum = 0. + 0*I;

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

double factorial(int n) {
  int c;
  double result = 1;

  for (c = 1; c <= n; c++)
    result = result * c;

  return result;
}

double prefac(int m, int l) {


  double p, num1, num2, denom1, denom2;

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

