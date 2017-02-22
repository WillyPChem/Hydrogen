#include<stdlib.h>
#include<stdio.h>
#include<math.h>


int main()    {

int i, MAX_I; 
double x, fxdx, fxpdx, dx, fx, xpdx, fp, m;
double  fx2 = fx*fx;
double pi;


pi = 4.*atan(1.0);

MAX_I = 1000;
dx = 10./MAX_I;




for (i=0; i<=MAX_I; i++) {

   x = dx*i;
   xpdx = x + dx;
   
   fx = sqrt(2./10.)*sin(pi*x/10.);
   fp = sqrt(2./10.)*(pi/10.)*cos(pi*x/10.);
   fxpdx = sqrt(2./10)*sin(pi*xpdx/10.);
   m = (fxpdx-fx)/(xpdx-x);

  printf(" %i %f  %f %f  \n", i, x, fp, m);
}
} 
 
 
