#include<stdlib.h>
#include<stdio.h>
#include<math.h>


int main()    {

int i, MAX_I; 
double x, dx, fx, fx2;
double pi;


pi = 4.*atan(1.0);

MAX_I = 6;
dx = 10./MAX_I;




for (i=0; i<=MAX_I; i++) {

   x = dx*i;
   fx = sqrt(2./10.)*sin(pi*x/10.);
   fx2 = fx*fx;
  printf(" %i %f  %f   \n", i, x, fx2);
}
} 
 
 
