#include<stdlib.h>
#include<stdio.h>
#include<math.h>


int main()    {

int i, MAX_I; 
double x, fx, dx, sum;
double  fx2 = fx*fx;
double pi;


pi = 4.*atan(1.0);

MAX_I = 100000000;
dx = 10./MAX_I;



sum = 0.;
for (i=0; i<=MAX_I; i++) {

   x = dx*i;
   
   fx = sqrt(2./10.)*sin(pi*x/10.);
   sum = sum +fx*dx;

  printf(" %i %f %f %f\n", i, x, fx, sum );
}
printf( " integral is equal to %f\n",sum);
} 
 
 
