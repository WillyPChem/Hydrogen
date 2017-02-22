#include<stdlib.h>
#include<stdio.h>

int main()    {

int i, MAX_I; 
double x, dx,fx;

MAX_I = 4;
dx = 4./MAX_I;
for (i=0; i<+MAX_I; i++) {

   x = dx*i;
   fx = -1*(x-2)*(x-2) + 4;
   printf(" %i %f  %f   \n", i, x, fx);
}
} 
 
 
