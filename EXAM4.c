#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<Radial.h>

int main()    {

int i, MAX_r; 
double r, j;
double frdr, dr, fr, sum;
double e;

e = 1


MAX_r = 100;
dr = 10./MAX_r;


 sum = 0.;
for (i=0; i<=MAX_r; i++) {
 r = dr*i;
 double R10 = Radial_Orb(1,0,r);
 double R20 = Radial_Orb(2,0,r);

   
 
   
   fr = R10 * e * r * r* r * R20;
   sum = sum + fr * dr

  printf(" %i %f  %f %f  \n", r, R, fr, );
}
} 
 
 
