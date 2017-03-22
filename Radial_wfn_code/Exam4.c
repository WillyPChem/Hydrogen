#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<Radial.h>

int main()    {

int i, MAX_r; 
double r, j;
double frdr, dr, fr, sum;
double e;

e = 1;
MAX_r = 1000;
dr = 100./MAX_r;


sum = 0.;
for (i=0; i<=MAX_r; i++) {
 r = dr*i;
 // Set R10 equal to Radial Function R_1,0 - 1 s orbital
 double R10 = Radial_Orb(n,l,r);
 // set R20 equal to Radial Function R_2,1 - 2 p orbital
 double R20 = Radial_Orb(np,lp,r);

 fr = R10 * e * r * r* r * R20;
 sum = sum + fr * dr;

}

printf("  Sum is %f\n",sum);
double dipole_integral = sum*AngularIntegral(l, m, lp, mp);
} 

 

double AngularIntegral( int l, int m, int lp, int mp) {

double fl, flp;  



 fl = 0;
 flp = 0;
  
 if (m==mp){


   if (l==lp+1){
     fl = sqrt(((l+1)*(l+1)-(m*m))/(4*(l+1)*(l+1)-1));
   }
   else { fl = 0;}

   
   if (l==lp-1) {
     flp = sqrt(((l*l-m*m)/(4*l*l-1));

    }
    else { flp = 0;}
  
  }
 else { fl = 0;
        flp = 0; }


integral = fl*flp;
                
printf(" Angular Integral is %f\n",integral);
return integral;

}
                
                
