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
dr = 20./MAX_r;


sum = 0.;
for (i=0; i<=MAX_r; i++) {
 r = dr*i;
 double R10 = Radial_Orb(1,0,r);
 double R20 = Radial_Orb(2,0,r);

 fr = R10 * e * r * r* r * R20;
 sum = sum + fr * dr;

}

printf("  Sum is %f\n",sum);

} 
int main {
 
double AngularIntegral( int l, int m, int lp, int mp) {


kroneckerDelta(m,mp)*(sqrt(((l+1)*(l+1)-(m*m))/(4*(l+1)*(l+1)-1))*kroneckerDelta(l,lp+1) + sqrt(((l*l-m*m)/(4*l*l-1))*kroneckerDelta(l,lp-1))

return integral;
                      }
}
