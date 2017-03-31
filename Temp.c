#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Radial_wfn_code/Radial.h"
int main() {
 int dim;
 int i, n, l, q, m, k, MAX_r;
 int np, lp, mp;
 int j, Norbs;
 double sum, r, dr, fr, e;
 
 e = 1.0;
 dim = 4;

 // Array B is to hold information about the quantum numbers n, l, and m
 // associated with a given orbital labelled with a single number.
 // E.g. orbital 1 is the 1s orbital, it has n=1, l=0, m=0
 int B[dim][3];
 double E[dim];
 // Consider n values from n=1 to n=4
 int orb_number = 0;
 for (n=1; n<=2; n++) {

   for (l=0; l<n; l++) {

     for (m=-l; m<=l; m++) {


       B[orb_number][0] = n;
       B[orb_number][1] = l;
       B[orb_number][2] = m;
       E[orb_number] = -13.6/(n*n);

       printf(" orbital number %i  n=%i, l=%i, m=%i\n",orb_number,n,l,m);
       printf(" Energy of orbital %i is %f\n",orb_number,E[orb_number]);
       orb_number++;
       


     }
   }
 }
 
double mur, mui;
double H[dim][dim];
for (i=0; i<dim; i++) {
 
  
  H[i][i] = E[i];

  for (j=0; j<dim; j++) {
    if (i!=j) {

       n = B[i][0];
       l = B[i][1];
       m = B[i][2]; 

       np = B[j][0];
       lp = B[j][1];
       mp = B[j][2]; 
     
     
     MAX_r = 1000;
     dr =  200/MAX_r;
     for (k=0; k<=MAX_r; k++) {
       r = dr*k;
       // Set R10 equal to Radial Function R_1,0 - 1 s orbital
       double R10 = Radial_Orb(n,l,r);
       // set R20 equal to Radial Function R_2,1 - 2 p orbital
       double R20 = Radial_Orb(np,lp,r);

       fr = R10 * e * r * r* r * R20;
       sum = sum + fr * dr;

     }

     N_AngularIntegral(l, m, lp, mp, &mur, &mui);
     H[i][j] = sum*mur;

     // get n, l, and m associated with orbital i
     // get n', l', and m' associated with orbital j [Use B vector for this]
     // duplicate code which calculated the radial integrals
     // call angular integral function
     // Multiply radial integral value by angular integral value
     //  here you need to calculate the radial integral
               //  and the angular integral for orbital i against orbital j
               //  H[i][j] = 0.;

    }
  }
}

// Making a small change
for (i=0; i<dim; i++) {

  for (j=0; j<dim; j++) {

    printf(" %f  ",H[i][j]);

  }
  printf("\n");
} 
  
}

/*int [i]= ( n,l,m )
int [j]= ( np, lp, mp )
  B[i][0] = n;
  B[i][1] = l;
  B[i][2] = m; 

  B[j][0] = np;
  B[j][1] = lp;
  B[j][2] = mp; 


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
*/


/*for (i=0; i<30; i++) {
for (n=1; n
 Norbs(n)=n*n;
 i=Norbs+1
 B[Norbs][0]=1
 E=
  H[i][i] = -13.6/(pow(n , 2));

  for (j=0; j<5; j++) {
    if (i!=j) H[i][j] = 0.;
  }
}

for (i=0; i<5; i++) {

  for (j=0; j<5; j++) {

    printf("  %f  ",H[i][j]);

  }
  printf("\n");
} 
  
}*/
