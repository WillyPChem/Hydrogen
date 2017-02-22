#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main() {
 
 int i, n, l, q, m;
 int j, Norbs;
 // Array B is to hold information about the quantum numbers n, l, and m
 // associated with a given orbital labelled with a single number.
 // E.g. orbital 1 is the 1s orbital, it has n=1, l=0, m=0
 int B[30][3];
 double E[30];
 // Consider n values from n=1 to n=4
 int orb_number = 0;
 for (n=1; n<=4; n++) {

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
} 

double H[30][30];

for (i=0; i<30; i++) {
 
  
  H[i][i] = E[i];

  for (j=0; j<30; j++) {
    if (i!=j) H[i][j] = 0.;
  }
}

for (i=0; i<30; i++) {

  for (j=0; j<30; j++) {

    printf(" %i, %j, %n  %f  ",H[i][j]);

  }
  printf("\n");
} 
  
}




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
