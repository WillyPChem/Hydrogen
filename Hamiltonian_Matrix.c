#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main() {

int i, j,  n;  
double H[5][5];

for (i=0; i<5; i++) {
 
  n = i + 1;
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
  
}

