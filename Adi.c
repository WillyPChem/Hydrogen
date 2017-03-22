
double AngularIntegral( int l, int m, int lp, int mp) {

double fl, flp;  


for (l=0; l<5; l++) {
for (lp=0; lp<5; lp++) {
for (m=0; m<5; m++) {
for (mp=0; mp<5; mp++) {
  
if (m==mp){


  if (l==lp+1){
    fl = sqrt(((l+1)*(l+1)-(m*m))/(4*(l+1)*(l+1)-1))
  }
  else { fl = 0}
if (l==lp-1) {
    flp = sqrt(((l*l-m*m)/(4*l*l-1))

  }
  else { flp = 0}
  
  }
 else { fl = 0;
        flp = 0 }

  }
  }
  }
  }


return integral;

}
