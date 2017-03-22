
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
