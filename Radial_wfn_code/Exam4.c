#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<Radial.h>

const int dim = 55;
int MAX_TIME;
double H[dim][dim];
double B[100][3];
double E[dim];
double complex C[dim];
double complex Cdot[dim];
double complex Ct[dim];
double complex D[dim][dim];
double mu[dim][dim];
double dt = 0.001;
// sigma is related to the electric field
// specifically, when the simulation time is equal to sigma, then
// the laser amplitude will be zero
// 804 atomic units of time is about 19 femto seconds
double sigma = 80.0;

// Good definition of pi
double pi = 4.*atan(1.0);
// Electric field amplitude in atomic units
double Emax=0.04;

void HdotC(double t);
void RK3(double t);
double EField(double t);
double DensityMatrix();

int main()    {

 int i, j, k;
 int n, l, m, np, lp, mp;
 double sum, r, dr, fr, e, R10, R20;
 double mur, mui;
 double *dpr, *dpi;
 double *spectrum_real, *spectrum_imag;
 int orb_number;

 
 // Maximum number of wavefunction updates is specified by MAX_TIME..
 // the actual amount of time in the simulation will be equal to MAX_TIME*dt (time in atomic units)
 // The simulation should proceed for at least enough time for the electrons to make their slowest orbital transition
 // which corresponds to the smallest energy gap in the system... this will be the gap between orbitals
 // with principle quantum number n=4 and n=5 -> delta E = 0.5 atomic units * (1/4^2 - 1/5^2) = 0.03125 atomic units
 // The frequency of this transition is delta E/hbar = 0.03125 inverse atomic units of time
 // The period is then 1/frequency = 1/0.03125 = 32 atomic units of time
 // Therefore, MAX_TIME*dt >= 32 -> MAX_TIME >= 32/dt... for good measure, we will sample ten periods
 double Min_deltaE = 0.5*(1./(4*4) - 1./(5*5));
 double LongPeriod = 1./Min_deltaE;
 MAX_TIME = 10*(int)(LongPeriod/dt); 
 e = 1.;
 int MAX_r = 1000;
 dr = 200./MAX_r;

 // Dipole moment vector needs to be MAX_TIME long
 dpr = (double *)malloc(MAX_TIME*sizeof(double));
 dpi = (double *)malloc(MAX_TIME*sizeof(double));
 spectrum_real = (double *)malloc(MAX_TIME*sizeof(double));
 spectrum_imag = (double *)malloc(MAX_TIME*sizeof(double));
 orb_number=0;
 // Initialize wavefunction vector as ground energy eigenstate
 C[orb_number] = 1. + 0.*I;
 // complex conjugate of initial wavefunction vector
 Ct[orb_number] = conj(C[orb_number]);

 // Build arrays that map from orbital number to quantum numnbers n, l, m
 // get orbital energy while we are at it
 for (n=1; n<=5; n++) {

   for (l=0; l<n; l++) {

     for (m=-l; m<=l; m++) {


       B[orb_number][0] = n;
       B[orb_number][1] = l;
       B[orb_number][2] = m;
       E[orb_number] = -0.5/(n*n);  // Energy in atomic units!

       orb_number++;



     }
   }
 }

 // Build Hamiltonian matrix and Dipole matrix
 for (i=0; i<dim; i++) {

   // Hamiltonian is diagonal - just orbital energies
   H[i][i] = E[i]-E[0];
  
   for (j=0; j<dim; j++) {
      if (i!=j) {
 
        n = B[i][0];
        l = B[i][1];
        m = B[i][2];

        np = B[j][0];
        lp = B[j][1];
        mp = B[j][2];
        //printf("  Going to compute (%i %i | %i %i)\n",n,l,np,lp);
        sum = 0.;
        for (k=0; k<=MAX_r; k++) {
       
          r = dr*k;
          // Set R10 equal to Radial Function R_1,0 - 1 s orbital
          R10 = Radial_Orb(n,l,r);
          // set R20 equal to Radial Function R_2,1 - 2 p orbital
          R20 = Radial_Orb(np,lp,r);

          fr = R10 * e * r * r* r * R20;
          sum = sum + fr * dr;
        }
        //printf("  (%i %i | %i %i) : %12.10e\n",n,l,np,lp,sum);
        N_AngularIntegral(l, m, lp, mp, &mur, &mui);
        //printf("  (%i %i %i | mu | %i %i %i): %12.10e\n",n,l,m,np,lp,mp,sum*mur);
        // Dipole is off-diagonal - elements are transition dipole integrals
        mu[i][j] = sum*mur;

      }
    }
  }

  //HdotC(0);

  // Loop over time;
  // Step 1:  Modify function EField to correspond to a laser pulse (as in 
  //          the J. Phys. Chem article)
  // Step 2:  Store each dipole moment in a vector with the same number
  //          of elements as there are timesteps... currently, there are
  //          50000 timesteps, so your vector must be 50000 elements long
  //          Note that the Fourier transform expects a real and imaginary 
  //          array of values, because your dipole moment is real, you will
  //          probably want a parner array of imaginary values (which will all 
  //          be zero) that must be 50000 elements long, as well
  // Step 3:  Implement the Fourier transform routine and use it on your
  //          array of dipole values.  The result can be used to compute 
  //          the absorption spectrum

  //          FourierTransform(double *dipole_real, double *dipole_imag, double *FT_real, double *FT_imag);
  //          And absorption spectrum = FT_real*FT_real + FT_imag*FT_imag
  double Tim, rho, dpm, ef;
  FILE *fp;
  fp = fopen("dipoleMoment_2.txt","w");
  for (int M=0; M<MAX_TIME; M++) {

    Tim=M*dt;
    RK3(Tim);
  
    rho = 0.;
    // Real part of dipole moment
    dpr[M] = DensityMatrix();
    // Dipole moment is real, imaginary part is zero
    dpi[M] = 0.;
    ef = EField(Tim);
    fprintf(fp,"%12.10e  %12.10e %12.10e\n",M*dt,dpr[M],ef);

  }

  fclose(fp);
/* 
 * Discrete Fourier transform
 * by Project Nayuki, 2017. Public domain.
 * https://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
 */
 
 compute_dft(dpr, dpi, spectrum_real, spectrum_imag, MAX_TIME);
 FILE *absfp;
 absfp = fopen("AbsorptionSpectrum.txt","w");
 
  double OMEGA_min = 0.050;
  double OMEGA_max = 1.;
  double dOMEGA = (OMEGA_max - OMEGA_min)/MAX_TIME;
 for (int i = 0; i < MAX_TIME; i++)
    double OMEGA = dOMEGA * i + OMEGA_min;
    double abs = spectrum_real[i]*spectrum_real[i] + spectrum_imag[i]*spectrum_imag[i];
    fprintf(absfp," %12.10e  %12.10e\n",OMEGA, abs);
}
fclose(absfp);
 // after you call compute_dft - print frequency, spectrum_real^2 + spectrum_imag^2
 //void compute_dft(const double inreal[], const double inimag[], double outreal[], double outimag[], int n



 
}
  
/* 
 * Computes the discrete Fourier transform (DFT) of the given vector.
 * All the array arguments must have the same length.
 */
void compute_dft(const double *inreal, const double *inimag, double *outreal, double *outimag, int n) {
	for (int k = 0; k < n; k++) {  /* For each output element */
		double sumreal = 0;
		double sumimag = 0;
  double OMEGA_min = 0.050;
  double OMEGA_max = 1.;
      double dOMEGA = (OMEGA_max - OMEGA_min)/MAX_TIME;
      double OMEGA = dOMEGA * k + OMEGA_min;

  double absor
		for (int t = 0; t < n; t++) {  /* For each input element */
			double angle = 2 * M_PI * t * OMEGA / n;
			sumreal +=  inreal[t] * cos(angle) + inimag[t] * sin(angle);
			sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
		}
		outreal[k] = sumreal;
		outimag[k] = sumimag;
	}
}
  // NOTE!!!  HERE YOU NEED TO CALL YOUR FOURIER TRANSFORM FUNCTION
  // AND YOU WILL GIVE IT THE dpr and dpi ARRAYS THAT WERE JUST COMPUTED ABOVE!

    /*for (i=0; i<dim; i++) {
 
      //printf("  %12.10e  %12.10e\n",creal(C[i]),cimag(C[i]));
      rho += creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
    }
    printf("  TRACE IS %12.10e\n",rho);
    */
  /*  Uncomment if you want to print the Hamiltonian/dipole matrix
  for (i=0; i<dim; i++) {

    for (j=0; j<dim; j++) {

      printf(" %f  ",H[i][j]);

    }
    printf("\n");
 }
 */

}

// i Cdot = (H - E(t)*mu) * C
void HdotC(double t) {
  double complex sum;  

  int i, j;
  for (i=0; i<dim; i++) {
    
    // Matrix-vector product of (H-E(t)*mu) * C
    sum = 0. + 0.*I;
    for (j=0; j<dim; j++) {

      sum += (H[i][j] - EField(t)*mu[i][j])*C[j];

    }
    // Cdot[i] = -I*sum
    Cdot[i] = -I*sum;
    //printf("  (%12.10e, %12.10e) -> (%12.10e, %12.10e)\n",creal(C[i]),cimag(C[i]),creal(Cdot[i]),cimag(Cdot[i]));
  }
}

//  Function that defines that laser field experienced
//  by the molecule at time point t
double EField(double t) {
  double freq, amplitude;

  // Envelope function for pulse - can change
  amplitude = Emax*sin(pi*t/(2*sigma))*sin(pi*t/(2*sigma)); 

 // Frequencies in pulse - can change these as well
 double  freq1 = -0.5/4. + 0.5/1.;
 double  freq2 = -0.5/9. + 0.5/4.;
 double  freq3 = -0.5/16. + 0.5/9.;
 double ef = amplitude*(sin(freq1*t) + sin(freq2*t) + sin(freq3*t));
 
return  ef;

  
}
// Updates wavefunction - aka solves TDSE 
void RK3(double t) {
  int i;
  // Need several arrays!
  double complex *k1, *k2, *k3, *c1;
  
  // allocate memory for them!
  k1 = (double complex *)malloc(dim*sizeof(double complex));
  k2 = (double complex *)malloc(dim*sizeof(double complex));
  k3 = (double complex *)malloc(dim*sizeof(double complex));
  c1 = (double complex *)malloc(dim*sizeof(double complex));

  // Get Cdot at time t for C
  HdotC(t);

  for (i=0; i<dim; i++) {

    // Euler step
    k1[i] = dt*Cdot[i];
    // Store C vec at time t
    c1[i] = C[i];
    // Partial update to C vec
    C[i] = c1[i] + k1[i]/2.;
  }  

  // Get Cdot at time t+dt/2 for C+k1/2
  HdotC(t+dt);

  for (i=0; i<dim; i++) {

    // Euler step
    k2[i] = dt*Cdot[i];
    // Partial update to C vec
    C[i] = c1[i] + k2[i]/2.;
  }
  
  // Get Cdot at time t+dt/2 for C+k2/2
  HdotC(t+dt);

  for (i=0; i<dim; i++) {

    // Euler step
    k3[i] = dt*Cdot[i];

    C[i] = c1[i] + k1[i]/6. + 2*k2[i]/3. + k3[i]/6;

  }
   
  

  free(k1);
  free(k2);
  free(k3);
  free(c1);


}

// Compute Density Matrix, Trace, and Trace(D Mu)
double DensityMatrix() {

  double Trace;
  double complex dipole;
  Trace = 0.;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {

      D[i][j] = C[i]*conj(C[j]); 
      dipole += D[i][j]*mu[i][j];
   
    }
    Trace += creal(D[i][i]);
  }

  printf("  TRACE DM IS %12.10e\n",Trace);

  return creal(dipole);

}
