#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

#ifndef Radial_h
#define Radial_h

double L_Z_A(double x);
double L_O_A(int alpha, double x);
double L_KP1_A(int k, int alpha, double Lk, double Lkm1, double x);
double Laguerre(int k, int alpha, double x);
double Radial_Orb(int n, int l, double r);
double factorial(int n);
void Spherical_Y(int l, int m, double theta, double phi, double *YR, double *YI);
double Legendre(int l, int m, double theta);
double prefac(int m, int l);
double factorial(int n);
double plgndr(int l, int m, double x);
void N_AngularIntegral(int li, int lf, int mi, int mf, double *mur, double *mui);

#endif
