#include<stdio.h>
#include<math.h>

#ifndef Radial_h
#define Radial_h

double L_Z_A(double x);
double L_O_A(int alpha, double x);
double L_KP1_A(int k, int alpha, double Lk, double Lkm1, double x);
double Laguerre(int k, int alpha, double x);
double Radial_Orb(int n, int l, double r);
double factorial(int n);
double AngularIntegral(int l, int m, int lp, int mp)
#endif
