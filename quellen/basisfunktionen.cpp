/*
 * basisfunctions.cpp
 *
 * Created on: Jun 24, 2013
 * Author: cfdickmann
 */

#include "AmericanOption.h"

using namespace std;

double dnd(double x) {
	return 1. / sqrt(2. * 3.141593) * exp(-x * x / 2.);
}

double cnd(double x) {
	return 0.5 * (1. + erf(x / sqrt(2)));
}

double trig(double x, int k) {
	if (k % 2 == 0)
		return cos((double) (k / 2) * x);
	else
		return sin((double) ((k + 1) / 2) * x);
}

double eingepasst(double x, double f) {
	double g = 100;
	double asymp = x < g ? 0 : 1;
	double gewicht = 1. / (1. + 0.1* pow(x - g, 2));
	return gewicht * f + (1. - gewicht) * asymp;
}

double AmericanOption::f1D(double x, int k, int d) {
return eingepasst(x,trig(x, k))   *x;
	}

double Sum(double* x, int d){
	double s=0;
	for(int i=0;i<d;++i)
		s+=x[i];
	return s;
}

double AmericanOption::f2D(double* x, int k, int d) {
	if (k < KpI / 4) { //Funktionen, die vom Maximum abhaengen
		if (x[d] == Max(x, D))
			return trig(x[d], k) * x[d];
		else
			return 0;
	}

	if (k < 2 * KpI / 4) { //Funktionen, die von den einzelnen assets abhaengen
		k -= 1 * KpI / 4;
		if ((k % KpI) % D == d)
			return trig(x[d], k / D) * x[d];
		else
			return 0;
	}

	if (k < 3 * KpI / 4) { //Funktionen, die von der Summe abhaengen
		k -= 2 * KpI / 4;
		double S = Sum(x, D) / (double) (D);
		return trig(S, k) * x[d];
	}

	if (k < 4 * KpI / 4) { //Funktionen, die von der Differenz abhaengen
		k -= 3 * KpI / 4;
		double faktor = -1;
		double max = Max(x, D);
		if (x[d] == max)
			faktor = 1;
		double S = max - (Sum(x, D) - max) / (double) (D - 1);
		return trig(S, k) * faktor * x[d];
	}

	printf("Error 477: etwas mit der Abzaehlung der k stimmt nicht.");
	return -1;
}

double AmericanOption::trueE(double x, int n){
	int ex=0;
	for(int e=0;e<number_of_Exercise_Dates;++e)
	if(Exercise_Dates[e]<=n)ex=e;

	double xx[1];
	xx[0]=x;
	double C=LSM_C_estimated(xx,ex);

	double P=0*(x-Strike)>0?(x-Strike):0;
		return C>P?C:P;
}

double pos(double x){
	return x>0?x:0;
}

//double line(double x, int k, int K){
//
//
//	return pos(x-);
//
//}

//double step(double x, double stair){
//if(x>)
//
//}

double AmericanOption::f(int k, double **x, int n, int d) {
// if (k == 2)
// return 0;
//	return 1.;
	//return 1;

	if (k / KpI !=	 n / NpI)
		return 0;

//	double y[D];
//	for (int dd = 0; dd < D; ++dd)
//		y[dd] = x[n][dd] / (2. * Strike);       //Achtung werden nicht immer benutzt

	int kAC = k % KpI;

	if (D == 1){
//		if(kAC==0) return pos(x[n][0]-Strike);
//		return trig(0.005*x[n][0],kAC);
if(kAC==0)		return EB.call_diff((double)n*dt,T,x[n][0],Strike,r,delta,sigma[0])*x[n][0];

//		return (x[n][0]>Strike*3.*(double)(kAC)/(double)(KpI))*x[n][0] ;


//return (trueE(x[n][0]+0.000001,n)-trueE(x[n][0],n))/0.000001*x[n][0];
//		if(kAC==0)return EB.call_diff(n*dt,1.0*T,x[n][0],Strike,r,delta,sigma[0])*x[n][0];
//		if(kAC==1)return EB.call_diff(n*dt,1.2*T,x[n][0],Strike,r,delta,0.7*sigma[0])*x[n][0];
//		if(kAC==2)return EB.call_diff(n*dt,1.8*T,x[n][0],Strike,r,delta,1.2*sigma[0])*x[n][0];

//		return (x[n][0]>Strike?1:0);
//		if(k==0) return EB.call_diff((double)dt*n,1.0*T,x[n][0],Strike,r,delta,sigma[0])*x[n][0];
//		if(k==1) return EB.call_diff((double)dt*n,1.2*T,x[n][0],Strike,r,delta,0.8*sigma[0])*x[n][0];
//		if(k==2) return EB.call_diff((double)dt*n,1.5*T,x[n][0],Strike,r,delta,1.2*sigma[0])*x[n][0];
//		if(k==3) return EB.call_diff((double)dt*n,1.8*T,x[n][0],Strike,r,delta,sigma[0])*x[n][0];
////		else
//		return f1D(x[n][0], kAC-1, d);
//		(2.*(double)(k)/(double)K)


	}


//	if (D >= 2)
//		return f2D(y, kAC, d);

	printf("Error 263: etwas mit dem Aufrufen der fD stimme nicht \n");
	return -1;
}
