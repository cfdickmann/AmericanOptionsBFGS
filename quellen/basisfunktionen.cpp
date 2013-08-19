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
//	double asymp=0;
	double gewicht = 1. / (1. + 0.1 * pow(x - g, 2));
	return gewicht * f + (1. - gewicht) * asymp;
}

double Sum(double* x, int d) {
	double s = 0;
	for (int i = 0; i < d; ++i)
		s += x[i];
	return s;
}
//
//double AmericanOption::F(double x, int k, double border, bool hauf) {
//	if (k == 0)
//		return 1;
//	double extra = (x >= border && hauf) ? 1 : 0;
//	if (k % 2 == 0)
//		return cos((double) k / 2. * x) * (fabs(x) < border) + extra;
//	else
//		return sin(((double) k + 1.) * x / 2.) * (fabs(x) < border) + extra;
//}
//
//double abgeschnitten(double x, double f){
//	if(x>120)return 1;
//	if(x<50)return 0;
//	return f;
//}

double AmericanOption::f1D(double x, int k,int n, int d) {
	if(k==0)return EB.call_diff(n*dt,1.00001*T,x,Strike,r,delta,sigma[0])*x;
//		if(k==1)return EB.call_diff(n*dt,1.1*T,x,Strike,r,delta,sigma[0])*x[n][0];

return eingepasst(x,trig(0.01*x,k))*x;

	//	return pow(x/X0[0],k);
	//		return grund;
//			if(kAC==0)return (x[n][0]>Strike?1:0)*x[n][0];
			//		return eingepasst(x[n][0],trig(0.01 * x[n][0], kAC)) * x[n][0];
//			return eingepasst(x[n][0],trig(0.01 * x[n][0], kAC)) * x[n][0];
//			return (x>Strike*3.*(double)(k)/(double)(KpI) )*x ;
	//				&& x[n][0]<=Strike*3.*(double)(kAC+1)/(double)(KpI) )*x[n][0] ;
	//		return (x[n][0]>Strike*0.8+0.4*Strike*(double)(kAC-10)/(double)(KpI-10))*x[n][0] ;
	//return (trueE(x[n][0]+0.000001,n)-trueE(x[n][0],n))/0.000001*x[n][0];
	//		if(kAC==0)
	//			return EB.call_diff(n*dt,1.00001*T,x[n][0],Strike,r,delta,sigma[0])*x[n][0];
	//		return 1./(1.+0.01*pow(x[n][0]-Strike,2))*trig(0.01*(x[n][0]-Strike),kAC);
	//		if(kAC==1)return EB.call_diff(n*dt,1.05*T,x[n][0]+1.,Strike,r,delta,sigma[0])*x[n][0];
}

double AmericanOption::f2D(double* x, int k, int n, int d) {
	if (k < KpI / 4) { //Funktionen, die vom Maximum abhaengen
		if (x[d] == Max(x, D))
			return   eingepasst(x[d],trig(0.01*x[d], k))   * x[d];
		else
			return 0;
	}

	if (k < 2* KpI / 4) { //Funktionen, die von den einzelnen assets abhaengen
		k -= KpI / 4;
		if ((k % KpI) % D == d)
			return  eingepasst(x[d],trig(0.01*x[d], k / D)) * x[d];
		else
			return 0;
	}

	if (k < 3 * KpI / 4) { //Funktionen, die von der Summe abhaengen
		k -= 2 * KpI / 4;
		double S = Sum(x, D) / (double) (D);
		return eingepasst(S, trig(0.01 * S, k)) * x[d];
	}

	if (k < 4 * KpI / 4) { //Funktionen, die von der Differenz abhaengen
//		k -= 3 * KpI / 4;
//		if (x[d] == Max(x, D))
//			return EB.call_diff(n * dt, 1.00001 * T + k * 0.1 * T, x[d], Strike,
//					r, delta, sigma[0]) * x[d];
//		else
//			return 0;

		k -= 3 * KpI / 4;
		double faktor = -1;
		double max = Max(x, D);
		if (x[d] == max)
			faktor = 1;
		double S = max - (Sum(x, D) - max) / (double) (D - 1);
		return eingepasst(S,trig(0.01*S, k) * faktor) * x[d];
	}

	printf("Error 477: etwas mit der Abzaehlung der k stimmt nicht.");
	return -1;
}

double AmericanOption::trueE(double x, int n) {
	int ex = 0;
	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		if (Exercise_Dates[e] <= n)
			ex = e;

	double xx[1];
	xx[0] = x;
	double C = LSM_C_estimated(xx, ex);

	double P = 0 * (x - Strike) > 0 ? (x - Strike) : 0;
	return C > P ? C : P;
}

double pos(double x) {
	return x > 0 ? x : 0;
}

double AmericanOption::grund(double* x, int n){
	double max = exp(-r * (double) n * dt) * payoff(x, n);
	double max_diff = x[0] > Strike ? exp(-r * (double) n * dt) * 1 : 0;
	for (int nn = n; nn < N; nn++) {
		double eur = EB.call((double) nn * dt, 1.001 * T, x[0], Strike, r,
				delta, sigma[0]);
		if (eur > max) {
			max = eur;
			max_diff = EB.call_diff((double) nn * dt, 1.001 * T, x[0],
					Strike, r, delta, sigma[0]);
		}
	}

	return sigma[0]*max_diff * x[0];
}


double AmericanOption::f(int k, double **x, int n, int d) {
	if (k / KpI != n / NpI)
		return 0;

	int kAC = k % KpI;

	if (D == 1)
		return f1D(x[n][0],kAC,n,d);

	if (D >= 2)
		return f2D(x[n], kAC,n, d);

	printf("Error 263: etwas mit dem Aufrufen der fD stimme nicht \n");
	return -1;
}
