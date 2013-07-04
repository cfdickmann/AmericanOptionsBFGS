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

double AmericanOption::f1D(double x, int k, int d) {
	return trig(x, k);
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

double AmericanOption::f(int k, double **x, int n, int d) {
// if (k == 2)
// return 0;
//	return 1.;
	if (k / KpI != n / NpI)
		return 0;

	double y[D];
	for (int dd = 0; dd < D; ++dd)
		y[dd] = x[n][dd] / (2. * Strike);

	int kAC = k % KpI;

	if (D == 1)
		return f1D(y[0], kAC, d);
	if (D >= 2)
		return f2D(y, kAC, d);

	printf("Error 263: etwas mit dem Aufrufen der fD stimme nicht \n");
	return -1;
}
