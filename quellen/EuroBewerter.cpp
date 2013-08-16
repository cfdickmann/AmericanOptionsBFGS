/*
 * EuroBewerter.cpp
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#include "EuroBewerter.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "Polynom.h"
#define pi 3.141592654

namespace std {

double cnd(double x) {
	return 0.5 * (1 + erf(x / sqrt(2)));
}

double min(double x, double y) {
	return x < y ? x : y;
}

double stammfkt2(int order, double x) {
	if (order % 2 != 0) {
		double p = 0;
		double fak = 1;
		for (int nn = order - 1; nn >= 0; nn -= 2) {
			p += fak * pow(x, nn);
			fak *= nn;
		}
		return -exp(-x * x / 2) * p;
	} else {

		double p = 0;
		double fak = 1;
		for (int nn = order - 1; nn >= 0; nn -= 2) {
			p += fak * pow(x, nn);
			if (nn != 0)
				fak *= nn;
		}
		return -exp(-x * x / 2) * p + fak * sqrt(pi / 2) * erf(x / sqrt(2));;
	}
	return 0;
}

double integralExpQ(Polynom* p, double g1, double g2) {
	if (p->length == 0)
		return 0;
	double summe = 0;
	for (int order = 0; order < p->length; ++order)
		summe += p->koeff[order]
				* (stammfkt2(order, g2) - stammfkt2(order, g1));

	if (g2 < g1) {
		return 0;
		printf("Error 200\n");
	}
	return summe;
}

double maxi(double x, double y) {
	return x < y ? y : x;
}


double EuroBewerter::call_diff(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double f=call(t,T,X0,Strike,r,delta,sigma);

	double fh=call(t,T,X0+0.0001,Strike,r,delta,sigma);
	return (fh-f)/0.0001;
	}


// mein call
double EuroBewerter::call(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	r = r - delta;
	T = T - t;
	double d1 = (log(X0 / Strike) + (r + sigma * sigma / 2) * T)
			/ (sigma * sqrt(T));
	double d2 = d1 - sigma * sqrt(T);
	return exp(-delta * T) * exp(-r * t)
			* (X0 * cnd(d1) - Strike * exp(-r * T) * cnd(d2));
}


////Call vom Nikolaus
//double EuroBewerter::call(double t, double T, double X0, double Strike,
//		double r, double delta, double sigma) {
//
//	T = T - t;
//	double d1 = (log(X0 / Strike) + (r + sigma * sigma / 2) * T)
//			/ (sigma * sqrt(T));
//	double d2 = d1 - sigma * sqrt(T);
//	return exp(r * t)
//			* (X0 * cnd(d1) - Strike  * cnd(d2));
//}
////
//double EuroBewerter::call_diff(double t, double T, double X0, double Strike,
//		double r, double delta, double sigma) {
//
//	T = T - t;
//	double d1 = (log(X0 / Strike) + (r + sigma * sigma / 2) * T)
//			/ (sigma * sqrt(T));
//
//	return exp(r * t)
//			* ( cnd(d1));
//}


//double EuroBewerter::call_diff(double t, double T, double X0, double Strike,
//		double r, double delta, double sigma) {
//return
//}


double EuroBewerter::put(double t, double T, double X0, double Strike, double r,
		double delta, double sigma) {
	r = r - delta;
	T = T - t;
	double d1 = (log(X0 / Strike) + (r + sigma * sigma / 2) * T)
			/ (sigma * sqrt(T));
	double d2 = d1 - sigma * sqrt(T);
	return exp(-delta * T) * (Strike * exp(-r * T) * cnd(-d2) - X0 * cnd(-d1));
}

double EuroBewerter::max_call(double t, double T, double* X0, int D,
		double Strike, double r, double delta, double sigma) {
	T = T - t;
	double summe = 0;

	for (int d = 0; d < D; ++d) {
		double d_minus = (log(X0[d] / Strike)
				+ (r - delta - sigma * sigma / 2.) * T) / (sigma * sqrt(T));
		double d_plus = d_minus + sigma * sqrt(T);
		Polynom ganz;
		double eins[1] = { 1. };
		ganz.set(eins, 1);
		double I;
		double max = -100000;
		double min = 10000000;
		for (int dd = 0; dd < D; ++dd)
			if (dd != d) {
//				printf("\nd_plus %f, d_minus %f,", d_plus, d_minus);
//				printf("verschiebefaktor %f\n",
//						sigma * sqrt(T)
//								+ log(X0[d] / X0[dd]) / (sigma * sqrt(T)));
//				double pp[18] = { 0.5, 0.3960985, 0, -0.061485, 0, 0.007456, 0,
//						-5.84946E-4, 0, 2.920034E-5, 0, -9.15823E-7, 0,
//						1.740319E-8, 0, -1.826093E-10, 0, 8.10495E-13 };
				double pp[10] = { 0.50000000000009, 0.38567951086190133, 0,
						-0.05010672697589501, 0, 0.004103597701237448, 0,
						-1.631010612321749E-4, 0, 2.4428290978202304E-6 };
				Polynom p;
				p.set(pp, 10);
				double v = sigma * sqrt(T)
						+ log(X0[d] / X0[dd]) / (sigma * sqrt(T));
				max = v > max ? v : max;
				min = v < min ? v : min;

				p.verschieben(v);
//				printf("I innen%f\n",integralExpQ(&p, maxi(-d_plus, -5. + min), 5. - max));
				ganz.multiply_with(&p);
			}
		I = integralExpQ(&ganz, maxi(-d_plus, -5. + min), 5. - max);
//		printf("I aussen%f\n", I);
		summe += X0[d] * exp(-delta * T) / sqrt(2 * 3.141592654) * I;
	}
	double prod = 1;
	for (int d = 0; d < D; ++d) {
		double d_minus = (log(X0[d] / Strike)
				+ (r - delta - sigma * sigma / 2.) * T) / (sigma * sqrt(T));
		prod *= (1 - (1 - cnd(-d_minus)));
	}

	double zu = -Strike * exp(-r * T) + Strike * exp(-r * T) * prod;
	double ergebnis = (summe + zu) * exp(-r * t); //return erf(0.1);

	double e1 = 0;
	for (int d = 0; d < D; ++d)
		e1 = maxi(e1, call(t, T, X0[d], Strike, r, delta, sigma));
	return maxi(ergebnis, e1);
}

//double EuroBewerter::min_put(double t, double T, double* X0, int D,
//		double Strike, double r, double delta, double sigma) {
//	return max_call(double t, double T, double* X0, int D,
//			double Strike, double r, double delta, double sigma)-;
//}

EuroBewerter::EuroBewerter() {
	// TODO Auto-generated constructor stub

}

EuroBewerter::~EuroBewerter() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
