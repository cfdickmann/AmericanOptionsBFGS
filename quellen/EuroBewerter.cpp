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

double dnd(double x) {
	return 1. / sqrt(2. * 3.14159265) * exp(-x * x / 2.);
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
	double f = call(t, T, X0, Strike, r, delta, sigma);

	double fh = call(t, T, X0 + 0.00001, Strike, r, delta, sigma);
	return (fh - f) / 0.00001;
}

//double EuroBewerter::put_diff(double t, double T, double X0, double Strike,
//		double r, double delta, double sigma) {
//	double f = put(t, T, X0, Strike, r, delta, sigma);
//	double fh = put(t, T, X0 + 0.000001, Strike, r, delta, sigma);
//	return (fh - f) / 0.000001;
//}

double EuroBewerter::put_diff(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2) * TStrich)
			/ (sigma * sqrt(TStrich));
	return exp(-r * t - delta * TStrich) * (cnd(d1) - 1.);
}

double EuroBewerter::put_diff2(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
//	double r_Strich = r - delta;
//		double T_Strich = T - t;
	double d1 = (log(X0 / Strike) + (r - delta + sigma * sigma / 2) * (T - t))
			/ (sigma * sqrt(T - t));
	return exp(-(r - delta) * t) * exp(-T * delta) * (-cnd(-d1));
}

//
//// mein call
//double EuroBewerter::call(double t, double T, double X0, double Strike,
//		double r, double delta, double sigma) {
//	double rStrich = r - delta;
//	double TStrich = T - t;
//	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2.) * TStrich)
//			/ (sigma * sqrt(TStrich));
//	double d2 = d1 - sigma * sqrt(TStrich);
//	double wert=exp(-delta * TStrich) * exp(-r * t)
//			* (X0 * cnd(d1) - Strike * exp(-r * TStrich) * cnd(d2));
//return wert*exp(-t*(r+delta));
//}

// Dieser hier ist getestet und gut, discounted und incl. dividend yield

double EuroBewerter::call(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2.) * TStrich)
			/ (sigma * sqrt(TStrich));
	double d2 = d1 - sigma * sqrt(TStrich);
	double wert = (X0 * cnd(d1) - Strike * exp(-rStrich * TStrich) * cnd(d2));
	wert *= exp(-(delta) * T);
	return wert * exp(-t * (r - delta));
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

//auch getestet und gut

//double EuroBewerter::put(double t, double T, double X0, double Strike, double r,
//		double delta, double sigma) {
//	r = r - delta;
//	T = T - t;
//	double d1 = (log(X0 / Strike) + (r + sigma * sigma / 2) * T)
//			/ (sigma * sqrt(T));
//	double d2 = d1 - sigma * sqrt(T);
//	return exp(-(r+delta)*t)*(exp(-delta * T) * (Strike * exp(-r * T) * cnd(-d2) - X0 * cnd(-d1)));
//}

double EuroBewerter::put(double t, double T, double X0, double Strike, double r,
		double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2) * TStrich)
			/ (sigma * sqrt(TStrich));
	double d2 = d1 - sigma * sqrt(TStrich);
	return exp(-r * t)
			* (exp(-delta * TStrich)
					* (Strike * exp(-rStrich * TStrich) * cnd(-d2)
							- X0 * cnd(-d1)));
}

double EuroBewerter::max_call_diff(double t, double T, double* X0, int D,
		double Strike, double r, double delta, double sigma, int d) {
	double f = max_call(t, T, X0, D, Strike, r, delta, sigma);
	X0[d] += 0.000001;
	double fh = max_call(t, T, X0, D, Strike, r, delta, sigma);
	X0[d] -= 0.000001;
	return (fh - f) / 0.000001;
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

//Margrabes Formel
double EuroBewerter::exchange_option(double x, double y, double t, double T,
		double r, double delta, double sigma) {
	// y to  be replaced
	double TStrich = T - t;

	double d1 = (log(x / y) + (sigma * sigma) * TStrich)
			/ (sqrt(2.) * sigma * sqrt(TStrich));
	double d2 = (log(x / y) - (sigma * sigma) * TStrich)
			/ (sqrt(2.) * sigma * sqrt(TStrich));

	double erg = x * exp(-delta * TStrich) * cnd(d1)
			- y * exp(-delta * TStrich) * cnd(d2);
	return erg * exp(-r * t);
}

//double EuroBewerter::exchange_option_diff(double x, double y, double t,
//		double T, double r, double delta, double sigma, int re) {
//	if (re == 0) {
//		double h = 0.00001;
//		double fh = exchange_option(x + h, y, t, T, r, delta, sigma);
//		double f = exchange_option(x, y, t, T, r, delta, sigma);
//		return (fh - f) / h;
//	} else {
//		double h = 0.00001;
//		double fh = exchange_option(x, y + h, t, T, r, delta, sigma);
//		double f = exchange_option(x, y, t, T, r, delta, sigma);
//		return (fh - f) / h;
//	}
//}


//auch getestet
double EuroBewerter::exchange_option_diff(double x, double y, double t,
		double T, double r, double delta, double sigma, int re) {

	double TStrich = T - t;
	double Wurzel = sqrt(2.) * sigma * sqrt(TStrich);

	double d1 = (log(x / y) + (sigma * sigma) * TStrich) / Wurzel;
	double d2 = (log(x / y) - (sigma * sigma) * TStrich) / Wurzel;

//if(rand()%100==0)
//printf("%f,\n",cnd(d1));
	if (re == 0)
		return exp(-r * t - delta * (T - t))
				* (cnd(d1) + (dnd(d1) - y / x * dnd(d2)) / Wurzel);
	else
		return exp(-r * t - delta * (T - t))
				* (-cnd(d2) + (dnd(d2) - x / y * dnd(d1)) / Wurzel);
}

EuroBewerter::EuroBewerter() {
	// TODO Auto-generated constructor stub

}

EuroBewerter::~EuroBewerter() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
