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
//	double asymp = x < g ? 0 : 1;

	double asymp = x < g ? -1 : 0;
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
double AmericanOption::abgeschnitten(double x, double f) {
	if (option == MAX_CALL) {
		if (x > 1.0)
			return 1;
		if (x < -1.0)
			return 0;
	}
	if (option == MIN_PUT) {
		if (x > 1.0)
			return 0;
		if (x < -1.0)
			return -1;
	}
	return f;
}

bool AmericanOption::inExerciseRegion(double x, int n) {
	double xx[1];
	xx[0] = x;
	return inExerciseRegion(xx, n);
}

bool AmericanOption::inExerciseRegion(double* x, int n) {
	int exTime = (n - 1) / (N / Dates);
	if (LSM_C_estimated(x, exTime) < payoff(x, Exercise_Dates[exTime]))

		return true;
	return false;
}

double AmericanOption::f1D(double x, double y, int k, int n) {
//	bool inExRegion = inExerciseRegion(x, n);

//	double d1 = (log(x / Strike)
//			+ (r - delta + sigma[0] * sigma[0] / 2)
//					* (T - (double) n * dt + 0.0001))
//			/ (sigma[0] * sqrt(T - (double) n * dt + 0.0001));
//
//	return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//			* (-(trig(-d1, k)))) * sigma[0] * x;

//				return (x>0.3*Strike+Strike*1.*(double)(k)/(double)(KpI) )*x ;
//
	return trig(y, k) * x;

//	if (option == MAX_CALL)
//		return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//				* (cnd(d1))) * x;
//	else
//		return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//				* (-cnd(-d1))) * x;

//	if (k == 0) {
//		if (inExRegion)
//			return 0;
//		else
//			return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//					* (-cnd(-d1))) * x;
//	}
//
//	if (k == 1) {
//		if (inExRegion)
//			return -x; // Achtung nur fuer 1D put
//		else
//			return 0;
//	}

//	return eingepasst(x,pow(0.01*x,k))  *x ;

//	if (k == 0) {
//		if (option == MIN_PUT)
//			return (payoff(&x, n) > 0 ? -1 : 0) * x;
//		else
//			return (payoff(&x, n) > 0 ? 1 : 0) * x;
//	}
//
//	if (k > 0)
//		return trig(0.01 * x, k - 1) * x;

//	double d1 = (log(x / Strike) + (r - delta + sigma[0] * sigma[0] / 2) *
//			(T - (double)n*dt+0.0001))/ (sigma[0] * sqrt(T - (double)n*dt+0.0001));
//	return(exp(-(r-delta)*(double)n*dt)*exp(-T*delta)* (-cnd(-d1)))*x;
//	return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//			* (abgeschnitten(d1, trig(d1, k)))) * x;

//return EB.put_diff2(n*dt,1.00001*T,x,Strike,r,delta,sigma[0])*x;

//double t=(double)n*dt;

//EB.put_diff(n*dt,1.00001*T,x,Strike,r,delta,sigma[0])*x;
//
//double t=(double) n*dt;
//double d1=(log(x/Strike)+(r-sigma[0]*sigma[0]/2.)*(T-t+0.00001))/(sigma[0]*sqrt(T-t+0.00001));
//return (exp(-r*t)*cnd(d1)-1)*x;

////		if(k==1)return EB.call_diff(n*dt,1.1*T,x,Strike,r,delta,sigma[0])*x[n][0];
//
//return eingepasst(x,trig(0.01*x,k))*x;

//return pow(x/X0[0],k);
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

//		if (x[d] == Max(x, D))
//			return pow(x[d] / Strike, k);
//		else
//			return 0;
//	} else {
//		k -= 5;
//		double S = Sum(x, D) / (double) (D);
//		return pow(S / Strike, k);
//	}

//double AmericanOption::f2D(double* x, int k, int n, int d) {
//	int d_extrem;
//	if(option==MAX_CALL)d_extrem=argMax(x, D);
//	else d_extrem=argMin(x,D);
//
//	if (d != d_extrem)
//		return 0;
//
//	double d1[D];
//	for (int d = 0; d < D; ++d)
//		d1[d] = (log(x[d] / Strike)
//				+ (r - delta + sigma[0] * sigma[0] / 2)
//						* (T - (double) n * dt + 0.0001))
//				/ (sigma[0] * sqrt(T - (double) n * dt + 0.0001));
//	if (option == MAX_CALL) {
//		return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//				* (cnd(d1[d_extrem]))) * x[d_extrem];
//	} else {
//		return (exp(-(r - delta) * (double) n * dt) * exp(-T * delta)
//				* (-cnd(-d1[d_extrem]))) * x[d_extrem];
//	}
//}

double AmericanOption::F(double x, int k, double border, bool hauf) {
	double sigma = 0.8;
	double bound = 0.45;
	double xs = x / sigma;
	if (k % 2 == 0)
		return cos(k / 2 * xs) * (fabs(xs) < bound) + (xs >= bound);
	else
		return sin((k + 1) * xs / 2) * (fabs(xs) < bound) + (xs >= bound);
}

double AmericanOption::F2(double x, int k, double border, bool hauf) {
	double sigma = 0.4;
	double bound = 0.55;
	double xs = x / sigma;
	if (k % 2 == 0)
		return cos(k / 2 * xs) * (fabs(xs) < bound) + (xs >= bound);
	else
		return sin((k + 1) * xs / 2) * (fabs(xs) < bound) + (xs >= bound);
}

double* AmericanOption::f2D_all(double* x, double*y, int k, int n) {
	double* ergs = new double[D];

	for (int d = 0; d < D; ++d)
		ergs[d] = 0;

	double* yy = y;

	double abs = fabs(x[0] - x[1]);

//	double yy[KpI];
//		for (int kk = 0; kk < KpI; ++kk)
//			yy[kk] = log(x[kk] / Strike) / sqrt(T - (double) n * dt+0.00001);
//
	if (abs < 15.) {

		if (k < 12) {
			ergs[0] = F(yy[0] + yy[1], k, 0, 0) * x[0];
			ergs[1] = F(yy[0] + yy[1], k, 0, 0) * x[1];
//			ergs[0] = F(yy[0], k, 0, 0) * (x[0] >= x[1]) * x[0];
//			ergs[1] = F(yy[1], k, 0, 0) * (x[0] <= x[1]) * x[1];
			return ergs;
		}

		if (k < 24) {
			ergs[0] = F(yy[0], k - 12, 0, 0) * x[0];
			ergs[1] = F(yy[1], k - 12, 0, 0) * x[1];
			return ergs;
		}
		return ergs;
	}

	if (abs >= 15.) {
		if (k < 36) {
			ergs[0] = F2(yy[0], k - 24, 0, 0) * (x[0] >= x[1]) * x[0];
			ergs[1] = F2(yy[1], k - 24, 0, 0) * (x[0] <= x[1]) * x[1];
			return ergs;
		}

		if (k < 48) {
			ergs[0] = F2(yy[0], k - 36, 0, 0) * x[0];
			ergs[1] = F2(yy[1], k - 36, 0, 0) * x[1];
			return ergs;
		}
		return ergs;
	}

//
//	if (k < 36) {
//		ergs[0] = F(yy[0], k - 12, 0, 0) * x[0];
//		ergs[1] = F(yy[1], k - 12, 0, 0) * x[1];
//		return ergs;
//	}

//	if (k < 24) {
//		ergs[0] = F(yy[0] + yy[1], k - 12, 0, 0) * x[0];
//		ergs[1] = F(yy[0] + yy[1], k - 12, 0, 0) * x[1];
//		return ergs;
//	}

//	if (k < 24) {
//		int groesser0 = -1;
//		if (yy[0] > yy[1])
//			groesser0 = 1;
//
//		ergs[0] = trig(fabs(yy[0] - yy[1]), k - 12) * groesser0 * x[0];
//		ergs[1] = trig(fabs(yy[0] - yy[1]), k - 12) * (-groesser0) * x[1];
//		return ergs;
//	}

//	if (k < KpI / 4) { //Funktionen, die vom Maximum abhaengen
//		int d = argMax(x, D);
//		ergs[d] = trig(y[d], k) * x[d];
//		return ergs;
//	}
//
//	if (k < 2 * KpI / 4) { //Funktionen, die von den einzelnen assets abhaengen
//		k -= 1 * KpI / 4;
//		for (int d = 0; d < D; ++d)
//			if ((k % KpI) % D == d)
//				ergs[d] = trig(y[d], k / D) * x[d];
//		return ergs;
//	}
//
//	if (k < 3 * KpI / 4) { //Funktionen, die von der Summe abhaengen
//		k -= 2 * KpI / 4;
//		double S = Sum(x, D) / (double) (D);
//
//		for (int d = 0; d < D; ++d)
//			ergs[d] = trig(S, k) * x[d];
//		return ergs;
//	}
//
//	if (k < 4 * KpI / 4) { //Funktionen, die von der Differenz abhaengen
//		k -= 3 * KpI / 4;
//		double faktor = -1;
//		double max = Max(x, D);
//
//		for (int d = 0; d < D; ++d) {
//			if (x[d] == max)
//				faktor = 1;
//			double S = max - (Sum(y, D) - max) / (double) (D - 1);
//			ergs[d] = trig(S, k) * faktor * x[d];
//		}
//		return ergs;
//	}

	printf("Error 845\n");
	return NULL;
}

////
//double* AmericanOption::f2D_all(double* x, int k, int n) {
//	double* ergs = new double[D];
//
////	double y[KpI];
////
////	for (int kk = 0; kk < KpI; ++kk)
////		y[kk] = log(x[kk] / Strike) / sqrt(T - (double) n * dt);
//
//	int reihe[D];
//	for (int j = 0; j < D; ++j)
//		reihe[j] = j;
//	BubbleSort(x, reihe, D);
//
//	int j_extrem;
//	if (option == MAX_CALL)
//		j_extrem = reihe[0];
//	if (option == MIN_PUT)
//		j_extrem = reihe[D - 1];
//
//	if (k < KpI / 4) { //Funktionen, die vom Maximum abhaengen
////			if (x[d] == Max(x, D))
////				//			return eingepasst(x[d], trig( x[d], k)) * x[d];
////				if (option == MAX_CALL)
////					j_extrem = reihe[0];
////			if (option == MIN_PUT)
////				j_extrem = reihe[D - 1];
////
////			if (j == j_extrem)
////				W = F(Xy[j], k, 0.8, true);
//		for (int d = 0; d < D; ++d)
//			ergs[d] = 0;
////						ergs[j_extrem]= trig(x[j_extrem], k) * x[j_extrem];
////		ergs[j_extrem] = F(y[j_extrem], k, 0.8, true) * x[j_extrem];
//		ergs[j_extrem] = trig(0.01*x[j_extrem], k) * x[j_extrem];
//		return ergs;
//	}
//
//	if (k < 2 * KpI / 4) { //Funktionen, die von den einzelnen assets abhaengen
//		k -= KpI / 4;
//
//
////		for (int d = 0; d < D; ++d)
////			ergs[d] = F(y[d], k, 0.8, false) * x[d];
//for(int d=0;d<D;++d)
//		if(k%D==d)
//		   ergs[d] = trig(0.01*x[d], k/D) * x[d];
//		return ergs;
//	}
//
//	if (k < 3 * KpI / 4) { //Funktionen, die von der Summe abhaengen
//		k -= 2 * KpI / 4;
//
//		double Summe = 0;
//		for (int d = 0; d < D; ++d)
//			Summe += x[d];
//		for (int d = 0; d < D; ++d)
////			ergs[d] = F(Summe, k, 1.6, false) * x[d];
//		ergs[d] = trig(0.01*Summe, k) * x[d];
//		return ergs;
//	}
//
//	if (k < 4 * KpI / 4) { //Funktionen, die von der Differenz abhaengen
//		k -= 3 * KpI / 4;
//		int j_1, j_2;
//		if (option == MAX_CALL) {
//			j_1 = reihe[0];
//			j_2 = reihe[1];
//		}
//		if (option == MIN_PUT) {
//			j_1 = reihe[D - 2];
//			j_2 = reihe[D - 1];
//		}
//
//		for (int d = 0; d < D; ++d)
//			ergs[d] = 0;
////
////		ergs[j_1] = +F(y[j_1] - y[j_2], k, 0.9, true) * x[j_1];
////		ergs[j_2] = -F(y[j_1] - y[j_2], k, 0.9, true) * x[j_2];
//
//		ergs[j_1] = trig(x[j_1] - x[j_2], k) * x[j_1];
//		ergs[j_2] = trig(x[j_1] - x[j_2], k) * x[j_2];
//		return ergs;
//	}
//
//	printf("Error 2343\n");
//	return NULL;
//}

//double* AmericanOption::f2D_all(double* x, int k, int n) {
//	double* ergs = new double[D];
//	for (int d = 0; d < D; ++d)
//		ergs[d] = 0;
//
//	int d_extrem = argMax(x, D);
////	ergs[d_extrem] = EB.max_call_diff((double) n * dt, T, x,D, Strike, r,
////			delta, sigma[d_extrem],d_extrem) * x[d_extrem] * sigma[d_extrem];
//
////	if( inExerciseRegion(x,n) && k<3){
////		return ergs;
////	}
//
//	if (k == 0) {
////	int d_extrem = argMax(x, D);
//		ergs[d_extrem] = EB.call_diff((double) n * dt, T, x[d_extrem], Strike,
//				r, delta, sigma[d_extrem]) * x[d_extrem] * sigma[d_extrem];
//
//	}
//	if (k == 1)
//		ergs[0] = EB.call_diff((double) n * dt, T, x[0], Strike, r, delta,
//				sigma[0]) * x[0] * sigma[0];
//
//	if (k == 2)
//		ergs[1] = EB.call_diff((double) n * dt, T, x[1], Strike, r, delta,
//				sigma[1]) * x[1] * sigma[1];
//
//	if (k == 3) {
////		if(inExerciseRegion(x,n))
////		ergs[d_extrem]=-x[d_extrem];
//		ergs[0] = EB.call_diff((double) n * dt, T, x[1] + x[0], Strike, r,
//				delta, sigma[0]) * x[0] * sigma[0];
//		ergs[1] = EB.call_diff((double) n * dt, T, x[1] + x[0], Strike, r,
//				delta, sigma[1]) * x[1] * sigma[1];
//	}
//
//	if (k >= 4) { //Funktionen, die von der Summe abhaengen
//		k -= 4;
//		double Summe = 0;
//		for (int d = 0; d < D; ++d)
//			Summe += x[d];
//		for (int d = 0; d < D; ++d)
//				ergs[d] = F(Summe/Strike, k, 1.6, false) * x[d];
////			ergs[d] = trig(0.01 * Summe, k) * x[d];
//		return ergs;
//	}
//
//	return ergs;
//}

//double AmericanOption::f2D(double* x, int k, int n, int d) {
//// double y[KpI];
////
//// for (int kk = 0; kk < KpI; ++kk)
//// y[kk] = log(x[kk] / Strike) / sqrt(T - (double) n * dt);
//
//
////
//// 	if (k < 5) {
//// 		if (x[d] == Max(x, D))
//// 			//			return pow(x[d] / Strike, k);
//// 						return pow(y[d], k);
////
//// 		else
//// 			return 0;
//// 	} else {
//// 		k -= 5;
////
//// 		double S = Sum(y, D) / (double) (D);
//// 		return pow(S, k);
//// //		double S = Sum(x, D) / (double) (D);
//// //		return pow(S / Strike, k);
//// 	}
//
// //	double y[KpI];
// //
// //	for (int kk = 0; kk < KpI; ++kk)
// //		y[kk] = log(x[kk] / Strike) / sqrt(T - (double) n * dt);
//
// 	if (k < KpI / 4) { //Funktionen, die vom Maximum abhaengen
// 		if (x[d] == Max(x, D))
// //			return eingepasst(x[d], trig( x[d], k)) * x[d];
// 			return trig(x[d], k) * x[d];
// 		else
// 			return 0;
// 	}
//
// 	if (k < 2 * KpI / 4) { //Funktionen, die von den einzelnen assets abhaengen
// 		k -= KpI / 4;
// 		if ((k % KpI) % D == d)
// 			return trig(x[d], k) * x[d];
// 		else
// 			return 0;
// 	}
//
// 	if (k < 3 * KpI / 4) { //Funktionen, die von der Summe abhaengen
// 		k -= 2 * KpI / 4;
// 		double S = Sum(x, D) / (double) (D);
// 		return trig(S, k) * x[d];
// 	}
//
// 	if (k < 4 * KpI / 4) { //Funktionen, die von der Differenz abhaengen
// //		k -= 3 * KpI / 4;
// //		if (x[d] == Max(x, D))
// //			return EB.call_diff(n * dt, 1.00001 * T + k * 0.1 * T, x[d], Strike,
// //					r, delta, sigma[0]) * x[d];
// //		else
// //			return 0;
//
// 		k -= 3 * KpI / 4;
// 		double faktor = -1;
// 		double max = Max(x, D);
// 		if (x[d] == max)
// 			faktor = 1;
// 		double S = max - (Sum(x, D) - max) / (double) (D - 1);
// 		return trig(S, k) * x[d];
// 	}
//
// 	printf("Error 477: etwas mit der Abzaehlung der k stimmt nicht.");
// 	return -1;
// }

//double AmericanOption::trueE(double x, int n) {
//	int ex = 0;
//	for (int e = 0; e < number_of_Exercise_Dates; ++e)
//		if (Exercise_Dates[e] <= n)
//			ex = e;
//
//	double xx[1];
//	xx[0] = x;
//	double C = LSM_C_estimated(xx, ex);
//
//	double P = 0 * (x - Strike) > 0 ? (x - Strike) : 0;
//	return C > P ? C : P;
//}
//
//double pos(double x) {
//	return x > 0 ? x : 0;
//}
//
//double AmericanOption::grund(double* x, int n){
//	double max = exp(-r * (double) n * dt) * payoff(x, n);
//	double max_diff = x[0] > Strike ? exp(-r * (double) n * dt) * 1 : 0;
//	for (int nn = n; nn < N; nn++) {
//		double eur = EB.call((double) nn * dt, 1.001 * T, x[0], Strike, r,
//				delta, sigma[0]);
//		if (eur > max) {
//			max = eur;
//			max_diff = EB.call_diff((double) nn * dt, 1.001 * T, x[0],
//					Strike, r, delta, sigma[0]);
//		}
//	}
//
//	return sigma[0]*max_diff * x[0];
//}
double* AmericanOption::f_all(int k, double **x, double** y, int n) {

	if (k / KpI != n / NpI)
		return NULL;

	int kAC = k % KpI;

	double* xx = new double[1];
	if (D == 1) {
		xx[0] = f1D(x[n][0], y[n][0], kAC, n);
		return xx;
	}
	//	return f2D(x[n], kAC, n, d);

	if (D >= 2)
		return f2D_all(x[n], y[n], kAC, n);

	printf("Error 263: etwas mit dem Aufrufen der fD stimme nicht \n");
	return NULL;
}

//
//double AmericanOption::f(int k, double **x, int n, int d) {
//	double* ergs = f_all(k, x, n);
//	if (ergs != NULL) {
//		double res = ergs[d];
//		delete[] ergs;
//		return res;
//	}
//	return 0;
//
////	if (k / KpI != n / NpI)
////		return 0;
////
////	int kAC = k % KpI;
////
////	if (D == 1)
////		return f1D(x[n][0], kAC, n, d);
////	//	return f2D(x[n], kAC, n, d);
////
////	if (D >= 2)
////		return f2D(x[n], kAC, n, d);
////
////	printf("Error 263: etwas mit dem Aufrufen der fD stimme nicht \n");
////	return -1;
//}
