#include "AmericanOption.h"
#include <ctype.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include "RNG.h"

using namespace std;
//using namespace alglib;

AmericanOption::AmericanOption() {
	BFGS_Nesterov_Intervals = 5;
	Iterations_Nummer = 0;
 nesterov = false;
	bfgs =  longstaffschwarz = false;
	Daten();
	neueExerciseDates(Dates);
	K = BFGS_Nesterov_Intervals * (KpI); //Anzahl der Basisfunktionen
	dt = T / (double) (N - 1);
	Anfangszeit = time(NULL);
	//(double*)malloc(sizeof(double)*Threadanzahl);
}

AmericanOption::~AmericanOption() {

}
//
//double AmericanOption::unif() {
//	//return MT();
//	return (double) (rand()) / (double) (RAND_MAX);
//}

double AmericanOption::BoxMuller(double U1, double U2) {
	double R = -2 * log(U1);
	double V = 2 * 3.1415926535 * U2;
	return sqrt(R) * cos(V);
}
//
//double AmericanOption::nextGaussian() {
//	//return qnorm(MT());//MersenneTwister benutzen
//	//return qnorm( (double)(random())/(double)(RAND_MAX) );
//
//	// Box-Muller algorithm
//	double U1 = unif();
//	double U2 = unif();
//	if (U1 == 0)
//		U1 = 0.00000001;
//	double R = -2 * log(U1);
//	double V = 2 * 3.1415926535 * U2;
//	return sqrt(R) * cos(V);
//
//	//----------------------- Die Polar-Methode von George Marsaglia
//	//	double a1,a2,q=2;
//	//	while(!((0 < q) && (q <= 1))){
//	//		a1 = 2 * MT() - 1.;  // "Zufallszahl" liefert in [0,1)
//	//		a2 = 2 * MT() - 1.;  //   gleichverteilte Werte
//	//		q = a1 * a1 + a2 * a2;
//	//	}
//	//	double p = sqrt (-2 * log(q) / q);
//	//	double z1 = a1 * p;
//	////	double z2 = a2 * p;
//	//	return z1;
//}

double AmericanOption::max(double d1, double d2) {
	if (d1 < d2)
		return d2;
	else
		return d1;
}

double AmericanOption::payoff(double** x, int time) {
	return payoff(x[time], time);
}

double AmericanOption::payoff(double* x, int time) {
	//	return payoff(x,time,D);
	//return x[0]-Strike;
	if (option == MIN_PUT)
		return max(Strike - Min(x, D), 0) * exp(-r * dt * (double) (time)); //Min Put
	if (option == MAX_CALL)
		return max(Max(x, D) - Strike, 0) * exp(-r * dt * (double) (time)); //Max Call
	printf("ERROR, option unknown!\n");
	return -1;
}

//double* AmericanOption::payoffAbl(double* x, int time) {
//	double* grad = DoubleFeld(D);
//	if (option == MIN_PUT) {
//		for (int d = 0; d < D; ++d)
//			if (x[d] == Min(x, D) && x[d] <= Strike)
//				grad[d] = -1. * exp(-r * dt * (double) (time));
//		return grad;
//	}
//	if (option == MAX_CALL) {
//		for (int d = 0; d < D; ++d)
//			if (x[d] == Max(x, D) && x[d] >= Strike)
//				grad[d] = 1. * exp(-r * dt * (double) (time));
//		return grad;
//	}
//	printf("ERROR, option unknown!\n");
//	return NULL;
//}

//
//double AmericanOption::payoff(double* x, int time, int D) {
//	if (option == MIN_PUT)
//		return max(Strike - Min(x, D), 0) * exp(-r * dt * (double) (time)); //Min Put
//	if (option == MAX_CALL)
//		return max(Max(x, D) - Strike, 0) * exp(-r * dt * (double) (time)); //Max Call
//	printf("ERROR, option unknown!\n");
//	return -1;
//}

void AmericanOption::neueExerciseDates(int n) {
	Exercise_Dates = new int[n + 1];
	number_of_Exercise_Dates = n;

	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		Exercise_Dates[e] =  (N - 1) *  e / (number_of_Exercise_Dates - 1);

	printf("neue Exercise dates: ");
	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		printf("%f, ", (double) Exercise_Dates[e]);

	printf("\n");
}
//
//void AmericanOption::Pfadgenerieren(double** X) {
//	double** wdiff = DoubleFeld(N, D);
//	for (int n = 0; n < N; ++n)
//		for (int j = 0; j < D; ++j)
//			wdiff[n][j] = sqrt(dt) * nextGaussian();
//	Pfadgenerieren(X, wdiff);
//	deleteDoubleFeld(wdiff, N, D);
//}
//
//void AmericanOption::Pfadgenerieren(double** X, int start, double * S) {
//	double** wdiff = DoubleFeld(N, D);
//
//	for (int n = 0; n < N; ++n)
//		for (int j = 0; j < D; ++j)
//			wdiff[n][j] = sqrt(dt) * nextGaussian();
//
//	Pfadgenerieren(X, wdiff, start, S);
//	deleteDoubleFeld(wdiff, N, D);
//}

void AmericanOption::Pfadgenerieren(double** X, double** wdiff) {
	Pfadgenerieren(X, wdiff, 0, X0);
}

void AmericanOption::Pfadgenerieren(double** X, int start, double* S,
		RNG* generator) {
	double** wdiff = DoubleFeld(N, D);
	for (int n = 0; n < N; ++n)
		for (int d = 0; d < D; ++d)
			wdiff[n][d] = sqrt(dt) * generator->nextGaussian();
	Pfadgenerieren(X, wdiff, 0, S);
	deleteDoubleFeld(wdiff, N, D);
}

void AmericanOption::Pfadgenerieren(double** X, double** wdiff, int start,
		double * S) {
	for (int d = 0; d < D; ++d)
		X[start][d] = S[d];

	for (int d = 0; d < D; ++d) {
		for (int n = start + 1; n < N; ++n) {
			if (PfadModell == ITO)
				X[n][d] = X[n - 1][d]
						* exp(
								(((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt
										+ sigma[d] * wdiff[n][d]));
			//X[n][d]=X[n-1][d]+wdiff[n][d];
//
//			if (PfadModell == EULER)
//				X[n][d] = X[n - 1][d] + (r - delta) * X[n - 1][d] * dt
//						+ sigma[d] * X[n - 1][d] * wdiff[n][d];
//			if (PfadModell == CIR)
//				X[n][d] = max(
//						X[n - 1][d] + kappa * (theta - X[n - 1][d]) * dt
//								+ sigma[d] * sqrt(X[n - 1][d]) * wdiff[n][d],
//						0); //mean reversion
//			if (PfadModell == ITOrho) {
//				if (D == 2) {
////				     [,1]      [,2]
////				[1,]  1.0 0.0000000
////				[2,]  0.3 0.9539392
//					double z[2];
//					z[1] = wdiff[n][0];
//					z[0] = 0.3 * wdiff[n][0] + 0.9539392 * wdiff[n][1];
//					X[n][d] = X[n - 1][d]
//							* exp(
//									(((r - delta) - 0.5 * sigma[d] * sigma[d])
//											* dt + sigma[d] * z[d]));
//				}
//				if (D == 3) {
////				     [,1]      [,2]     [,3]
////				[1,]  1.0 0.0000000 0.000000
////				[2,]  0.3 0.9539392 0.000000
////				[3,]  0.3 0.2201398 0.928191
//					double z[3];
//					z[0] = wdiff[n][0];
//					z[1] = 0.3 * wdiff[n][0] + 0.9539392 * wdiff[n][1];
//					z[2] = 0.3 * wdiff[n][0] + 0.2201398 * wdiff[n][1]
//							+ 0.928191 * wdiff[n][2];
//					X[n][d] = X[n - 1][d]
//							* exp(
//									(((r - delta) - 0.5 * sigma[d] * sigma[d])
//											* dt + sigma[d] * z[d]));
//				}
//			}
//			//if (PfadModell == JDI)
//			//	X[n][j] = X[n - 1][j] * exp(((r - delta) - 0.5 * sigma[j] * sigma[j]) * dt + sigma[j] * wdiff[n][j]) * exp(sprue[n][j]);
//			//	//			X[n][j] = max( X[n - 1][j] + (r-delta) *X[n-1][j]*dt + sigma[j] *X[n-1][j]*wdiff[n][j] +X[n - 1][j] *sprue[n][j],0);
		}
		if (X[N - 1][0] <= 0)
			printf("Error0\n");
	}
}

//double AmericanOption::newSprung() {
//	return nextGaussian() * eta/*-0.5*eta*eta*/;
//}
//
//int AmericanOption::Poisson(double theta) {
//	double p = exp(-theta);
//	double F = p;
//	double N = 0;
//	double U = MT();
//	while (U > F) {
//		N++;
//		p = p * theta / N;
//		F += p;
//	}
//	return N;
//}
//
//double AmericanOption::EuropeanPut1D_discounted(double t, double T, double S,
//		double strike) {
//	if (t >= T)
//		return max(strike - S, 0);
//	double d1 = (log(S / strike) + (r + sigma[0] * sigma[0] / 2.) * (T - t))
//			/ (sigma[0] * sqrt(T - t));
//	double d2 = d1 - sigma[0] * sqrt(T - t);
//	double P = CumulativeNormalDistribution(-d2) * strike * exp(-r * (T - t))
//			- CumulativeNormalDistribution(-d1) * S;
//	return P * exp(-t * r);
//}
//
//double AmericanOption::EuropeanCall1D_discounted(double t, double T, double S,
//		double strike) {
//	if (t >= T)
//		return max(S - strike, 0);
//	double P = EuropeanPut1D(t, T, S, strike);
//	double B = exp(-r * (T - t));
//	double C = P + S - strike * B;
//	return C * exp(-t * r);
//}
//
//double AmericanOption::EuropeanPut1D(double t, double T, double S,
//		double strike) {
//	if (t >= T)
//		return max(strike - S, 0);
//	double d1 = (log(S / strike) + (r + sigma[0] * sigma[0] / 2.) * (T - t))
//			/ (sigma[0] * sqrt(T - t));
//	double d2 = d1 - sigma[0] * sqrt(T - t);
//	double P = CumulativeNormalDistribution(-d2) * strike * exp(-r * (T - t))
//			- CumulativeNormalDistribution(-d1) * S;
//	return P;
//}
//
//double AmericanOption::EuropeanCall1D(double t, double T, double S,
//		double strike) {
//	if (t >= T)
//		return max(S - strike, 0);
//	double P = EuropeanPut1D(t, T, S, strike);
//	double B = exp(-r * (T - t));
//	double C = P + S - strike * B;
//	return C;
//}
//
//double AmericanOption::EuropeanOption1D_discounted(double t, double T, double S,
//		double strike) {
//	if (option == MIN_PUT)
//		return EuropeanPut1D_discounted(t, T, S, strike);
//	if (option == MAX_CALL)
//		return EuropeanCall1D_discounted(t, T, S, strike);
//	printf("Error 378");
//	return 0;
//}
//
//double AmericanOption::european_MaxCall_ND(double* x, double t, double T) {
//	/*//if(k==0) return 1;
//	 //	double a=((double)(k)+0.5)/(double)K*180.;
//	 double d_minus[D];
//	 double d_plus[D];
//	 for (int d = 0; d < D; ++d) {
//	 d_minus[d] = (log(x[d] / Strike)+(r - delta - sigma[d] * sigma[d] / 2.)*(T - t)) / sigma[d] / sqrt(T - t + 0.000001);
//	 d_plus[d] = d_minus[d] + sigma[d] * sqrt(T - t);
//	 }
//
//	 double erg = 0;
//	 for (int l = 0; l < D; ++l) {
//	 double integralSumme = 0;
//	 double dz = 0.01;
//	 //	printf("%f, %f, %f\n", d_plus, d_minus, t);
//	 for (double z = -8; z < d_plus[l]; z += dz) {
//	 double df = exp(-0.5 * z * z);
//	 for (int l_Strich = 0; l_Strich < D; ++l_Strich)
//	 if (l_Strich != l)df *= CumulativeNormalDistribution(log(x[l] / x[l_Strich]) / sigma[l] / sqrt(T - t + 0.000001) - z + sigma[l] * sqrt(T - t));
//	 //		if(z==-3)printf("%f, %f\n",z, df);
//	 integralSumme += df*dz;
//	 //			if(df<0.0001 && rand()%1000==0)
//	 //				printf("sehr klein bei %f\n",z);
//	 }
//
//	 erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654) * integralSumme;
//	 }
//	 double prod = 1;
//	 for (int l = 0; l < D; ++l)
//	 prod *= (1 - CumulativeNormalDistribution(d_minus[l]));
//	 erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
//	 return exp(-r * t)*(erg);*/
//	//if(k==0) return 1;
//	//	double a=((double)(k)+0.5)/(double)K*180.;
//	double d_minus[D];
//	double d_plus[D];
//	for (int d = 0; d < D; ++d) {
//		d_minus[d] = (log(x[d] / Strike)
//				+ (r - delta - sigma[d] * sigma[d] / 2.) * (T - t)) / sigma[d]
//				/ sqrt(T - t + 0.000001);
//		d_plus[d] = d_minus[d] + sigma[d] * sqrt(T - t);
//	}
//
//	double erg = 0;
//	for (int l = 0; l < D; ++l) {
//		double integralSumme = 0;
//		double dz = 0.01;
//		//	printf("%f, %f, %f\n", d_plus, d_minus, t);
//		for (double z = -8; z < d_plus[l]; z += dz) {
//			double df = exp(-0.5 * z * z);
//			for (int l_Strich = 0; l_Strich < D; ++l_Strich)
//				if (l_Strich != l)
//					df *= CumulativeNormalDistribution(
//							log(x[l] / x[l_Strich]) / sigma[l]
//									/ sqrt(T - t + 0.000001) - z
//									+ sigma[l] * sqrt(T - t));
//			//		if(z==-3)printf("%f, %f\n",z, df);
//			integralSumme += df * dz;
//			//			if(df<0.0001 && rand()%1000==0)
//			//				printf("sehr klein bei %f\n",z);
//		}
//
//		erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654)
//				* integralSumme;
//	}
//	double prod = 1;
//	for (int l = 0; l < D; ++l)
//		prod *= (1 - CumulativeNormalDistribution(d_minus[l]));
//	erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
//	return exp(-r * t) * (erg);
//}
//
//double AmericanOption::europeanValue(double* x, double t, double T) {
//	if (D >= 2) {
//		if (option == MAX_CALL)
//			return european_MaxCall_ND(x, t, T);
//	}
//	if (D == 1) {
//		if (option == MAX_CALL)
//			return EuropeanCall1D_discounted(t, T, x[0], Strike);
//		if (option == MIN_PUT)
//			return EuropeanPut1D_discounted(t, T, x[0], Strike);
//	}
//	printf("Error17\n");
//	return 0;
//}

