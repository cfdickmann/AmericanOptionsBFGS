/*
 * StochIntErstellen.cpp
 *
 *  Created on: Aug 22, 2013
 *      Author: cfdickmann
 */
#include "AmericanOption.h"

using namespace std;

void AmericanOption::StochInt_erstellen() {
	StochIntegrals = DoubleFeld(K, M, number_of_Exercise_Dates);
	double ** STI = DoubleFeld(K, N);
	for (int m = 0; m < M; ++m) {
		for (int k = 0; k < K; ++k)
			STI[k][0] = 0;
		if (m % 10 == 0) {
			printf("StochInt erstellen ... %d %% \r", m * 100 / M);
			cout.flush();
		}
		for (int k = 0; k < K; ++k) {
			bool schon=false;
			for (int n = 1; n < N; ++n) {
				if(X[m][n-1][0]<Strike || X[m][n-1][1]<Strike)schon=true;
				double s = 0;
					double* ergs = f_all(k, X[m], Y[m], reihe[m], n - 1);
					if (ergs != NULL) {
						for (int d = 0; d < D; ++d)
							s += ergs[d] * WDiff[m][n][d]*schon;;
						delete[] ergs;
					}
				STI[k][n] = STI[k][n - 1] + s;
			}
			for (int e = 1; e < number_of_Exercise_Dates; ++e)
				StochIntegrals[k][m][e] = STI[k][Exercise_Dates[e]];
		}
	}
	deleteDoubleFeld(STI, K, N);
	printf("StochInt erstellen ... fertig\n");
}

//		for (int k = 0; k < K; ++k) {
//			int ex = -1;
//			for (int n = 1; n < N; ++n) {
//				if (ex == -1)
//					if (payoff(X[m], n - 1) > 0)
//						ex = n;
//
//				double s = 0;
//				if (ex != -1)
//					for (int d = 0; d < D; ++d)
//						s += f(k, X[m], n - 1, d) * WDiff[m][n][d];
//
//				//					s = EB.put_diff((double) n * dt, T, X[m][n - 1][0], Strike,
//				//							r, delta, sigma[0])*X[m][n-1][0]* WDiff[m][n][0];
//
//				//										s = EB.put((double) (n) * dt, T, X[m][n][0], Strike,
//				//												r, delta, sigma[0])
//				//										- EB.put((double) (n-1) * dt, T, X[m][n - 1][0], Strike,
//				//												r, delta, sigma[0]);
//
////						int ex=  (n-1)/(N/Dates);
////						if( LSM_C_estimated(X[m][Exercise_Dates[ex]],ex)*1.25<payoff(X[m],Exercise_Dates[ex])  )
////							s=-X[m][n-1][0]* WDiff[m][n][0];
//
//				StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + s;
//			}
//		}

//		for (int k = 0; k < K; ++k) {
////			int ex = -1;
//			for (int n = 1; n < N; ++n) {
////				if (ex == -1)
////					if (payoff(X[m], n - 1) > 0)
////						ex = n;
//				double summe = 0;
////				if (ex != -1)
//				for (int d = 0; d < D; ++d)
//					summe += f(k, X[m], n - 1, d) * WDiff[m][n][d];
//				StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + summe;
//			}
//		}

//		for (int k = 0; k < K; ++k) {
//			for (int n = 1; n < N; ++n) {
//				double summe = 0;
//				if (payoff(X[m], n - 1) > LSM_C_estimated(X[m],(n-1)/(N/Dates))) {
//					int d_max = argMax(X[m][n - 1], D);
////
//					summe = -1.*exp(-(double)(n-1)*dt) * X[m][n - 1][d_max] * WDiff[m][n][d_max]
//							* sigma[d_max];
//				} else
////					for (int d = 0; d < D; ++d)
////						summe += f(k, X[m], n - 1, d) * WDiff[m][n][d]
////								* sigma[d];
//					summe=EB.put_diff((double) n * dt, T, X[m][n - 1][0], Strike,
//												r, delta, sigma[0])*sigma[0]*X[m][n-1][0]*WDiff[m][n][0];
//
//				StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + summe;
//			}
//		}

//		for (int k = 0; k < K; ++k) {
//					int ex = -1;
//					for (int n = 1; n < N; ++n) {
//						if (ex == -1)
//							if (payoff(X[m], n - 1) > 0)
//								ex = n;
//						double summe = 0;
//						if (ex != -1)
//						for (int d = 0; d < D; ++d)
//							summe += (cndf((log(X[m][n-1][0]/Strike) +(r+sigma[0]*sigma[0]/2.)*(T-(double)(n-1)*dt+0.0001)   )
//									/(sigma[0]*sqrt(T-(double)(n-1)*dt+0.0001)))-1.) *X[m][n-1][0]* WDiff[m][n][d];
//						StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + summe;
//					}
//				}

//		// Urusov
//		 for (int k = 0; k < K; ++k) {
//			int zzz = -1;
//			if(payoff(X[m], N/3)==0)	zzz = 0;
//			int indikator = -1;
//			int indikator2 = -1;
//			int indikator3 = -1;
//
//			for (int n = 1; n < N; ++n) {
//				double s = 0;
//				if (indikator == -1)
//					if (payoff(X[m], n - 1) > 0) {
//						indikator = n;
//						if (indikator < (N/4))
//							indikator2 = 0;
//					}
//
////				if (indikator != -1)
////					s = EB.put_diff((double) n * dt, T, X[m][n - 1][0], Strike,
////							r, delta, sigma[0])*X[m][n-1][0]*WDiff[m][n][0];
//
//				if (indikator != -1){
//					s = EB.put((double) (n) * dt, T, X[m][n][0], Strike, r, delta,
//						sigma[0])
//						- EB.put((double) (n - 1) * dt, T, X[m][n - 1][0],
//								Strike, r, delta, sigma[0]);
//				}
//				if (indikator3 == -1 && zzz == 0 && indikator2 == 0 && n > (N/3))
//				{
//					if (payoff(X[m], n - 1) > 0)
//						indikator3 = n;
//					else s = 0;
//				}
//
//				StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + s;
//			}
//		}

//	for (int k = 0; k < K; ++k) {
//			for (int n = 1; n < N; ++n) {
//				double summe = 0;
//				for (int d = 0; d < D; ++d)
//					summe += grund(X[m][n - 1], n - 1) * WDiff[m][n][d];
//				GrundIntegrals[m][n] = GrundIntegrals[m][n - 1] + summe;
//			}
//		}

