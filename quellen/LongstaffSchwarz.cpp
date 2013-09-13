//#include <iostream>
//#include <math.h>
//#include <stdlib.h>
#include "MTRand.h"
//#include "math.h"
//#include "stdlib.h"
#include "AmericanOption.h"
#include "Hilfsmittel.h"

using namespace std;

AmericanOption* zeiger;

void* DELEGATE_LSM_mittelwert(void* data) {
	zeiger->LSM_mittelwert(((int*) data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void AmericanOption::LongstaffSchwartz() { //TODO
	LSM_Mtraining = 1000;

	LSM_setting();

	X = DoubleFeld(LSM_Mtraining, N, D);
	V = DoubleFeld(LSM_Mtraining, number_of_Exercise_Dates);
	B = DoubleFeld(Mphi, Mphi);
	BV = DoubleFeld(Mphi);
	betas = DoubleFeld(N, Mphi);

	RNG generator;
	for (int m = 0; m < LSM_Mtraining; ++m)
		Pfadgenerieren(X[m], 0, X0, &generator);

	betas = new double*[number_of_Exercise_Dates];

	for (int m = 0; m < LSM_Mtraining; ++m)
		V[m][number_of_Exercise_Dates - 1] = payoff(X[m],
				Exercise_Dates[number_of_Exercise_Dates - 1]);

	for (int lauf = number_of_Exercise_Dates - 2; lauf >= 1; --lauf) {
		printf("Longstaff und Schwartz ... %d      \r", lauf);
		cout.flush();

		int Elauf = Exercise_Dates[lauf];

		for (int r = 0; r < Mphi; ++r) {
			for (int q = r; q < Mphi; ++q) {
				double erg = 0;
				for (int m = 0; m < LSM_Mtraining; ++m)
					erg = erg
							+ LSM_phi(X[m][Elauf], r, Elauf)
									* LSM_phi(X[m][Elauf], q, Elauf);
				B[r][q] = erg / (double) LSM_Mtraining;
				B[q][r] = B[r][q];
			}

			double erg2 = 0;

			for (int m = 0; m < LSM_Mtraining; ++m)
				erg2 = erg2 + LSM_phi(X[m][Elauf], r, Elauf) * V[m][lauf + 1];
			BV[r] = erg2 / (double) LSM_Mtraining;
		}

		// betas[lauf]=gausseidel(B,Mphi,BV);
		betas[lauf] = gauss(B, BV, Mphi);

//		for (int i = 0; i < Mphi; ++i)
//			printf("%f, ", betas[lauf][i]);
//		printf("\n");

		for (int m = 0; m < LSM_Mtraining; ++m) {
			V[m][lauf] = max(payoff(X[m], Elauf),
					(double) LSM_C_estimated(X[m][Elauf], lauf));
		}

	} // Schritte zu Ende
	printf("Longstaff und Schwartz ... fertig        \n");

	double high = 0;
	for (int m = 0; m < LSM_Mtraining; ++m)
		high += V[m][1] / (double) LSM_Mtraining;
	printf("high: %f\n", high);

	deleteDoubleFeld(X, LSM_Mtraining, N, D);
	deleteDoubleFeld(V, LSM_Mtraining, number_of_Exercise_Dates);
	deleteDoubleFeld(B, Mphi, Mphi);
	deleteDoubleFeld(BV, Mphi);
}

double AmericanOption::LSM_C_estimated(double* x, int ex) {
	double erg = 0;
	for (int k = 0; k < Mphi; ++k)
		erg = erg + betas[ex][k] * LSM_phi(x, k, Exercise_Dates[ex]);
	return erg;
}

double AmericanOption::LSM_phi(double* x, int j, int time) {
	if (D == 1) {
		if (j == 0)
			return 1;
		if (j == 1)
			return x[0];
		if (j == 2)
			return x[0] * x[0];
		if (j == 3)
			return x[0] * x[0] * x[0];
		if (j == 4)
			return x[0] * x[0] * x[0] * x[0];
		if (j == 5)
			return payoff(x, time);
	}

	if (D == 2) {
		if (j == 0)
			return 1;
		if (j == 1)
			return x[0];
		if (j == 2)
			return x[0] * x[1];
		if (j == 3)
			return x[1];
		if (j == 4)
			return x[0] * x[0];
		if (j == 5)
			return x[1] * x[1];
		if (j == 6)
			return payoff(x, time);
	}

	if (mlsm) {
		if (D > 2) {
			int reihe[D];
			for (int jj = 0; jj < D; ++jj)
				reihe[jj] = jj;
			BubbleSort(x, reihe, D);

			if (j < LSM_K0)
				return 1;

			if (j < LSM_K0 + LSM_K1 && j >= LSM_K0)
				return pow(x[reihe[0]], j - LSM_K0 + 3);

			if (j < LSM_K0 + LSM_K1 + LSM_K2 && j >= LSM_K0 + LSM_K1) {
				int a = (j - LSM_K0 - LSM_K1) % D;
				int b = ((j - LSM_K0 - LSM_K1) - a) / D;
				//if(verbose)printf("gemischt %d,%d\n",a,b);
				return pow(x[reihe[a]], b + 1);
			}

			if (j < LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3
					&& j >= LSM_K0 + LSM_K1 + LSM_K2)
				return x[reihe[j - LSM_K1 - LSM_K2 - LSM_K0]]
						* x[reihe[j - LSM_K1 - LSM_K2 - LSM_K0 + 1]];

			if (j < LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4
					&& j >= LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3) {
				double product = 1;
				for (int jj = 0; jj < D; ++jj)
					product *= x[reihe[jj]];
				return product;
			}
		}
	} else { // not mlsm
		if (D > 2) {
			int reihe[D];
			for (int jj = 0; jj < D; ++jj)
				reihe[jj] = jj;
			BubbleSort(x, reihe, D);

			if (j < LSM_K0)
				return 1;

			if (j < LSM_K0 + LSM_K1 && j >= LSM_K0) {
				if (j - LSM_K0 == 0)
					return x[reihe[0]];
				if (j - LSM_K0 == 1)
					return x[reihe[1]];
				if (j - LSM_K0 == 2)
					return pow(x[reihe[0]], 2);
				if (j - LSM_K0 == 3)
					return pow(x[reihe[1]], 2);
				if (j - LSM_K0 == 4)
					return pow(x[reihe[0]], 3);
				if (j - LSM_K0 == 5)
					return pow(x[reihe[1]], 3);
				if (j - LSM_K0 == 6)
					return pow(x[reihe[0]], 4);
				if (j - LSM_K0 == 7)
					return pow(x[reihe[1]], 4);
				if (j - LSM_K0 == 8)
					return pow(x[reihe[0]], 5);
				if (j - LSM_K0 == 9)
					return pow(x[reihe[1]], 5);
			}

			if (j < LSM_K0 + LSM_K1 + LSM_K2 && j >= LSM_K0 + LSM_K1) {
				if (j - LSM_K0 - LSM_K1 == 0)
					return pow(x[reihe[0]], 1) * pow(x[reihe[1]], 1);
				if (j - LSM_K0 - LSM_K1 == 1)
					return pow(x[reihe[0]], 1) * pow(x[reihe[1]], 2);
				if (j - LSM_K0 - LSM_K1 == 2)
					return pow(x[reihe[0]], 2) * pow(x[reihe[1]], 1);
				if (j - LSM_K0 - LSM_K1 == 3)
					return pow(x[reihe[0]], 2) * pow(x[reihe[1]], 2);
				if (j - LSM_K0 - LSM_K1 == 4)
					return pow(x[reihe[0]], 3) * pow(x[reihe[1]], 1);
				if (j - LSM_K0 - LSM_K1 == 5)
					return pow(x[reihe[0]], 3) * pow(x[reihe[1]], 2);
				if (j - LSM_K0 - LSM_K1 == 6)
					return pow(x[reihe[0]], 2) * pow(x[reihe[1]], 3);
				if (j - LSM_K0 - LSM_K1 == 7)
					return pow(x[reihe[0]], 1) * pow(x[reihe[1]], 3);
				if (j - LSM_K0 - LSM_K1 == 8)
					return pow(x[reihe[0]], 4) * pow(x[reihe[1]], 1);
				if (j - LSM_K0 - LSM_K1 == 9)
					return pow(x[reihe[0]], 1) * pow(x[reihe[1]], 4);
			}
		}
	}

	printf("ERROR56 %d\n", j);
	return 0;
}

double GBM(double drift, double sigma, double t, double normal) {
	return exp((drift - 0.5 * sigma * sigma) * t + sqrt(t) * sigma * normal);
}

void AmericanOption::LSM_mittelwert(int threadnummer) {
//	int mAnfang = M * threadnummer / Threadanzahl;
//	int mEnde = M * (threadnummer + 1) / Threadanzahl;
//
//	if (verbose)
//		printf("Thread %d: von %d bis %d\n", threadnummer, mAnfang, mEnde);
//
//	for (int r = 0; r < Mphi; ++r) {
//		for (int q = 0; q < Mphi; ++q) {
//			double erg = 0;
//			for (int m = mAnfang; m < mEnde; ++m)
//				//				if(1 ==0 && LSlauf!=0 && LSlauf!=1 && LSlauf!=N-1 && LSlauf!=N-2 && Cestimated(X[m][LSlauf+1],LSlauf+1)<= payoff( X[m][LSlauf+1],LSlauf+1))
//				//					ausgesetzt[threadnummer]++; else
//				erg = erg
//						+ LSM_phi(X[m][LSlauf], r, LSlauf)
//								* LSM_phi(X[m][LSlauf], q, LSlauf);
//			B[r][q][threadnummer] = erg / (double) M;
//		}
//
//		double erg2 = 0;
//		for (int m = mAnfang; m < mEnde; ++m)
//			//			if(1==0 && LSlauf!=0 && LSlauf!=1 && LSlauf!=N-1 && LSlauf!=N-2 && Cestimated(X[m][LSlauf+1],LSlauf+1)<= payoff( X[m][LSlauf+1],LSlauf+1))
//			//				ausgesetzt[threadnummer]++; else
//			erg2 = erg2 + LSM_phi(X[m][LSlauf], r, LSlauf) * V[m][LSlauf + 1];
//		BV[r][threadnummer] = erg2 / (double) M;
//	}
}

void AmericanOption::LSM_setting() {
	zeiger = this;  // damit threada auf diese objekt zugreifen koennen
	//Problemdaten laden

	if (D == 1)
		Mphi = 6;
	if (D == 2)
		Mphi = 7;
	if (D > 2) {
		if (mlsm)   //Anzahl der Basisfunktionen initialisieren
		{
			LSM_K0 = 1;   //Konstante
			LSM_K1 = 3;   // 3,4,5 polynom des teuersten assets
			LSM_K2 = D * 2; //Polynome zweiter ordnung in allen einzelnen assets
			if (D > 2)
				LSM_K3 = (D - 1);   //Produkte von Verfolgern
			else
				LSM_K3 = 0;
			LSM_K4 = 1;   //Produkte aller assets
		} else {
			LSM_K0 = 1;   //Konstante
			LSM_K1 = 10;
			LSM_K2 = 10;   //Polynome zweiter ordnung in allen einzelnen assets
			LSM_K3 = 0;
			LSM_K4 = 0;   //Produkte aller assets
		}
		Mphi = LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4;
	}

	//Basisfunktionen testen
	//	double tt[Mphi];
	//	for(int i=0;i<Mphi;++i)
	//		tt[i]=i%D+2;

	//	for(int t=0;t<Mphi;++t){
	//		if( t==K_0 || t==K_1+K_0 || t==K_1+K_2+K_0 || t==K_1+K_2+K_0+K_3 || t==K_1+K_2+K_0+K_3+K_4)printf("\n");
	//		printf("t=%d, %f\n",t,phi(tt,t,1));
	//	}

	//N = Testing_Dates; // fuer longstaff schwartz keine hoehere genauigkeit noetig, da explizite formel fuer pfade verwendet wurde
	//dt = T / (double(N - 1));
	//M = LSM_Mtraining;
}

