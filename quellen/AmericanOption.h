#ifndef AMERICANOPTION_H_
#define AMERICANOPTION_H_

//#include <iostream>
//#include <stdio.h>
//#include <stdio.h>
//#include <vector>
//#include <stdio.h>
//#include <math.h>
//#include <iostream>
//#include <stdlib.h>
//#include <time.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <math.h>
//#include <iostream>
//#include <stdlib.h>
//#include <cstring>

#include "../alglib/optimization.h"
#include "../alglib/stdafx.h"
#include "../alglib/ap.h"
#include "../alglib/optimization.h"

#include "EuroBewerter.h"

#include "RNG.h"
#include "MTRand.h"
#include "Hilfsmittel.h"

#define MAX_CALL 1  //bedeutet Put im Fall D=1
#define MIN_PUT 0   //bedeutet Put im Fall D=1
#define ITO 1
#define ITOrho 11
#define EULER 2
#define CIR 3
#define JDI 4
#define MTYP double    //Datentyp double oder float fuer die Speicherung angeben

//using namespace alglib;

namespace std {
class AmericanOption {
public:
	AmericanOption();
	virtual ~AmericanOption();

	int Iterations_Nummer;
	time_t Anfangszeit;
	int BFGS_Iterations;

	EuroBewerter EB;
	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	double* X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	double* sigma; //Volatility
	double T; //Gesamtzeit

	bool inExerciseRegion(double x, int n);
	bool inExerciseRegion(double* x, int n);


	double rho;
int Dates;
	int N; //time discretization
	int D;
	int M; // numer of training paths for BFGS
	int BFGS_Nesterov_Intervals;

	double abgeschnitten(double x, double f);

	double obj(double * alpha);
	double* obj_diff(double * alpha);
	double static const p = 3.; //Glättungsparameter

	void Wdiff_und_X_erstellen();
	void StochInt_erstellen();

	int* stoppzeiten;
	void stoppzeiten_erstellen();

	int K; //K=NN*5+1; //Anzahl der Basisfunktionen
	double dt;
	MTYP *** StochIntegrals;
	double** GrundIntegrals;
	double*** X;
	double*** Y;
	double*** WDiff;

	int* Exercise_Dates;
	int number_of_Exercise_Dates;
	int lauf;

	bool longstaffschwarz;
	bool bfgs;
	bool nesterov;
	int Threadanzahl;
	int PfadModell;

void testen();

int KpI;
int NpI;
int*** reihe;

	double f1D(double x,double y, int k, int n);
	double f2D(double* x, int k, int n, int d);
	double f(int k, double **x, int n, int d);
	double* f2D_all(double* x,double *y,int* reihe, int k, int n);
	double* f_all(int k, double **x,double** y,int** reihe, int n);
	double F(double x, int k, double border, bool hauf);
	double F2(double x, int k, double border, bool hauf);
	//MTRand MT; //MersenneTwister
	void Pfadgenerieren(double** X, double** wdiff, int start, double* S);
	void Pfadgenerieren(double** X,  int start, double* S, RNG* generator);
	void Pfadgenerieren(double** X, double** wdiff);
	void Daten();
	void neueExerciseDates(int n);
	double BoxMuller(double U1, double U2);

	double max(double d1, double d2);

	double payoff(double* x, int time);
	double payoff(double** x, int time);

	void objfs_aufrufen(double* x, double &func, double* grad);
	void objfs(double* x, double& func, double* grad);
	void BFGS(double* alpha);
	void CG(double* alpha);
	void BFGS_aufrufen();

	//Nesterov members
	int Nesterov_Iterations;
	double Nesterov_L;
	double* fx;
	double** x;
	double** y;
	double** z;
	double** Dfx;
	double* argminProblem(int lauf);
	double* TQx(int lauf);
	void Nesterov_aufrufen();
	void Nesterov(double* alpha, double L, int iterations);
	void objfsNesterov(double* alphas, double& func, double* gradient, int M);

	//Longstaff and Schwarz members
	void LSM_setting();
	void LongstaffSchwartz();
	void LSM_mittelwert(int threadnummer);
	double LSM_C_estimated(double* x, int ex);
	double LSM_phi(double* x, int j, int time);
	double **V;
	int Mphi;
	double** betas;
	double** B;
	double* BV;
	int LSlauf;
	bool mlsm;
	int LSM_K0;
	int LSM_K1;
	int LSM_K2;
	int LSM_K3;
	int LSM_K4;
	int LSM_Mtraining;
	int LSM_Mtesting;

};

} /* namespace std */
#endif /* AMERICANOPTION_H_ */
