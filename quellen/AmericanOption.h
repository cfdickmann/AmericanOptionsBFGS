#ifndef AMERICANOPTION_H_
#define AMERICANOPTION_H_

#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
//#include "../src/optimization.h"

#include <time.h>
#include "RNG.h"
#include "MTRand.h"
//#include "../src/ap.h"
#include "math.h"
#include "Hilfsmittel.h"
//#include "../src/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "../src/optimization.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

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

	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	double* X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	double* sigma; //Volatility
	double T; //Gesamtzeit

	double lambdaJump; //for jumps
	double eta; //for jumps
	double kappa; // for mean reversion
	double theta; // for mean reversion

	double rho;
	int Testing_Dates;
	int Training_Dates;
	int N; //time discretization
	int D;
	int K1;
	int K2;
	int K3;
	int K4;
	int K5;
	int M; // numer of training paths for BFGS
	int BFGS_Nesterov_Intervals;

	int K; //K=NN*5+1; //Anzahl der Basisfunktionen
	double dt;
	MTYP *** StochIntegrals;
	double*** X;
	double*** WDiff;
	double*** Sprue;
	double* alpha;

	int* Exercise_Dates;
	int number_of_Exercise_Dates;
	int lauf;
	bool antithetics;
	bool verbose;
	bool testing;
	bool parallelTest;
	bool semiinf;
	bool extremTest;
	bool longstaffschwarz;
	bool andersenbroadie;
	bool speedup;
	bool bfgs;
	bool loadAlphas;
	bool nesterov;
	double* zwischenwert_ohne_glaettung;
	double zwischenwert_pene; //penelization term
	int Threadanzahl;
	int PfadModell; //Ito or euler or CIR
	double** gradient;
	double* sup_glatt;

	MTRand MT; //MersenneTwister
	void Pfadgenerieren(double** X, int start, double * S);
	void Pfadgenerieren(double** X, double** wdiff, int start, double* S);
	void Pfadgenerieren(double** X,  int start, double* S, RNG* generator);
	void Pfadgenerieren(double** X, double** wdiff);
	void Pfadgenerieren(double** X);
	void Daten();
	void neueExerciseDates(int n);
	double newSprung();
	double BoxMuller(double U1, double U2);
	void stuetzerwartung_ausrechnenThread(int k);

	double euklidMetrik(double* s1, double* s2);

	double max(double d1, double d2);
	double nextGaussian();
	int Poisson(double theta);
	double payoff(double* x, int time);
	double* payoffAbl(double* x, int time);

	double Upayoff(double* x, int time);
	double unif();
	double EuropeanPut1D_discounted(double t, double T, double S, double Strike);
	double EuropeanCall1D_discounted(double t, double T, double S, double Strike);
	double EuropeanPut1D(double t, double T, double S, double Strike);
	double EuropeanCall1D(double t, double T, double S, double Strike);
	double europeanValue(double* x, double t, double T);
	double EuropeanOption1D_discounted(double t, double T, double S, double Strike);

	//BFGS members
	void BFGS_StochIntgenerieren(int threadNummer);
	void BFGS_setting();
	double BFGS_testing(double number_of_replications);
	void BFGS_parallelTesting(double number_of_replications);
	void BFGS_extremeTesting(int l, double number_of_replications);
	void objfs_aufrufen(double* x, double &func, double* grad);
	void objfs(double* x, double& func, double* grad);
	void BFGS();
	void BFGS_mittelwert(int threadNummer);
	void StochInt(double** STi, double** X, double** WDiff, double** Sprue);
	//double* Basis(double x, int j, int n);
	//double BelomestFunktion(double* x,int k, int t, int j);
	double F(double x, int k, double border, bool hauf);
	//double F2(double x, int k, double border, bool hauf);
	//double F2(double x, int k, int t, int j);
	//double F3(double* x, int k, int t, int j);
	//double BelomestFunktion2(double* x,int k, int t, int j);
	//double BelomestFunktion3(int k, int m, int t, int j, int tre);
	//double f_Basisfunktionen(int k, int t, double *x, int j);
	double european_MaxCall_ND(double* x, double t, double T);

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
	void Nesterov();
	void objfsNesterov(double* alphas, double& func, double* gradient, int M);
	//void argMinfunction1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *ptr);
	//void TQfunction1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *ptr);

	//Longstaff and Schwarz members
	void LSM_setting();
	void LongstaffSchwartz();
	void LSM_mittelwert(int threadnummer);
	double LSM_C_estimated(double* x, int time);
	double LSM_phi(double* x, int j, int time);
	double **V;
	int Mphi;
	double** betas;
	double*** B;
	double** BV;
	int LSlauf;
	bool mlsm;
	int LSM_K0;
	int LSM_K1;
	int LSM_K2;
	int LSM_K3;
	int LSM_K4;
	int LSM_Mtraining;
	int LSM_Mtesting;

	//AndersenBroadie members
	void AndersenBroadie();
	void AndersenBroadieTest();
	double AndersenBroadieEinzel(double ** x, int nsubpaths);
	double AndersenBroadieEinzel(double ** x, int nsubpaths, int startzeit);
	double AndersenBroadieEinzelRand(double ** x, int nsubpaths, int startzeit);
	double AndersenBroadieMulti(double ** x, int nsubpaths);
	void AndersenBroadie5000(int threadnummer);
	double** betasAusDateiLaden();
	int testlaufactual;
	int AB_Level;
	int* AB_k;
	int* AB_n;
	double*** AB_mittelwerte;
	double*** AB_quadratsummen;

	//semi members
	void  stuetzpunkte_ausrichten();
	bool* stuetzstelle_active;
	double * actualkoeff;
	double koeff_testen(double* koeff);
	void koeff_testen_THREAD(int threadnummer);
	double *** koeff_testingpaths;
	int durchlaeufe;
	double ***** semi_inner_paths;
	double** semi_betas;
	int ** semi_betas_index;
	int * semi_betas_index_max;
	//        double* LP_mitR_Loesen();
	//        double* LP_mitALGLIB_Loesen();
	double* LP_mitGLPK_Loesen();
	void semi_testing();
	void semi_testThread(int threadnummer);
	int semi_testingpaths;
	void semi_mehrere_S0_testen();
	void semi_ergebnisse_ausgeben();
	void stuetzpunkte_ausgeben();
	void semi_inner_paths_erzeugen();
	void inner_paths_erzeugen_THREAD(int threadnummer);
	int durchlaufactual;
	void stuetzerwartungen_ausrechnen();
	void stuetzerwartung_ausrechnen();

	void stuetzpunkte_setzen(int n);

	void lp_ausgeben();
	int J;
	double** Matrix;
	double* RS;
	double* C;
	double** stuetzpunkte;
	double* stuetzerwartung;
	void StuetzErwartung(int t);
	double semi_Basisfunktionen(int zeit, int j, double* x);
	double* semi_BasisfunktionenAbl(int zeit, int j, double* x);
	double semi_Basisfunktionen1D(int zeit, int j, double* x);
	double* semi_Basisfunktionen1DAbl(int zeit, int j, double* x);
	double semi_Basisfunktionen2D(int zeit, int j, double* x);
	double* semi_Basisfunktionen2DAbl(int zeit, int j, double* x);
	double semi_BasisfunktionenHigherD(int zeit, int j, double* x);
	double* semi_BasisfunktionenHigherDAbl(int zeit, int j, double* x);
	void semi();
	//	void Test();
	int nactual;
	double linearCombinationOfBasis(int zeit, double* x);
	double* linearCombinationOfBasis_Abl(int zeit, double* x);
	double linearCombination(double* koeff, double* x);
	double semi_f(int n, double* x);
	double semi_f_Abl(int n, double* x, int d);
};

} /* namespace std */
#endif /* AMERICANOPTION_H_ */
