#include "AmericanOption.h"
using namespace std;

#include "Hilfsmittel.h"

int argMinlauf;
AmericanOption* zeiger2;

double* AmericanOption::argminProblem(int lauf) {
	double* erg;
	erg = (double*) malloc(sizeof(double) * K);
	for (int k = 0; k < K; ++k) {
		double S = 0;
		for (int l = 0; l <= lauf; ++l)
			S += ((double) l + 1.) / 2. * Dfx[l][k];
		erg[k] = -S / Nesterov_L;
	}
	return erg;
}

int TQlauf;

double* AmericanOption::TQx(int lauf) {
	double* erg;
	erg = (double*) malloc(sizeof(double) * K);
	for (int k = 0; k < K; ++k)
		erg[k] = x[lauf][k] - Dfx[lauf][k] / Nesterov_L;
	return erg;
}

void AmericanOption::Nesterov_aufrufen() {
	LongstaffSchwartz();
	Wdiff_und_X_erstellen();
	stoppzeiten_erstellen();
	StochInt_erstellen();

	double alpha[K];
	for (int k = 0; k < K; ++k)
		alpha[k] = 0.;

	M /= 10;
	Nesterov(alpha, 1000., 100);
	M *= 10;
	Nesterov(alpha, 1000., 100);

}

void AmericanOption::Nesterov(double* alpha, double L, int iterations) {
	Nesterov_Iterations = iterations;
	Nesterov_L = L;

	zeiger2 = this;

	fx = new double[Nesterov_Iterations + 100];
	x = DoubleFeld(Nesterov_Iterations + 100, K);
	y = DoubleFeld(Nesterov_Iterations + 100, K);
	z = DoubleFeld(Nesterov_Iterations + 100, K);
	Dfx = DoubleFeld(Nesterov_Iterations + 100, K);

	for (int k = 0; k < K; ++k)
		x[0][k] = alpha[k];

	printf("Nesterov\n");

	for (int lauf = 0; lauf < Nesterov_Iterations - 1; ++lauf) {
		printf("\n\nNesterov-Iteration: %d: \n", lauf);

		//Schritt 1
		TQlauf = argMinlauf = lauf;
		objfs_aufrufen(x[lauf], fx[lauf], Dfx[lauf]);

		printf("Df_%d =  ", lauf);
		ausgeben(Dfx[lauf], K);
		y[lauf] = TQx(lauf);

		//Schritt 3
		z[lauf] = argminProblem(lauf);
//		if (verbose) {
//			printf("z_%d =   ", lauf);
//			ausgeben(z[lauf], K);
//		}

		//Schritt 4
		for (int k = 0; k < K; ++k)
			x[lauf + 1][k] = 2. / ((double) lauf + 3.) * z[lauf][k]
					+ ((double) lauf + 1.) / ((double) lauf + 3.) * y[lauf][k];
//		if (verbose) {
//			printf("x_%d+1 = ", lauf);
//			ausgeben(x[lauf + 1], K);
//		}
		for (int k = 0; k < K; ++k)
			alpha[k] = y[lauf][k];
	}
}

//void AmericanOption::Nesterov()
//{
//	Nesterov_Iterations=4000;
//	Nesterov_L=3000;
//
//	zeiger2=this;
//
//	Daten();
//
//	fx=(double*)malloc(sizeof(double)*max(Nesterov_Iterations,705));
//	x=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		x[i]=(double*)malloc(sizeof(double)*K);
//	y=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		y[i]=(double*)malloc(sizeof(double)*K);
//	z=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		z[i]=(double*)malloc(sizeof(double)*K);
//	Dfx=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		Dfx[i]=(double*)malloc(sizeof(double)*K);
//
//	for(int k=0;k<K;++k)
//		x[0][k]=0.000;
//
//	if(verbose)printf("Alphas aus Datei laden\n");
//	if(loadAlphas){
//		double* alpha=alphasLaden(K);
//		for (int k = 0; k < K; ++k) x[k][0]=alpha[k];
//	}
//
//	int M_store=M;
//	if(verbose)printf("Nesterov\n");
//
//	for(int lauf=0;(lauf<Nesterov_Iterations-1)||(speedup);++lauf)
//	{
//		printf("\n\nNesterov-Iteration: %d: \n",lauf);
//
//		TQlauf=argMinlauf=lauf;
//		objfs_aufrufen(x[lauf],fx[lauf],Dfx[lauf]);
//
//	{printf("Df_%d =  ",lauf);ausgeben(Dfx[lauf],K);}
//
//		//Schritt 2
//
//		y[lauf]=TQx(lauf);
//
//
//		//Schritt 3
//		z[lauf]=argminProblem(lauf);
//		if(verbose){printf("z_%d =   ",lauf);ausgeben(z[lauf],K);}
//
//		//Schritt 4
//		for(int k=0;k<K;++k)
//			x[lauf+1][k]=2./((double)lauf+3.)*z[lauf][k]+((double)lauf+1.)/((double)lauf+3.)*y[lauf][k];
//		if(verbose){printf("x_%d+1 = ",lauf);ausgeben(x[lauf+1],K);}
//                for(int k=0;k<K;++k)
//                    alpha[k]=y[lauf][k];
//	}
//        alphasSchreiben(alpha,K);
//}
