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

AmericanOption::AmericanOption() {
	BFGS_Nesterov_Intervals = 5;
	Iterations_Nummer = 0;
	antithetics = verbose = loadAlphas = false;
	parallelTest = nesterov = extremTest = false;
	semiinf = speedup = andersenbroadie = false;
	bfgs = testing = longstaffschwarz = false;
	Daten();
	K = BFGS_Nesterov_Intervals * (K1 + K2 + K3 + K4 + K5); //Anzahl der Basisfunktionen
	dt = T / (double) (N - 1);
	Anfangszeit = time(NULL);
	zwischenwert_ohne_glaettung = DoubleFeld(Threadanzahl);
	//(double*)malloc(sizeof(double)*Threadanzahl);
}

AmericanOption::~AmericanOption() {

}

double AmericanOption::unif() {
	//return MT();
	return (double)(rand())/(double)(RAND_MAX);
}

double AmericanOption::BoxMuller(double U1, double U2) {
	double R = -2 * log(U1);
	double V = 2 * 3.1415926535 * U2;
	return sqrt(R) * cos(V);
}

double AmericanOption::nextGaussian() {
	//return qnorm(MT());//MersenneTwister benutzen
	//return qnorm( (double)(random())/(double)(RAND_MAX) );

	// Box-Muller algorithm
	double U1 = unif();
	double U2 = unif();
	if(U1==0)U1=0.00000001;
	double R = -2 * log(U1);
	double V = 2 * 3.1415926535 * U2;
	return sqrt(R) * cos(V);

	//----------------------- Die Polar-Methode von George Marsaglia
	//	double a1,a2,q=2;
	//	while(!((0 < q) && (q <= 1))){
	//		a1 = 2 * MT() - 1.;  // "Zufallszahl" liefert in [0,1)
	//		a2 = 2 * MT() - 1.;  //   gleichverteilte Werte
	//		q = a1 * a1 + a2 * a2;
	//	}
	//	double p = sqrt (-2 * log(q) / q);
	//	double z1 = a1 * p;
	////	double z2 = a2 * p;
	//	return z1;
}

double AmericanOption::max(double d1, double d2) {
	if (d1 < d2) return d2;
	else return d1;
}

double AmericanOption::payoff(double* x, int time) {
	//	return payoff(x,time,D);
	if (option == MIN_PUT)
		return max(Strike - Min(x, D), 0) * exp(-r * dt * (double) (time)); //Min Put
	if (option == MAX_CALL)
		return max(Max(x, D) - Strike, 0) * exp(-r * dt * (double) (time)); //Max Call
	printf("ERROR, option unknown!\n");
	return -1;
}

double* AmericanOption::payoffAbl(double* x, int time) {
	double* grad=DoubleFeld(D);
	if (option == MIN_PUT)
	{
		for(int d=0;d<D;++d)
			if(x[d]==Min(x,D) && x[d]<=Strike)grad[d]=-1.* exp(-r * dt * (double) (time));
		return grad;
	}
	if (option == MAX_CALL)
	{
		for(int d=0;d<D;++d)
			if(x[d]==Max(x,D) && x[d]>=Strike)grad[d]=1.* exp(-r * dt * (double) (time));
		return grad;
	}
	printf("ERROR, option unknown!\n");
	return NULL;
}

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
	Exercise_Dates = (int*) malloc(sizeof (int) *(n + 1));
	number_of_Exercise_Dates = n;
	if (verbose && !parallelTest)printf("Exercise dates: ");
	for (int e = 0; e < number_of_Exercise_Dates; ++e) {
		Exercise_Dates[e] = (int) ((double) (N - 1)*(double) (e) / (double) (number_of_Exercise_Dates - 1));
		if (verbose && !parallelTest)printf("%f, ", (double) Exercise_Dates[e] * dt);
	}
	if (verbose && !parallelTest)printf("\n");
}

void AmericanOption::Pfadgenerieren(double** X) {
	double** wdiff = DoubleFeld(N, D);
	for (int n = 0; n < N; ++n)
		for (int j = 0; j < D; ++j)
			wdiff[n][j] = sqrt(dt) * nextGaussian();
	Pfadgenerieren(X, wdiff);
	deleteDoubleFeld(wdiff,N,D);
}

void AmericanOption::Pfadgenerieren(double** X, int start, double * S) {
	double** wdiff = DoubleFeld(N, D);

	for (int n = 0; n < N; ++n)
		for (int j = 0; j < D; ++j)
			wdiff[n][j] = sqrt(dt) * nextGaussian();

	Pfadgenerieren(X, wdiff, start, S);
	deleteDoubleFeld(wdiff,N,D);
}

void AmericanOption::Pfadgenerieren(double** X, double** wdiff) {
	Pfadgenerieren(X, wdiff, 0, X0);
}

void AmericanOption::Pfadgenerieren(double** X,  int start, double* S, RNG* generator) {
	double** wdiff =DoubleFeld(N,D);
	for(int n=0;n<N;++n)
		for(int d=0;d<D;++d)
			wdiff[n][d]=sqrt(dt)*generator->nextGaussian();
	Pfadgenerieren(X, wdiff, 0, S);
	deleteDoubleFeld(wdiff,N,D);
}

void AmericanOption::Pfadgenerieren(double** X, double** wdiff, int start, double * S) {
	for (int d = 0; d < D; ++d)
		X[start][d] = S[d];

	for (int d = 0; d < D; ++d) {
		for (int n = start + 1; n < N; ++n) {
			if (PfadModell == ITO)
				X[n][d] = X[n - 1][d] * exp((((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt + sigma[d] * wdiff[n][d]));
			if (PfadModell == EULER)
				X[n][d] = X[n - 1][d] + (r - delta) * X[n - 1][d] * dt + sigma[d] * X[n - 1][d] * wdiff[n][d];
			if (PfadModell == CIR)
				X[n][d] = max(X[n - 1][d] + kappa * (theta - X[n - 1][d]) * dt + sigma[d] * sqrt(X[n - 1][d]) * wdiff[n][d], 0); //mean reversion
			if (PfadModell == ITOrho)
			{
				if(D==2){
//				     [,1]      [,2]
//				[1,]  1.0 0.0000000
//				[2,]  0.3 0.9539392
					double z[2];
					z[1]=wdiff[n][0];
					z[0]=0.3*wdiff[n][0]+0.9539392*wdiff[n][1];
					X[n][d] = X[n - 1][d] * exp((((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt + sigma[d] * z[d]));
				}
				if(D==3){
//				     [,1]      [,2]     [,3]
//				[1,]  1.0 0.0000000 0.000000
//				[2,]  0.3 0.9539392 0.000000
//				[3,]  0.3 0.2201398 0.928191
									double z[3];
									z[0]=wdiff[n][0];
									z[1]=0.3*wdiff[n][0]+0.9539392*wdiff[n][1];
									z[2]=0.3*wdiff[n][0]+0.2201398*wdiff[n][1]+0.928191*wdiff[n][2];
									X[n][d] = X[n - 1][d] * exp((((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt + sigma[d] * z[d]));
								}
			}
			//if (PfadModell == JDI)
			//	X[n][j] = X[n - 1][j] * exp(((r - delta) - 0.5 * sigma[j] * sigma[j]) * dt + sigma[j] * wdiff[n][j]) * exp(sprue[n][j]);
			//	//			X[n][j] = max( X[n - 1][j] + (r-delta) *X[n-1][j]*dt + sigma[j] *X[n-1][j]*wdiff[n][j] +X[n - 1][j] *sprue[n][j],0);
		}
		if (X[N - 1][0] <= 0)printf("Error0\n");
	}
}

double AmericanOption::newSprung() {
	return nextGaussian() * eta/*-0.5*eta*eta*/;
}

int AmericanOption::Poisson(double theta) {
	double p = exp(-theta);
	double F = p;
	double N = 0;
	double U = MT();
	while (U > F) {
		N++;
		p = p * theta / N;
		F += p;
	}
	return N;
}

double AmericanOption::EuropeanPut1D_discounted(double t, double T, double S, double strike) {
	if (t >= T)return max(strike - S, 0);
	double d1 = (log(S / strike) + (r + sigma[0] * sigma[0] / 2.)
			* (T - t)) / (sigma[0] * sqrt(T - t));
	double d2 = d1 - sigma[0] * sqrt(T - t);
	double P = CumulativeNormalDistribution(-d2) * strike * exp(-r * (T - t)) - CumulativeNormalDistribution(-d1) * S;
	return P * exp(-t * r);
}

double AmericanOption::EuropeanCall1D_discounted(double t, double T, double S, double strike) {
	if (t >= T)return max(S - strike, 0);
	double P = EuropeanPut1D(t, T, S, strike);
	double B = exp(-r * (T - t));
	double C = P + S - strike * B;
	return C * exp(-t * r);
}

double AmericanOption::EuropeanPut1D(double t, double T, double S, double strike) {
	if (t >= T)return max(strike - S, 0);
	double d1 = (log(S / strike) + (r + sigma[0] * sigma[0] / 2.)
			* (T - t)) / (sigma[0] * sqrt(T - t));
	double d2 = d1 - sigma[0] * sqrt(T - t);
	double P = CumulativeNormalDistribution(-d2) * strike * exp(-r * (T - t)) - CumulativeNormalDistribution(-d1) * S;
	return P;
}

double AmericanOption::EuropeanCall1D(double t, double T, double S, double strike) {
	if (t >= T)return max(S - strike, 0);
	double P = EuropeanPut1D(t, T, S, strike);
	double B = exp(-r * (T - t));
	double C = P + S - strike * B;
	return C;
}

double AmericanOption::EuropeanOption1D_discounted(double t, double T, double S, double strike) {
	if (option == MIN_PUT)return EuropeanPut1D_discounted(t, T, S, strike);
	if (option == MAX_CALL)return EuropeanCall1D_discounted(t, T, S, strike);
	printf("Error 378");
	return 0;
}

double AmericanOption::european_MaxCall_ND(double* x, double t, double T) {
	/*//if(k==0) return 1;
    //	double a=((double)(k)+0.5)/(double)K*180.;
    double d_minus[D];
    double d_plus[D];
    for (int d = 0; d < D; ++d) {
        d_minus[d] = (log(x[d] / Strike)+(r - delta - sigma[d] * sigma[d] / 2.)*(T - t)) / sigma[d] / sqrt(T - t + 0.000001);
        d_plus[d] = d_minus[d] + sigma[d] * sqrt(T - t);
    }

    double erg = 0;
    for (int l = 0; l < D; ++l) {
        double integralSumme = 0;
        double dz = 0.01;
        //	printf("%f, %f, %f\n", d_plus, d_minus, t);
        for (double z = -8; z < d_plus[l]; z += dz) {
            double df = exp(-0.5 * z * z);
            for (int l_Strich = 0; l_Strich < D; ++l_Strich)
                if (l_Strich != l)df *= CumulativeNormalDistribution(log(x[l] / x[l_Strich]) / sigma[l] / sqrt(T - t + 0.000001) - z + sigma[l] * sqrt(T - t));
            //		if(z==-3)printf("%f, %f\n",z, df);
            integralSumme += df*dz;
            //			if(df<0.0001 && rand()%1000==0)
            //				printf("sehr klein bei %f\n",z);
        }

        erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654) * integralSumme;
    }
    double prod = 1;
    for (int l = 0; l < D; ++l)
        prod *= (1 - CumulativeNormalDistribution(d_minus[l]));
    erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
    return exp(-r * t)*(erg);*/
	//if(k==0) return 1;
	//	double a=((double)(k)+0.5)/(double)K*180.;
	double d_minus[D];
	double d_plus[D];
	for (int d = 0; d < D; ++d) {
		d_minus[d] = (log(x[d] / Strike)+(r - delta - sigma[d] * sigma[d] / 2.)*(T - t)) / sigma[d] / sqrt(T - t + 0.000001);
		d_plus[d] = d_minus[d] + sigma[d] * sqrt(T - t);
	}

	double erg = 0;
	for (int l = 0; l < D; ++l) {
		double integralSumme = 0;
		double dz = 0.01;
		//	printf("%f, %f, %f\n", d_plus, d_minus, t);
		for (double z = -8; z < d_plus[l]; z += dz) {
			double df = exp(-0.5 * z * z);
			for (int l_Strich = 0; l_Strich < D; ++l_Strich)
				if (l_Strich != l)df *= CumulativeNormalDistribution(log(x[l] / x[l_Strich]) / sigma[l] / sqrt(T - t + 0.000001) - z + sigma[l] * sqrt(T - t));
			//		if(z==-3)printf("%f, %f\n",z, df);
			integralSumme += df*dz;
			//			if(df<0.0001 && rand()%1000==0)
			//				printf("sehr klein bei %f\n",z);
		}

		erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654) * integralSumme;
	}
	double prod = 1;
	for (int l = 0; l < D; ++l)
		prod *= (1 - CumulativeNormalDistribution(d_minus[l]));
	erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
	return exp(-r * t)*(erg);
}

double AmericanOption::europeanValue(double* x, double t, double T) {
	if (D >= 2) {
		if (option == MAX_CALL)return european_MaxCall_ND(x, t, T);
	}
	if (D == 1) {
		if (option == MAX_CALL)return EuropeanCall1D_discounted(t, T, x[0], Strike);
		if (option == MIN_PUT)return EuropeanPut1D_discounted(t, T, x[0], Strike);
	}
	printf("Error17\n");
	return 0;
}


//
////Basisfunktionen vom Typ 4
//	if(k>=K1+K2+K3 && k<K1+K2+K3+K4)
//	{
//
//		double faktor=1;
//				if(j==0 && x[0]<x[1])faktor=-1;
//				if(j==1 && x[0]>x[1])faktor=-1;
//				return F(fabs(y[1]-y[0]+0.001),k-K1-K2-K3-K4,t)*faktor;
//		//		int j_extrem=0;
//		//		if(option==MAX_CALL)j_extrem=argMax(x,J);
//		//		if(option==MIN_PUT)j_extrem=argMin(x,J);
//		//		if(j==j_extrem){
//		//			double summe=0;
//		//			for(int jj=0;jj<J;jj++)
//		//				summe+=y[jj];
//		//			return F(summe,k-K1-K2,t);    //M1
//		//
//		//		}
//		//		else
//		//			return 0;
//		//		int j_extrem_minimum=0;
//		//		if(option==MAX_CALL)j_extrem_minimum=argMin(x,J);
//		//		if(option==MIN_PUT) j_extrem_minimum=argMax(x,J);
//		//		if(j==j_extrem_minimum){
//		//			double summe=0;
//		//			for(int jj=0;jj<J;jj++)
//		//				summe+=y[jj];
//		//			double y=log( fabs(x[j]-x[(j+1)%J])/x[(j+1)%J]+0.001)/sqrt(T-t*dt+0.00001);
//		//			return F(y,k-K1-K2-K3,t);    //M1
//		//		}
//		//		else return 0;
//
//		//		if(j==0){
//		//			double y=log(x[0]/Strike)/sqrt(T-t*dt+0.00001); //bringt ein bisschen
//		//		return F(y,k-K1-K2-K3,t);
//		//		}else
//		//			return 0;
//
//		//		double y=log(x[1]/Strike)/sqrt(T-t*dt+0.00001);
//		//		return F(y,k-K1-K2-K3,t,j);
//
//		//		double summe=0;
//		//		for(int jj=0;jj<J;jj++)
//		//			summe+=fabs(y[jj]*pow(-1,jj));
//		//		return F(summe,k-K1-K2,t,j);    //M1
//
//		//		double produkt=0;
//		//		for(int jj=0;jj<J;jj++)
//		//			produkt*=y[jj];
//		//		return F(produkt,k-K1-K2-K3,t,j);    //M1
////		if(J==1)return F2(x[0],k-K1-K2-K3,t,j); /// Soeren Christensen Funktionen
////		else return F3(x,k-K1-K2-K3,t,j);
//		//		//		double a=((double)(k-K1-K2-K3)+0.5)/(double)K*180.;
//		//		////		printf("%f\n",a);
//		//		//		return exp(-r*(T-dt*t))*																			//Vorfaktor
//		//		//				glockenfunktion(-(log(x[0]/a)+(r-sigma*sigma/2.)*(T-t*dt))/sigma/sqrt(T-t*dt+0.000001))* 	//Aeussere Ableitung
//		//		//				-a/x[0]/sigma/sqrt(T-t*dt+0.000001)* 														//innere Ableitung
//		//		//				1./a; 																						//zweite innere Ableitung
//	}
//
//	//Basisfunktionen vom Typ 5
//	if(k>=K1+K2+K3+K4){
//
////
////		int j_extrem=0;
////		if(option==MAX_CALL)j_extrem=argMax(x,J);
////		if(option==MIN_PUT) j_extrem=argMin(x,J);
////		if(j==j_extrem){
////
////		}
////		double summe=0;
////		for(int jj=0;jj<J;jj++)
////			summe+=y[jj];
////		return F(abs(),k-K1-K2,t);    //M1
//		//			return exp(-pow((y[0]-y[1]),2)+k*sqrt(0.001+fabs(y[0]+y[1])))*(-2*y[0]+2*y[1]+(k-K1-K2-K3-K4)/(2*sqrt(0.001+fabs(y[0]+y[1]))));
//
//		//		int j_extrem=0;
//		//						if(option==MAX_CALL)j_extrem=argMax(x,J);
//		//						if(option==MIN_PUT)j_extrem=argMin(x,J);
//		//						if(j!=j_extrem){
//		//							double summe=0;
//		//									for(int jj=0;jj<J;jj++)
//		//										summe+=y[jj];
//		//									return F(summe,k-K1-K2,t);    //M1
//		//
//		//								}
//		//								else
//		//									return 0;
//
//		//			if(j==1){
//		//				double y=log(x[1]/Strike)/sqrt(T-t*dt+0.00001); // bringt ein Bisschen
//		//				return F(y,k-K1-K2-K3-K4,t);
//		//			}
//		//		else return 0;
//		//		double summe=0;
//		//		for(int jj=0;jj<J;jj++)
//		//			summe+=y[jj];
//		//		return F(abs(y[j]-summe)*sqrt(T-t*dt+0.00001),k-K1-K2-K3-K4,t);
//	}
//
//
//											//				if(k>=K1+K2+K3){
//											//					double ET=EuropeanPut1D_discounted(0,T,S0 ,Strike/2.+Strike*(double)k/(double)K);
//											//					double tt=dt*(double)n;
//											//					double E =EuropeanPut1D_discounted(   tt,T,X[m][n-1][0],Strike/2.+Strike*(double)k/(double)K   );
//											//					StochIntegrals[k][m][n]=E-ET;
//											//				}else{



//double AmericanOption::F2(double x, int k, int t, int j){  // Soeren Christensen Funktionen
//	if(k==0) return 1;
//	double border=0.8; //0.8
//	double y=log(x/Strike)/sqrt(T-t*dt+0.00001);
//	double extra=(y>=border)?1:0;
//	double a=((double)(k)+0.5)/(double)K*180.;
//	//	double erg=exp(-r*(T-dt*t))*																			//Vorfaktor
//	//					glockenfunktion(-(log(x/a)+(r-sigma*sigma/2.)*(T-t*dt))/sigma/sqrt(T-t*dt+0.000001))* 	//Aeussere Ableitung
//	//					-a/x/sigma/sqrt(T-t*dt+0.000001)* 														//innere Ableitung
//	//					1./a;
//	double d1=(log(x/a)+(r-sigma*sigma/2.)*(T-t*dt))/sigma/sqrt(T-t*dt+0.000001);
//	double erg=CumulativeNormalDistribution(d1)-1;
//	return erg*(fabs(y)<border)+extra;
//}


//double AmericanOption::BelomestFunktion2(double* x,int k, int t, int j)
//{
//	//	if(delta!=0)
//	//	double Strike=max((r/delta)*Strike,Strike);
//	//K1 Basisfunktionen von Typ 1
//	if(k<K1){
//		int j_extrem=0;
//		if(option==MAX_CALL)j_extrem=argMax(x,D);
//		if(option==MIN_PUT) j_extrem=argMin(x,D);
//		if(j==j_extrem){
//			//			double y[J];
//			//			for(int jj=0;jj<J;++jj)
//			//			double y=(log(x[j]/Strike))/sqrt(T-t*dt+0.00001);
//			double y=(log(x[j]/Strike))/sqrt(T-t*dt+0.00001);
//			return F(y,k);
//			//return 0.1;
//		}
//		else
//			return 0;
//	}
//
//	//K2 Basisfunktionen von Typ 2
//	if(k>=K1 && k<K1+K2){
//		double y=log(x[j]/Strike)/sqrt(T-t*dt+0.00001);
//		return F(y,k-K1);
//	}
//
//	//M3 Basisfunktionen von Typ 3
//	double y[D];
//	for(int jj=0;jj<D;++jj)
//		y[jj]=log(x[jj]/Strike)/sqrt(T-t*dt+0.00001);
//	if(k>=K2+K1 && k<K1+K2+K3)
//	{
//		int j_1=argMax(x,D);
//		int j_3=argDritter(x,D);
//		int j_2=argZweiter(x,D);
//		double summe=0;
//		if(j==j_1 || j==j_2 || j==j_3)
//		{
//			summe=y[j_1]+y[j_2]+y[j_3];
//			//					for(int jj=0;jj<J;jj++)
//			//				summe+=y[jj];
//			return F(summe,k-K1-K2);    //M1
//		}else
//			return 0;
//		//		return F(y[J-j-1],k-K1-K2,t,j);     //M2
//	}
//
//	//Basisfunktionen vom Typ 4
//	if(k>=K1+K2+K3 && k<K1+K2+K3+K4)
//	{
//		//		if(j==0 && x[0]<x[1])faktor=-1;
//		//		if(j==1 && x[0]>x[1])faktor=-1;
//		//		return F(fabs(y[1]-y[0]+0.001),k-K1-K2-K3-K4,t)*faktor;
//		//		int r[J];
//		//		for(int jj=0;jj<J;++jj)r[jj]=j;
//		//		BubbleSort(x,r,J);
//		int j_1=argMax(x,D);
//		//		int j_2=argMin(x,J);
//		int j_2=argZweiter(x,D);
//		//		int j_1=r[0];
//		//		int j_2=r[J-1];
//		double faktor=0;
//		if(j==j_1)faktor=1;
//		if(j==j_2)faktor=-1;
//		if(faktor==0)return 0;
//		return F(fabs(y[j_1]-y[j_2])+0.001,k-K1-K2-K3)*faktor;
//	}
//
//	//Basisfunktionen vom Typ 5
//	if(k>=K1+K2+K3+K4){
//		//		if(j==0 && x[0]<x[1])faktor=-1;
//		//		if(j==1 && x[0]>x[1])faktor=-1;
//		//		return F(fabs(y[1]-y[0]+0.001),k-K1-K2-K3-K4,t)*faktor;
//		int j_1=argMax(x,D);
//		//		int j_2=argMin(x,J);
//		double faktor=0;
//		if(j==j_1)faktor=1;else
//			faktor=-1;
//		double summe=0;
//		for(int jj=0;jj<D;jj++)
//			summe+=y[jj];
//		//		if(j==j_2)faktor=-1;
//		return F(fabs(y[j_1]-summe)+0.01,k-K1-K2-K3-K4)*faktor;
//	}
//
//	printf("ERROR134");
//	return 0;
//}


//neu ende
//	for (int m = 0; m < M; ++m) { //Summieren  ueber alle replications
//		double d[number_of_Exercise_Dates];
//		for(int ex=0;ex<number_of_Exercise_Dates;++ex)
//		{
//			d[ex]=payoff(Exercise_Dates[ex], X[m][Exercise_Dates[ex]]);
//			for (int k = 0; k < K; ++k)
//				d[ex] -= x[k] * StochIntegrals[k][m][Exercise_Dates[ex]];
//		}
//		double s=0;
//		for(int ex=0;ex<number_of_Exercise_Dates;++ex)
//			s+=exp((p * d[ex]));
//
//		for (int k = 0; k < K; ++k) {
//			gradient[k][m]=0;
//			double S =0;
//			for(int ex=0;ex<number_of_Exercise_Dates;++ex)
//				S+=StochIntegrals[k][m][Exercise_Dates[ex]] * exp(p * d[ex]);
//			gradient[k][m] = -S / s;
//		}
//
//		sup_glatt[m]= log(s) / p;
//		mean_glatt +=sup_glatt[m] /(double)(M);
//		mean_unglatt += Max(d,number_of_Exercise_Dates)/ (double)(M);
//	}
//if(verbose)printf("Zwischenwert Ohne glaettung: %f\n",mean_unglatt/(double)(M));
//zwischenwert_ohne_glaettung=mean_unglatt;


/*if(n*dt>T/2.){
if(k-K1-K2-K3-K4>=0 && k-K1-K2-K3-K4<K1){
        int j_extrem=0;
        if(D>1){
                if(option==MAX_CALL)j_extrem=reihe[0];
                if(option==MIN_PUT) j_extrem=reihe[D-1];
        }
        if(j==j_extrem)
                W=F(Xy[j],k-K1-K2-K3-K4,0.8,true);
}

if(k-K1-K2-K3-K4>=K1 && k-K1-K2-K3-K4<K1+K2){
        W=F(Xy[j],k-K1-K1-K2-K3-K4,0.8,false);
}

if(k-K1-K2-K3-K4>=K2+K1 && k-K1-K2-K3-K4<+K1+K2+K3){
        bool nehmen=false;
        if(option==MAX_CALL)
                for(int jjj=0; jjj<D/2+1;++jjj)
                        if(j==reihe[jjj])nehmen=true;
        if(option==MIN_PUT)
                for(int jjj=D-1; jjj>=0;--jjj)
                        if(j==reihe[jjj])nehmen=true;
        if(nehmen)
                W=F(s,k-K1-K2-K1-K2-K3-K4,1.6,false);
}

if(k-K1-K2-K3-K4>=K3+K2+K2 ){
        int j_1,j_2;
        if(option==MAX_CALL){
                j_1=reihe[0];
                j_2=reihe[1];
        }
        if(option==MIN_PUT){
                j_1=reihe[D-2];
                j_2=reihe[D-1];
        }
        double faktor=0;
        if(j==j_1)faktor=1;
        if(j==j_2)faktor=-1;

        W=F2(Xy[j_1]-Xy[j_2],k-K1-K2-K3-K4-K1-K2-K3,0.9,true)*faktor;
}
}*/



//W=F(pow(Xy[0]*Xy[1],2),k-K1-K2-K3-K4,2.8,false)*pow(Xy[1-j],2)*Xy[j];
//					if(j==0)
//						W=F(fabs(Xy[0]-Xy[1])+Xy[0]+Xy[1],k-K1-K2-K3-K4,2.8,false)*(1.+1.*(Xy[0]>Xy[1])-1.*(Xy[1]>Xy[0]));
//
//					if(j==1)
//						W=F(fabs(Xy[0]-Xy[1])+Xy[0]+Xy[1],k-K1-K2-K3-K4,2.8,false)*(1.+1.*(Xy[1]>Xy[0])-1.*(Xy[0]>Xy[1]));
//					if(j==0)W=F(Xy[0],k-K1-K2-K3-K4,0.8,true);

//
//					int j_extrem=0;
//
//									if(option==MAX_CALL)j_extrem=reihe[D-1];
//									if(option==MIN_PUT) j_extrem=reihe[0];
//
//								if(j==j_extrem)
//									W=F(Xy[j]*Xy[1-j],k-K1-K2-K3-K4,0.8,true)*Xy[1-j];
//											else W=F(Xy[j]*Xy[1-j],k-K1-K2-K3-K4,0.8,true);



//					if(k-K1-K2-K3-K4==0)W=1;
//					if(k-K1-K2-K3-K4==1)if(j==0)W  =1;
//					if(k-K1-K2-K3-K4==2)if(j==1)W  =1;
//					if(k-K1-K2-K3-K4==3)if(j==0)W  =X[n-1][0]/sqrt(T-n*dt+0.00001);
//					if(k-K1-K2-K3-K4==4)if(j==1)W  =X[n-1][1]/sqrt(T-n*dt+0.00001);
//					if(k-K1-K2-K3-K4==5){if(j==1)W =X[n-1][0]/sqrt(T-n*dt+0.00001); else W=X[n-1][1]/sqrt(T-n*dt+0.00001);}

//					printf("Error 431\n");
//					int r_jetzt=(k-K1-K2-K3-K4)%D;
//					int p=      (k-K1-K2-K3-K4)/D;
//					if(j==reihe[r_jetzt])
//						W= F(Xy[reihe[0]]-Xy[j],p);

//					double prod=1;
//					for(int d=0;d<D;++d)
//						prod*=X[n-1][d];
//					W= F(prod,(k-K1-K2-K3-K4),2.2,false)*prod/X[n-1][j];
//
//					int j_extrem=0;
//										if(D>1){
//											if(option==MAX_CALL)j_extrem=reihe[0];
//											if(option==MIN_PUT) j_extrem=reihe[D-1];
//										}
//										if(j==j_extrem)

//					int j_extrem=0;
//					if(D>1){
//						if(option==MAX_CALL)j_extrem=reihe[0];
//						if(option==MIN_PUT) j_extrem=reihe[D-1];
//					}
//					if(j==j_extrem)
//						W=F(0.5*Xy[j]+0.25*(Xy[0]+Xy[1]),k-K1-K2-K3-K4,2.8,true)*(0.25+0.5);
//					else
//						W=F(0.5*Xy[j]+0.25*(Xy[0]+Xy[1]),k-K1-K2-K3-K4,2.8,true)*0.25;


