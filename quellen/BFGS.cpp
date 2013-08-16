#include <ctype.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>

#include "../alglib/optimization.h"
#include "../alglib/linalg.h"
#include "../alglib/ap.h"

#include "AmericanOption.h"

using namespace std;
using namespace alglib;

static AmericanOption* zeiger;

double* zx;

void* DELEGATE_BFGS_mittelwert(void* data) {
	zeiger->BFGS_mittelwert(((int*) data)[0]);
	pthread_exit(NULL);
	return NULL;
}

double AmericanOption::obj(double * alpha) {
	double ergE = 0;
	double erg1 = 0;
	double ergQ1 = 0;

	double p = 2.; //Glättungsparameter
	for (int m = 0; m < M / 2; ++m) { //Summieren  ueber alle replications
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][Exercise_Dates[ex]];
		}

		double s = 0;
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
			s += exp((p * d[ex]));

		double u = log(s) / p;
		erg1 += u / (double) (M / 2);
		ergQ1 += u * u / (double) (M / 2);

		ergE += Max(d, number_of_Exercise_Dates) / (double) M;
	}

	double erg2 = 0;
	double ergQ2 = 0;
//	for (int m = 0; m < M/2; ++m) {
//		double lincomb = 0;
////		printf("stopp %d\n",stoppzeiten[m]);
//		for (int k = 0; k < K; ++k)
//			lincomb += alpha[k] * StochIntegrals[k][m][stoppzeiten[m]];
//		double u = payoff(X[m], N - 1) - lincomb;
//		erg2 += u / (double) (M/2);
//		ergQ2 += u * u / (double) (M/2);
//	}
	double var2 = (ergQ2 - erg2 * erg2);

	for (int m = M / 2; m < M; ++m) { //Summieren  ueber alle replications
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][Exercise_Dates[ex]];
		}
		ergE += Max(d, number_of_Exercise_Dates) / (double) M;
	}

	double V = 1. * erg1 + 0. * var2;

	static int zaehler = 0;
	zaehler++;
	if (zaehler % zeiger->K == 0) {
		for (int k = 0; k < K; ++k)
			printf("%.2lf,", alpha[k]);
		printf("\nV=%f\n", V);
		printf("dual(testing): %f\n", ergE);
	}

	if (zaehler % zeiger->K == 0) {
		printf("schreiben\n");
		fstream f;
		f.open("mart.txt", ios::out);
		for (int n = 0; n < N; ++n) {
			for (double x = 10; x < 300; ++x) {
				double** xx = DoubleFeld(N, D);
				xx[n][0] = x;
				double func = 0;
				for (int k = 0; k < K; ++k)
					func += alpha[k] * this->f(k, xx, n, 0) / x;
				f << func << endl;
			}
			f << endl;
		}
		f.close();
		system("gnuplot plot.gnuplot");
	}
	return V;
}

void AmericanOption::objfs(double* x, double &func, double* grad) {
	double mean_glatt = 0;
//	double mean_unglatt=0;
	zwischenwert_pene = 0;
	double GRAD[K];
	for (int k = 0; k < K; ++k)
		GRAD[k] = 0;

//neu
	zx = x;
	zeiger = this;
	pthread_t threads[Threadanzahl];
// Alle Threads starten
	for (int i = 0; i < Threadanzahl; i++) {
		int* z = (int*) malloc(sizeof(int));
		z[0] = i;
		pthread_create(&threads[i], NULL, DELEGATE_BFGS_mittelwert, (void*) z);
		//                       Argument für thread_function ---^
	}

// Hier wird gewartet bis alle Threads beendet sind.
	for (int i = 0; i < Threadanzahl; i++) {
		pthread_join(threads[i], NULL);
	}

	mean_glatt = 0;
	for (int m = 0; m < M; ++m)
		mean_glatt += sup_glatt[m] / (double) (M);

//geglaettetes Ergebnis
	for (int k = 0; k < K; ++k) {
		GRAD[k] = 0;
		for (int m = 0; m < M; ++m)
			GRAD[k] += gradient[k][m] / (double) (M);
	}

	if (true) {
		double h = 0.00001;
		zx[4] += h;
		pthread_t threads[Threadanzahl];
		// Alle Threads starten
		for (int i = 0; i < Threadanzahl; i++) {
			int* z = (int*) malloc(sizeof(int));
			z[0] = i;
			pthread_create(&threads[i], NULL, DELEGATE_BFGS_mittelwert,
					(void*) z);
			//                       Argument für thread_function ---^
		}

		// Hier wird gewartet bis alle Threads beendet sind.
		for (int i = 0; i < Threadanzahl; i++) {
			pthread_join(threads[i], NULL);
		}

		double fh = 0;
		for (int m = 0; m < M; ++m)
			fh += sup_glatt[m] / (double) (M);
		zx[4] -= h;

		double gr = (fh - mean_glatt) / h;

		printf("\nunterschied:%f\n", GRAD[4] - gr);
	}

	double erg = mean_glatt;

//penelization...
//...fuer f
	double lambda = 0.;
	if (lambda > 0) {
		double W = 0;
		for (int m = 0; m < M; ++m)
			W += pow(sup_glatt[m] - mean_glatt, 2);
		erg += lambda / sqrt(M) * sqrt(W);
		zwischenwert_pene = erg - mean_glatt;

		//...fuer den gradient
		double mean_gradient[K];
		for (int k = 0; k < K; ++k) {
			mean_gradient[k] = 0;
			for (int m = 0; m < M; ++m)
				mean_gradient[k] += gradient[k][m] / (double) M;
		}

		for (int k = 0; k < K; ++k)
			for (int m = 0; m < M; ++m)
				GRAD[k] += lambda / sqrt(M) * (sup_glatt[m] - mean_glatt)
						* (gradient[k][m] - mean_gradient[k]) / sqrt(W);
	}
//arcus tangens
	bool atanAnwenden = false;
	if (atanAnwenden) {
		double q = 10.;
		for (int k = 0; k < K; ++k)
			GRAD[k] *= 1. / (1. + erg * erg / q / q) / q;
		erg = atan(erg / q);
	}

//Werte zurueckgeben
	func = erg;
	for (int k = 0; k < K; k++)
		grad[k] = GRAD[k]; // GRAD[i] wurden durch objfs gesetzt
}

//void AmericanOption::BFGS_mittelwert(int threadNummer) {
//	int mAnfang = M * threadNummer / Threadanzahl;
//	int mEnde = M * (threadNummer + 1) / Threadanzahl;
//	double p = 3; //  use 2 or 3
//	zwischenwert_ohne_glaettung[threadNummer] = 0;
//	for (int m = mAnfang; m < mEnde; ++m) { //Summieren  ueber alle replications
//		double d[number_of_Exercise_Dates];
//		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
//			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
//			for (int k = 0; k < K; ++k)
//				d[ex] -= zx[k] * StochIntegrals[k][m][ex];
//		}
//
//		double s = 0;
//		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
//			s += exp((p * d[ex]));
//
//		for (int k = 0; k < K; ++k) {
//			gradient[k][m] = 0;
//			double S = 0;
//			for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
//				S += StochIntegrals[k][m][ex] * exp(p * d[ex]);
//			gradient[k][m] = -S / s;
//		}
//
//		sup_glatt[m] = log(s) / p;
//		//			mean_glatt +=sup_glatt[m] /(double)(M);
//		zwischenwert_ohne_glaettung[threadNummer] += Max(d,
//				number_of_Exercise_Dates) / (double) (-mAnfang + mEnde);
//	}
//}

void AmericanOption::BFGS_mittelwert(int threadNummer) {
	int mAnfang = M * threadNummer / Threadanzahl;
	int mEnde = M * (threadNummer + 1) / Threadanzahl;
	double p = 3; //  use 2 or 3
	zwischenwert_ohne_glaettung[threadNummer] = 0;
	for (int m = mAnfang; m < mEnde; ++m) { //Summieren  ueber alle replications
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= zx[k] * StochIntegrals[k][m][ex];
		}

		double s = 0;
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
			s += exp((p * d[ex]));

		for (int k = 0; k < K; ++k) {
			gradient[k][m] = 0;
			double S = 0;
			for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
				S += StochIntegrals[k][m][ex] * exp(p * d[ex]);
			gradient[k][m] = -S / s;
		}

		sup_glatt[m] = log(s) / p;
		//			mean_glatt +=sup_glatt[m] /(double)(M);
		zwischenwert_ohne_glaettung[threadNummer] += Max(d,
				number_of_Exercise_Dates) / (double) (-mAnfang + mEnde);
	}
}

//void AmericanOption::objfs_aufrufen(double* x, double &func, double* grad) {
//	if (!speedup)
//		Iterations_Nummer++;
//	zeiger = this;
//	objfs(x, func, grad);
//
//	printf("\n\n");
//	if (nesterov)
//		printf("Nesterov ");
//	else
//		printf("BFGS ");
//	printf("Iteration: %d:", Iterations_Nummer);
//	if (speedup)
//		printf("(speedup) ");
//	time_t t = time(NULL);
//	printf("             Bisherige Rechenzeit: %ld Sekunden", t - Anfangszeit);
//	printf("\n");
//	printf("f(alpha) =  %.10lf\n", func);
//	printf("alpha = ");
//	ausgeben(x, K);
//	double Zw = 0;
//	for (int z = 0; z < Threadanzahl; ++z)
//		Zw += zwischenwert_ohne_glaettung[z] / (double) Threadanzahl;
//	printf("Zwischenwert Ohne glaettung: %.4lf\n", Zw);
//	printf("penelization term: %f\n", zeiger->zwischenwert_pene);
//
//	if (!speedup && Iterations_Nummer % 50 == 48) {
//		alphasSchreiben(x, K);
//		printf("------gespeichert------\n");
//	}
//
//	if (BFGS_Iterations <= Iterations_Nummer && ((testing) || parallelTest)
//			&& bfgs) //Abbrechen
//			{
//		for (int k = 0; k < K; ++k)
//			grad[k] = 0;
//	}
//}

void AmericanOption::stoppzeiten_erstellen() {
	stoppzeiten = IntFeld(M);
	double erg = 0;
	for (int m = 0; m < M; ++m) {
		if (m % 10 == 0) {
			printf("Stoppzeiten erstellen ... %d %% \r", m * 100 / M);
			cout.flush();
		}

		int stopp = number_of_Exercise_Dates - 1;
		for (int n = 1; n < number_of_Exercise_Dates - 1; ++n)
			if (payoff(X[m], Exercise_Dates[n])
					> LSM_C_estimated(X[m][Exercise_Dates[n]], n)) { //TODO Achtung
				stopp = n;
				break;
			}
		erg += payoff(X[m], Exercise_Dates[stopp]);
		stoppzeiten[m] = Exercise_Dates[stopp];
	}
	printf("Stoppzeiten erstellen ... fertig\n");
	printf("stoppzeiten durchschnitt %f \n", erg / (double) M);

}

//ALT
//static void Delegate_static_objfs(const real_1d_array &x, double &func,
//		real_1d_array &grad, void *ptr) {
//	int K = zeiger->K;
//	double* xx = (double*) malloc(sizeof(double) * K);
//	double* gradient = (double*) malloc(sizeof(double) * K);
//	for (int k = 0; k < K; ++k) {
//		gradient[k] = grad[k];
//		xx[k] = x[k];
//	}
//	zeiger->objfs_aufrufen(xx, func, gradient);
//
////	zeiger->BFGS_lauf++;
////	if(zeiger->BFGS_lauf<zeiger->BFGS_MAX_lauf)
//	for (int k = 0; k < K; ++k)
//		grad[k] = gradient[k];
////		elsefor(int k=0;k<K;++k)grad[k]=0;
////x[k]=xx[k];
//
//	bool numericDifferentiation = false;
//	if (numericDifferentiation) {
//		double f = func;
//		int K = zeiger->K;
//		double gr[K];
//		for (int k = 0; k < K; ++k) {
//			double h = 0.00000001;
//			//				real_1d_array b;
//			//				b.setlength(K);
//			double b[K];
//			for (int kk = 0; kk < zeiger->K; ++kk)
//				b[kk] = x[kk];
//			b[k] += h;
//			zeiger->objfs((double*) b, func, (double*) gradient);
//			double a = func;
//			gr[k] = (a - f) / h;
//		}
//
//		zeiger->objfs((double*) xx, func, (double*) gradient);
//		for (int k = 0; k < K; ++k) {
//			printf("unterschied: % f", grad[k] - gr[k]);
//			grad[k] = gr[k];
//		}
//	}
//}

//Mit numerischer Differentiation
static void Delegate_static_objfs(const real_1d_array &x, double &func,
		real_1d_array &grad, void *ptr) {

	int K = zeiger->K;
	double alpha[K];
	for (int k = 0; k < K; ++k)
		alpha[k] = x[k];

	func = zeiger->obj(alpha);

	for (int k = 0; k < K; ++k) {
		double h = 0.001;
		alpha[k] += h;
		double fh = zeiger->obj(alpha);
		grad[k] = (fh - func) / h;
		alpha[k] -= h;
	}
}

double AmericanOption::F(double x, int k, double border, bool hauf) {
	if (k == 0)
		return 1;
	double extra = (x >= border && hauf) ? 1 : 0;
	if (k % 2 == 0)
		return cos((double) k / 2. * x) * (fabs(x) < border) + extra;
	else
		return sin(((double) k + 1.) * x / 2.) * (fabs(x) < border) + extra;
}

//void AmericanOption::BFGS_StochIntgenerieren(int threadNummer) {
//	int mAnfang = M * threadNummer / Threadanzahl;
//	int mEnde = M * (threadNummer + 1) / Threadanzahl;
//	if (verbose)
//		printf("%d bis %d\n", mAnfang, mEnde);
//
//	if (verbose)
//		printf("Pfade erzeugen\n");
//	for (int m = mAnfang; m < mEnde; ++m) {
//		Pfadgenerieren(X[m], WDiff[m]);
//	}
//
//	if (verbose)
//		printf("StochInt generieren\n");
//	double** STi = DoubleFeld(K, N);
//	for (int m = mAnfang; m < mEnde; ++m) {
//		for (int k = 0; k < K; ++k)
//			STi[k][0] = 0;
//		StochInt(STi, X[m], WDiff[m], Sprue[m]);
//		if (m % 30 == 0)
//			printf("StochInt generieren: %d Prozent fertig\r",
//					(m - mAnfang) * 100 / (mEnde - mAnfang));
//		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
//			for (int k = 0; k < K; ++k)
//				StochIntegrals[k][m][ex] = STi[k][Exercise_Dates[ex]]; // Achtung
//			//		if(m%100==0)printf("%.0lf\%",m/(double)(mEnde-mAnfang));
//	}
//}
//
//void* DLEGATE_BFGS_StochIntegrieren(void* data) {
//	zeiger->BFGS_StochIntgenerieren(((int*) data)[0]);
//	pthread_exit(NULL);
//	return NULL;
//}

void AmericanOption::StochInt_erstellen() {
	StochIntegrals = DoubleFeld(K, M, N);
	for (int k = 0; k < K; ++k)
		for (int m = 0; m < M; ++m)
			StochIntegrals[k][m][0] = 0;

	for (int m = 0; m < M; ++m) {
		if (m % 10 == 0) {
			printf("StochInt erstellen ... %d %% \r", m * 100 / M);
			cout.flush();
		}
		for (int k = 0; k < K; ++k) {
			for (int n = 1; n < N; ++n) {
				double summe = 0;
				for (int d = 0; d < D; ++d)
					summe += f(k, X[m], n - 1, d) * WDiff[m][n][d];
				//	* (stoppzeiten[m] >= n); //TODO Achtung
				StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + summe;
			}
		}
	}
	printf("StochInt erstellen ... fertig\n");
}

//void AmericanOption::BFGS_setting() {
//
//	sup_glatt = (double*) malloc(sizeof(double) * M);
//	gradient = DoubleFeld(K, M);
//
//	if (verbose)
//		printf("Felder allocieren\n");
//	alpha = (double*) malloc(sizeof(double) * (K));
//
//	if (verbose)
//		printf("Zufallszahlen allocieren\n");
//	WDiff = DoubleFeld(M, N, D);
//
//	if (verbose)
//		printf("Spruenge allocieren\n");
//	Sprue = DoubleFeld(M, N, D);
//
//	if (verbose)
//		printf("Pfade allocieren\n");
//	X = DoubleFeld(M, N, D);
//
//	if (verbose)
//		printf("StochInt allocieren\n");
//	StochIntegrals = (MTYP***) malloc(sizeof(MTYP**) * K);
//	for (int k = 0; k < K; ++k) {
//		StochIntegrals[k] = (MTYP**) malloc(sizeof(MTYP*) * M);
//		for (int m = 0; m < M; ++m)
//			StochIntegrals[k][m] = (MTYP*) malloc(
//					sizeof(MTYP) * number_of_Exercise_Dates);
//	}
//
//	srand(time(NULL));
//	MT.seed(time(NULL));
//	if (verbose)
//		printf("Zufallszahlen erstellen\n");
//	for (int m = 0; m < M; ++m)
//		for (int n = 0; n < N; ++n) {
//			for (int j = 0; j < D; ++j) {
//				if (antithetics && m % 2 == 1)
//					WDiff[m][n][j] = WDiff[m - 1][n][j] * -1.;
//				else
//					WDiff[m][n][j] = sqrt(dt) * nextGaussian();
//				int NumberOfJumps = Poisson(lambdaJump * dt);
//				Sprue[m][n][j] = 0;
//				for (int jump = 0; jump < NumberOfJumps; ++jump)
//					Sprue[m][n][j] += newSprung();
//			}
//		}
//
//	zeiger = this;
//	pthread_t threads[Threadanzahl];
//// Alle Threads starten
//	int* z = new int[Threadanzahl];
//	for (int i = 0; i < Threadanzahl; i++) {
//		z[i] = i;
//		pthread_create(&threads[i], NULL, DLEGATE_BFGS_StochIntegrieren,
//				(void*) &(z[i]));
//		//                       Argument für thread_function ---^
//	}
//
//// Hier wird gewartet bis alle Threads beendet sind.
//	for (int i = 0; i < Threadanzahl; i++) {
//		pthread_join(threads[i], NULL);
//	}
//
//	delete[] z;
//}

void AmericanOption::BFGS_parallelTesting(double number_of_replications) {
	Daten();
	int ergebnispipe[Threadanzahl][2];
	for (int z = 0; z < Threadanzahl; ++z)
		pipe(ergebnispipe[z]);

	int pid;
	for (int f = 0; f < Threadanzahl; ++f) {
		pid = fork();
		if (pid == 0) {
			close(ergebnispipe[f][0]);
			char string[20];
			srand(f);
			double e = BFGS_testing(50000. / (double) Threadanzahl);
			sprintf(string, "%f", e);
			write(ergebnispipe[f][1], string, (strlen(string) + 1));
			exit(0);
		}
	}
	double erg = 0;
	double ergebnisse[Threadanzahl];
	for (int f = 0; f < Threadanzahl; ++f) {
		close(ergebnispipe[f][1]);
		char readbuffer[2000];
		read(ergebnispipe[f][0], readbuffer, sizeof(readbuffer));
		ergebnisse[f] = atof(readbuffer);
		erg += ergebnisse[f] / (double) (Threadanzahl);
		printf("Ergebnis %d: %f\n", f, atof(readbuffer));
	}
	double s = 0;
	for (int f = 0; f < Threadanzahl; ++f)
		s += pow(ergebnisse[f] - erg, 2);
	printf("geschaetzte Varianz: %f\n", s / double(Threadanzahl - 1));
	printf("Gesamtergebnis: %f\n", erg);
	ErgebnisAnhaengen(erg);
}

double AmericanOption::BFGS_testing(double number_of_replications) {
//	neueExerciseDates(Testing_Dates); //Achtung fehlt
	MT.seed(time(NULL) + getpid());
	srand(time(NULL) + getpid());
	if (!nesterov && !bfgs) {
		alpha = alphasLaden(K);
		if (verbose && !parallelTest)
			printf("alphas geladen\n");
	} else if (verbose && !parallelTest)
		printf("alphas vorhanden");

//Alphas ausgeben
	if (verbose && !parallelTest) {
		printf("alpha=[");
		for (int k = 0; k < K; ++k)
			printf("%f, ", alpha[k]);
		printf("]\n");
	}

//Felder fuer je 1 replication
	double** stochintegrals = DoubleFeld(K, N);
	double** x = DoubleFeld(N, D);
	double** wdiff = DoubleFeld(N, D);
	double** sprue = DoubleFeld(N, D);

	double mean = 0;

	int Npfade = 15;
	double w[N * Npfade]; //fuer die ausgabe von pfaden

	for (int lauf = 1; lauf <= number_of_replications; lauf++) {
		for (int n = 0; n < N; ++n)
			for (int j = 0; j < D; ++j)
				wdiff[n][j] = sqrt(dt) * nextGaussian();

		Pfadgenerieren(x, wdiff);

		if (lauf < Npfade + 1) {
			for (int n = 0; n < N; ++n)
				w[n + (lauf - 1) * N] = x[n][0];
		}
		if (lauf == Npfade + 1) {
			werteSchreiben(w, N * Npfade, N);
		}

		for (int k = 0; k < K; ++k)
			stochintegrals[k][0] = 0;
		StochInt(stochintegrals, x, wdiff, sprue);

		//for(int n=0;n<N;++n)
		//	printf("%.2lf ",stochintegrals[299][n]);

		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(x[Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * stochintegrals[k][Exercise_Dates[ex]];
		}

		mean += Max(d, number_of_Exercise_Dates);

		if (lauf % 100 == 0)
			printf(
					"Zwischenergebnis (kumuliert): %f bei %d Pfaden (%d Prozent fertig)\r",
					mean / ((double) (lauf)), lauf,
					(int) floor(
							(double) (lauf) / number_of_replications * 100.));
	}

	if (!parallelTest) {
		printf("\nresult: %f \n", mean / number_of_replications);
		ErgebnisAnhaengen(mean / number_of_replications);
	}
	return mean / number_of_replications;
}

void AmericanOption::StochInt(double** STi, double** X, double** WDiff,
		double** Sprue) {
	for (int n = 1; n < N; ++n)
		for (int k = 0; k < K; ++k) {
			double summe = 0;
			for (int d = 0; d < D; ++d)
				summe += f(k, X, n - 1, d) * WDiff[n][d];
			STi[k][n] = STi[k][n - 1] + summe;
		}
}

void AmericanOption::Wdiff_und_X_erstellen() {
	RNG generator;
	//generator.setSeed(seed);
	WDiff = DoubleFeld(M, N, D);
	X = DoubleFeld(M, N, D);

	for (int m = 0; m < M; ++m) {
//		if (m % 10 == 0) {
//			printf("X erstellen ... %d %%            \r", m * 100 / M);
//			cout.flush();
//		}

		{
			for (int n = 0; n < N; ++n)
				for (int d = 0; d < D; ++d)
					WDiff[m][n][d] = generator.nextGaussian() * sqrt(dt);
			Pfadgenerieren(X[m], WDiff[m], 0, X0);
		}
	}

	printf("X erstellen ... fertig\n");
}

void AmericanOption::testen() {
	double erg = 0;
	int lauf = 0;
	double** xx = DoubleFeld(N, D);
	double** wdiff = DoubleFeld(N, D);
	RNG generator;

	while (true) {
		for (int n = 0; n < N; ++n)
			for (int d = 0; d < D; ++d)
				wdiff[n][d] = sqrt(dt) * generator.nextGaussian();

		Pfadgenerieren(xx, wdiff);
		lauf++;
		int stopp = number_of_Exercise_Dates - 1;
		for (int n = 1; n < number_of_Exercise_Dates - 1; ++n)
			if (payoff(xx, Exercise_Dates[n])
					> LSM_C_estimated(xx[Exercise_Dates[n]], n)) { //TODO Achtung
				stopp = n;
				break;
			}
		erg = (erg * (double) (lauf - 1) + payoff(xx, Exercise_Dates[stopp]))
				/ (double) (lauf);

		if (lauf % 100 == 0) {
			printf("Stoppzeit testen ... %f %% \r", erg);
			cout.flush();
		}

	}
}

void AmericanOption::BFGS() {
	LongstaffSchwartz();
//	testen();
//	exit(0);
	Wdiff_und_X_erstellen();
//	stoppzeiten_erstellen();
//exit(0);
	StochInt_erstellen();

	BFGS_Iterations = 200;
	zeiger = this;

//	BFGS_setting();

	real_1d_array x;
	x.setlength(K);
	for (int k = 0; k < K; ++k)
		x[k] = 0;

	{
		printf("L-BFGS starten\n");
		M /= 50;
		double epsg = 0.000001;
		double epsf = 0;
		double epsx = 0;
		ae_int_t maxits = 20;
		mincgstate state;
		mincgreport rep;
		mincgcreate(x, state);
		mincgsetcond(state, epsg, epsf, epsx, maxits);
		alglib::mincgoptimize(state, Delegate_static_objfs);
		mincgresults(state, x, rep);
		printf("Termination Type(expected 4): %d\n", int(rep.terminationtype));
		M *= 50;
	}

	printf("L-BFGS starten\n");
	double epsg = 0.000001;
	double epsf = 0;
	double epsx = 0;
	ae_int_t maxits = 0;
	mincgstate state;
	mincgreport rep;
	mincgcreate(x, state);
	mincgsetcond(state, epsg, epsf, epsx, maxits);
	alglib::mincgoptimize(state, Delegate_static_objfs);
	mincgresults(state, x, rep);
	printf("Termination Type(expected 4): %d\n", int(rep.terminationtype));

	for (int k = 0; k < K; ++k)
		alpha[k] = x[k];

}

void AmericanOption::BFGS_extremeTesting(int l, double number_of_replications) {
	printf("extreme Testing ausgeschaltet\n");
//	if(D!=2)printf("ERROR2d");
//	Daten();
//
//	//Alphas ausgeben
//	neueExerciseDates(Testing_Dates);
//	MT.seed(time(NULL));
//	srand(time(NULL));
//	if(!nesterov && ! bfgs){
//		DateiLeser DL;
//		alpha=DL.alphasLaden(K);
//		if(verbose)printf("alphas geladen\n");
//		for(int k=0;k<K;++k)
//			if(verbose)printf("alpha[%d]: %f\n",k,alpha[k]);
//	}else if(verbose) printf("alphas vorhanden");
//
//	int n=l; //zeit
//	double** DV = (double**)malloc(sizeof(double*) * 100);
//	for(int i=0;i<100;++i)
//		DV[i]=(double*)malloc(sizeof(double)* 100);
//	for(double z=50;z<150;z+=1)
//		for(double y=50;y<150;y+=1)
//			DV[(int)y-50][(int)z-50]=0;
//
//	for(double z=50;z<150;z+=25)
//		for(double y=50;y<150;y+=25)
//		{
//			int ergebnispipe[Threadanzahl][2];
//			for(int t=0;t<Threadanzahl;++t)
//				pipe(ergebnispipe[t]);
//			int pid;
//			for(int f=0;f<Threadanzahl;++f)
//			{
//				pid = fork();
//				if (pid == 0)
//				{
//					close(ergebnispipe[f][0]);
//					char string[20];
//					srand(f);
//					MT.seed(getpid()+f+time(NULL));
//					double e=0;
//
//					double ** x=DoubleFeld(N,D);
//					double* start;
//					start=(double*)malloc(sizeof(double)*2);
//					start[0]=y;start[1]=z;
//
//					int durchlaeufe=500./Threadanzahl;
//					for (int k = 0; k < durchlaeufe; ++k){
//						if(verbose)printf("1\n");
//						Pfadgenerieren(x,l,start);
//						if(verbose)printf("2\n");
//						e+=AndersenBroadieEinzel(x,400,l)/(double)durchlaeufe;
//					}
//					sprintf(string,"%f",e);
//					write(ergebnispipe[f][1],string,(strlen(string)+1));
//					exit(0);
//				}
//			}
//
//			double erg=0;
//			double ergebnisse[Threadanzahl];
//			for(int f=0;f<Threadanzahl;++f)
//			{
//				close(ergebnispipe[f][1]);
//				char readbuffer[2000];
//				read (ergebnispipe[f][0],readbuffer, sizeof(readbuffer));
//				ergebnisse[f]=atof(readbuffer);
//				erg+=ergebnisse[f]/(double)(Threadanzahl);
//			}
//			printf("V_%d(%f,%f)=%f\n",l,z,y,erg);
//			DV[(int)y-50][(int)z-50]+=erg;
//		}
//
//	for(double z=50;z<150;z+=1)
//		for(double y=50;y<150;y+=1)
//		{
//			if(DV[(int)y-50][(int)z-50]!=0)continue;
//			DV[(int)y-50][(int)z-50]=DV[50][50];
//
//			double DT=0.1;
//			for(double t=0;t<1.;t+=DT){
//
//				double abl[2];
//				abl[0]=abl[1]=0;
//
//				double* x =(double*)malloc(sizeof(double)*2);
//				double* Xy=(double*)malloc(sizeof(double)*2);
//				int* reihe=(int*)   malloc(sizeof(int)*2);
//				x[0]=100.*(1-t)+y*t;x[1]=100.*(1-t)+z*t;
//
//				reihe[0]=0; reihe[1]=1;
//
//				BubbleSort(x,reihe,D);
//
//				for(int jj=0;jj<D;++jj)
//					Xy[jj]=log(x[jj]/Strike)/sqrt(T-n*dt+0.00001);
//
//				//				printf("%f,%f\n",x[0],sqrt(T-n*dt+0.00001));
//				//				printf("%f,%f\n",Xy[0],Xy[1]);
//
//				double s=0;
//
//				if(option==MAX_CALL)
//					for(int jjj=0; jjj<D/2+1;++jjj)
//						s+=Xy[reihe[jjj]];
//				if(option==MIN_PUT)
//					for(int jjj=D-1; jjj>=0;--jjj)
//						s+=Xy[reihe[jjj]];
//
//				for (int k = 0; k < K; ++k){
//					for(int j = 0; j < D; ++j){
//						double W=0;
//						//ANFANG Basisfunktionen //TODO
//						if(k<K1){//Basisfunktionen von Typ 1
//							int j_extrem=0;
//							if(D>1){
//								if(option==MAX_CALL)j_extrem=reihe[0];
//								if(option==MIN_PUT) j_extrem=reihe[D-1];
//							}
//							if(j==j_extrem)
//								W=F(Xy[j],k,0.8,true);
//						}
//
//						if(k>=K1 && k<K1+K2)//Basisfunktionen von Typ 2
//							//							if(Xy[0]>0. && Xy[1]>0.)
//							W=F(Xy[j],k-K1,0.8,false);
//
//						if(k>=K2+K1 && k<K1+K2+K3)//Basisfunktionen von Typ 3
//						{
//							bool nehmen=false;
//							if(option==MAX_CALL)
//								for(int jjj=0; jjj<D/2+1;++jjj)
//									if(j==reihe[jjj])nehmen=true;
//							if(option==MIN_PUT)
//								for(int jjj=D-1; jjj>=0;--jjj)
//									if(j==reihe[jjj])nehmen=true;
//							if(nehmen)
//								W=F(s,k-K1-K2,1.6,false);
//						}
//
//						if(k>=K1+K2+K3 && k<K1+K2+K3+K4)//Basisfunktionen vom Typ 4
//						{
//							int j_1,j_2;
//							if(option==MAX_CALL){
//								j_1=reihe[0];
//								j_2=reihe[1];
//							}
//							if(option==MIN_PUT){
//								j_1=reihe[D-2];
//								j_2=reihe[D-1];
//							}
//							double faktor=0;
//							if(j==j_1)faktor=1;
//							if(j==j_2)faktor=-1;
//
//							W=F(Xy[j_1]-Xy[j_2],k-K1-K2-K3,0.9,true)*faktor;
//						}
//
//						if(k>=K1+K2+K3+K4)//Basisfunktionen vom Typ 5
//						{
//							int t=(int)((double)n/(double)N*(double)K5/3.);
//							//printf("t=%d\n",t);
//							if( k-K1-K2-K3-K4>=5*t && k-K1-K2-K3-K4<5*(t+1))
//							{
//								//	printf("nr %d\n",k-K1-K2-K3-K4-5*t);
//								if(k-K1-K2-K3-K4-5*t==0)
//								{
//									if(j==0)W=1;
//								}
//								if(k-K1-K2-K3-K4-5*t==1)
//								{
//									if(j==1)W=1;
//								}
//								if(k-K1-K2-K3-K4-5*t==2)
//								{
//									if(j==0)W=Xy[1];
//									if(j==1)W=Xy[0];
//								}
//								//						if(k-K1-K2-K3-K4-5*t==3)
//								//						{
//								//							if(j==0)W=pow(Xy[1],2)*Xy[0];
//								//							if(j==1)W=pow(Xy[0],2)*Xy[1];
//								//						}
//								//						if(k-K1-K2-K3-K4-5*t==4)
//								//						{
//								//							W=1;
//								//						}
//							}
//						}
//						//ENDE Basisfunktionen //TODO
//						abl[j]+=alpha[k]*W;
//					}
//				}
//				DV[(int)y-50][(int)z-50]+=DT*(abl[0]*(y-100.) + abl[1]*(z-100.));
//			}
//		}
//
//	DateiLeser DL;
//	DL.ausgeben(l,DV,100,100);
//	printf("%d\n",n);
}
