#include "../alglib/optimization.h"
#include "../alglib/linalg.h"
#include "../alglib/ap.h"

#include "AmericanOption.h"

using namespace std;
using namespace alglib;

static AmericanOption* zeiger;

double AmericanOption::obj(double * alpha) {
	//Training
	double ergETr = 0;
	double erg = 0;
	double ergQ = 0;
	double p = 5.; //Glättungsparameter
	for (int m = 0; m < M / 2; ++m) { //Summieren  ueber alle replications
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
//			d[ex] -= GrundIntegrals[m][Exercise_Dates[ex]];
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][Exercise_Dates[ex]];
		}
		ergETr += Max(d, number_of_Exercise_Dates) / (double) (M / 2);
		double s = 0;
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
			s += exp((p * d[ex]));
		double u = log(s) / p;
		erg += u / (double) (M / 2);
		ergQ += u * u / (double) (M / 2);
	}
	double var = (ergQ - erg * erg);
	double V = erg + 0.001 * var;

	static int zaehler = 0;
	zaehler++;
	if (zaehler % K == 0) {
		//Testing
	double ergETe = 0;
	int stat[number_of_Exercise_Dates];
	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		stat[e] = 0;

	for (int m = M / 2; m < M; ++m) {
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
//			d[ex] -= GrundIntegrals[m][Exercise_Dates[ex]];
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][Exercise_Dates[ex]];
		}
		ergETe += Max(d, number_of_Exercise_Dates) / (double) (M / 2);
		stat[argMax(d, number_of_Exercise_Dates)]++;
	}

	//Ausgabe
		for (int k = 0; k < K; ++k)
			printf("%.6lf,", alpha[k]);
		printf("\nV=%f(ergGlatt: %f + lambda* var %f) Training_exakt:%f:\n", V,erg,var,ergETr);
		printf("dual(testing): %f\n", ergETe);

		for (int e = 0; e < number_of_Exercise_Dates; ++e)
			printf("stat %d: %d,", e, stat[e]);
		printf("\n;");

		if ((zaehler % K )%5 == 0) {printf("schreiben\n");
		fstream f;
			f.open("mart.txt", ios::out);

		double** xx = DoubleFeld(N, D);
		for (int n = 0; n < number_of_Exercise_Dates; ++n) {
			for (double x = 1; x < 200; x+=3) {
				xx[Exercise_Dates[n]][0] = x;
				double func = 0;
//				func= grund(xx[n],n)/x;
				for (int k = 0; k < K; ++k)
					func += alpha[k] * this->f(k, xx, Exercise_Dates[n], 0) / x;
				f << func << endl;
			}
			f << endl;
		}
		f.close();
		deleteDoubleFeld(xx,N,D);
//		system("gnuplot plot.gnuplot");
		}
	}

	return V;
}

void AmericanOption::stoppzeiten_erstellen() {
	stoppzeiten = IntFeld(M);
	int stat[N];
	for (int n = 0; n < N; ++n)
		stat[n] = 0;
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
		stat[Exercise_Dates[stopp]]++;
		erg += payoff(X[m], Exercise_Dates[stopp]);
		stoppzeiten[m] = Exercise_Dates[stopp];
	}

	printf("Stoppzeiten erstellen ... fertig\n");
	printf("stoppzeiten durchschnitt %f \n", erg / (double) M);
//
//	for (int n = 0; n < N; ++n)
//		printf("Stoppzeit %d: %d\n", n, stat[n]);
}

//Mit numerischer Differentiation
static void Delegate_static_objfs(const real_1d_array &x, double &func,
		real_1d_array &grad, void *ptr) {
	int K = zeiger->K;
	double alpha[K];
	for (int k = 0; k < K; ++k)
		alpha[k] = x[k];

	double gradient[K];
	zeiger->objfs_aufrufen(alpha, func, gradient);

	for (int k = 0; k < K; ++k)
		grad[k] = gradient[k];

}

void AmericanOption::objfs_aufrufen(double* alpha, double &func, double* grad) {
	func = obj(alpha);
	for (int k = 0; k < K; ++k) {
		double h = 0.0000001;
		alpha[k] += h;
		double fh = obj(alpha);
		grad[k] = (fh - func) / h;
		alpha[k] -= h;
	}
}

//double AmericanOption::F(double x, int k, double border, bool hauf) {
//	if (k == 0)
//		return 1;
//	double extra = (x >= border && hauf) ? 1 : 0;
//	if (k % 2 == 0)
//		return cos((double) k / 2. * x) * (fabs(x) < border) + extra;
//	else
//		return sin(((double) k + 1.) * x / 2.) * (fabs(x) < border) + extra;
//}

void AmericanOption::StochInt_erstellen() {
	StochIntegrals = DoubleFeld(K, M, N);
	GrundIntegrals = DoubleFeld(M, N);

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
				StochIntegrals[k][m][n] = StochIntegrals[k][m][n - 1] + summe;
			}
		}

//		for (int k = 0; k < K; ++k) {
//			for (int n = 1; n < N; ++n) {
//				double summe = 0;
//				for (int d = 0; d < D; ++d)
//					summe += grund(X[m][n - 1], n - 1) * WDiff[m][n][d];
//				GrundIntegrals[m][n] = GrundIntegrals[m][n - 1] + summe;
//			}
//		}
	}
	printf("StochInt erstellen ... fertig\n");
}

//void AmericanOption::BFGS_parallelTesting(double number_of_replications) {
//	Daten();
//	int ergebnispipe[Threadanzahl][2];
//	for (int z = 0; z < Threadanzahl; ++z)
//		pipe(ergebnispipe[z]);
//
//	int pid;
//	for (int f = 0; f < Threadanzahl; ++f) {
//		pid = fork();
//		if (pid == 0) {
//			close(ergebnispipe[f][0]);
//			char string[20];
//			srand(f);
//			double e = BFGS_testing(50000. / (double) Threadanzahl);
//			sprintf(string, "%f", e);
//			write(ergebnispipe[f][1], string, (strlen(string) + 1));
//			exit(0);
//		}
//	}
//	double erg = 0;
//	double ergebnisse[Threadanzahl];
//	for (int f = 0; f < Threadanzahl; ++f) {
//		close(ergebnispipe[f][1]);
//		char readbuffer[2000];
//		read(ergebnispipe[f][0], readbuffer, sizeof(readbuffer));
//		ergebnisse[f] = atof(readbuffer);
//		erg += ergebnisse[f] / (double) (Threadanzahl);
//		printf("Ergebnis %d: %f\n", f, atof(readbuffer));
//	}
//	double s = 0;
//	for (int f = 0; f < Threadanzahl; ++f)
//		s += pow(ergebnisse[f] - erg, 2);
//	printf("geschaetzte Varianz: %f\n", s / double(Threadanzahl - 1));
//	printf("Gesamtergebnis: %f\n", erg);
//	ErgebnisAnhaengen(erg);
//}

//double AmericanOption::BFGS_testing(double number_of_replications) {
////	neueExerciseDates(Testing_Dates); //Achtung fehlt
//	MT.seed(time(NULL) + getpid());
//	srand(time(NULL) + getpid());
//	if (!nesterov && !bfgs) {
//		alpha = alphasLaden(K);
//		if (verbose && !parallelTest)
//			printf("alphas geladen\n");
//	} else if (verbose && !parallelTest)
//		printf("alphas vorhanden");
//
////Alphas ausgeben
//	if (verbose && !parallelTest) {
//		printf("alpha=[");
//		for (int k = 0; k < K; ++k)
//			printf("%f, ", alpha[k]);
//		printf("]\n");
//	}
//
////Felder fuer je 1 replication
//	double** stochintegrals = DoubleFeld(K, N);
//	double** x = DoubleFeld(N, D);
//	double** wdiff = DoubleFeld(N, D);
//	double** sprue = DoubleFeld(N, D);
//
//	double mean = 0;
//
//	int Npfade = 15;
//	double w[N * Npfade]; //fuer die ausgabe von pfaden
//
//	for (int lauf = 1; lauf <= number_of_replications; lauf++) {
//		for (int n = 0; n < N; ++n)
//			for (int j = 0; j < D; ++j)
//				wdiff[n][j] = sqrt(dt) * nextGaussian();
//
//		Pfadgenerieren(x, wdiff);
//
//		if (lauf < Npfade + 1) {
//			for (int n = 0; n < N; ++n)
//				w[n + (lauf - 1) * N] = x[n][0];
//		}
//		if (lauf == Npfade + 1) {
//			werteSchreiben(w, N * Npfade, N);
//		}
//
//		for (int k = 0; k < K; ++k)
//			stochintegrals[k][0] = 0;
////		StochInt(stochintegrals, x, wdiff, sprue);
//
//		//for(int n=0;n<N;++n)
//		//	printf("%.2lf ",stochintegrals[299][n]);
//
//		double d[number_of_Exercise_Dates];
//		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
//			d[ex] = payoff(x[Exercise_Dates[ex]], Exercise_Dates[ex]);
//			for (int k = 0; k < K; ++k)
//				d[ex] -= alpha[k] * stochintegrals[k][Exercise_Dates[ex]];
//		}
//
//		mean += Max(d, number_of_Exercise_Dates);
//
//		if (lauf % 100 == 0)
//			printf(
//					"Zwischenergebnis (kumuliert): %f bei %d Pfaden (%d Prozent fertig)\r",
//					mean / ((double) (lauf)), lauf,
//					(int) floor(
//							(double) (lauf) / number_of_replications * 100.));
//	}
//
//	if (!parallelTest) {
//		printf("\nresult: %f \n", mean / number_of_replications);
//		ErgebnisAnhaengen(mean / number_of_replications);
//	}
//	return mean / number_of_replications;
//}

//void AmericanOption::StochInt(double** STi, double** X, double** WDiff,
//		double** Sprue) {
//	for (int n = 1; n < N; ++n)
//		for (int k = 0; k < K; ++k) {
//			double summe = 0;
//			for (int d = 0; d < D; ++d)
//				summe += f(k, X, n - 1, d) * WDiff[n][d];
//			STi[k][n] = STi[k][n - 1] + summe;
//		}
//
//	for (int n = 1; n < N; ++n)
//		for (int k = 0; k < K; ++k) {
//			double summe = 0;
//			for (int d = 0; d < D; ++d)
//				summe += grund(x[n-1],n-1)* WDiff[n][d];
//			 GrundIntegrals[k][n] = GrundIntegrals[k][n - 1] + summe;
//		}
//}

void AmericanOption::Wdiff_und_X_erstellen() {
	RNG generator;
	//generator.setSeed(seed);
	WDiff = DoubleFeld(M, N, D);
	X = DoubleFeld(M, N, D);

	for (int m = 0; m < M; ++m) {
		if (m % 10 == 0) {
			printf("X erstellen ... %d %%            \r", m * 100 / M);
			cout.flush();
		}

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
	stoppzeiten_erstellen();
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
		M /= 100;
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
		M *= 100;
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
