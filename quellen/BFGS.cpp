#include <fstream>

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
	for (int m = 0; m < M / 2; ++m) { //Summieren  ueber alle replications
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
//			d[ex] -= GrundIntegrals[m][Exercise_Dates[ex]];
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][ex];
		}
		ergETr += Max(d, number_of_Exercise_Dates) / (double) (M / 2);
//		double s = 0;
//		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
//			s += exp((p * d[ex]));
//		double u = log(s) / p;
		double u = Max(d, number_of_Exercise_Dates);
		erg += u / (double) (M / 2);
		ergQ += u * u / (double) (M / 2);
	}
	double var = (ergQ - erg * erg);
	double V = erg + 0.00000 * var;
//double V=ergETr;

	static int zaehler = 0;
	zaehler++;
//	if (zaehler % K == 0)
	{
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
					d[ex] -= alpha[k] * StochIntegrals[k][m][ex];
			}
			ergETe += Max(d, number_of_Exercise_Dates) / (double) (M / 2);
			stat[argMax(d, number_of_Exercise_Dates)]++;
		}

		//Ausgabe
		for (int k = 0; k < K; ++k)
			printf("%.6lf,", alpha[k]);
		printf("\nV=%f(ergGlatt: %f + lambda* var %f) Training_exakt:%f:\n", V,
				erg, var, ergETr);
		printf("Testing_exakt: %f\n", ergETe);
		cout.flush();

//		for (int e = 0; e < number_of_Exercise_Dates; ++e)
//			printf("stat %d: %d,", e, stat[e]);
//		printf("\n;");

////		if ((zaehler / K) % 5 == 0)
//		static int zzz = 0;
//		zzz++;
//		if (zzz % 10 == 0) {
//			printf("schreiben\n");
//			fstream f;
//			f.open("mart.txt", ios::out);
//
//			double** xx = DoubleFeld(N, D);
//			double** yy = DoubleFeld(N, D);
//			int** rr = IntFeld(N, D);
//			for (int n = 0; n < number_of_Exercise_Dates; ++n) {
//				for (double x = 20; x < 170; x += 3) {
//					xx[Exercise_Dates[n]][0] = x;
//					rr[Exercise_Dates[n]][0] = 0;
//					yy[Exercise_Dates[n]][0] = log(x / Strike)
//							/ sqrt(T - n * dt + 0.001);
//					double func = 0;
////				func= grund(xx[n],n)/x;
//
//					for (int k = 0; k < K; ++k) {
//						double* all = this->f_all(k, xx, yy, rr, Exercise_Dates[n]); //Hier lauert ein segmentation fault
//						if (all != NULL) {
//							func += alpha[k] * all[0] / x;
//							delete[] all;
//						}
//					}
//					f << func << endl;
//				}
//				f << endl;
//			}
//			f.close();
//			deleteDoubleFeld(xx, N, D);
//			deleteIntFeld(rr, N, D);
//			deleteDoubleFeld(yy, N, D);
////		system("gnuplot plot.gnuplot");
//		}
	}

	return V;
}

double* AmericanOption::obj_diff(double * alpha) { // mit Glaettungsparameter ohne varianz
	double* diffs = new double[K];

	for (int k = 0; k < K; ++k)
		diffs[k] = 0;

	for (int m = 0; m < M / 2; ++m) { //Summieren  ueber alle replications
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][ex];
		}

		double s = 0;
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
			s += exp((p * d[ex]));

		for (int k = 0; k < K; ++k) {
			double S = 0;
			for (int ex = 0; ex < number_of_Exercise_Dates; ++ex)
				S += StochIntegrals[k][m][ex] * exp(p * d[ex]);
			diffs[k] += -S / s / (double) (M / 2);
		}
	}
	return diffs;
}
//
//double* AmericanOption::obj_diff(double * alpha) { //ohne Glaettungsparameter
//	double* diffs = new double[K];
//
//	for (int k = 0; k < K; ++k)
//		diffs[k] = 0;
//
//	for (int m = 0; m < M / 2; ++m) { //Summieren  ueber alle replications
//		double d[number_of_Exercise_Dates];
//		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
//			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
////			d[ex] -= GrundIntegrals[m][Exercise_Dates[ex]];
//			for (int k = 0; k < K; ++k)
//				d[ex] -= alpha[k] * StochIntegrals[k][m][ex];
//		}
//		int n_max = argMax(d, number_of_Exercise_Dates);
//		for (int k = 0; k < K; ++k)
//			diffs[k] += -StochIntegrals[k][m][n_max]
//					/ (double(M) * 2);
//	}
//	return diffs;
//}

void AmericanOption::stoppzeiten_erstellen() {
	stoppzeiten = IntFeld(M);
	int stat[number_of_Exercise_Dates];
	for (int n = 0; n < number_of_Exercise_Dates; ++n)
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
		stat[stopp]++;
		erg += payoff(X[m], Exercise_Dates[stopp]);
		stoppzeiten[m] = Exercise_Dates[stopp];
	}

	printf("Stoppzeiten erstellen ... fertig\n");
	printf("stoppzeiten durchschnitt %f \n", erg / (double) M);
//
	printf("Stoppzeiten: \n");
	for (int n = 0; n < number_of_Exercise_Dates; ++n)
		printf("%d,", stat[n]);
	printf("\n");
//	exit(0);
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

//	double grad2[K];
//	for (int k = 0; k < K; ++k) {
//		double h = 0.0000001;
//		alpha[k] += h;
//		double fh = obj(alpha);
//		grad[k] = (fh - func) / h;
//		alpha[k] -= h;
//	}

	double* diffs = obj_diff(alpha);
	for (int k = 0; k < K; ++k)
		grad[k] = diffs[k];
	delete[] diffs;
//
//	for (int i = 0; i < K; ++i)
//		printf("%f und %f,\n", grad[i], grad2[i]);
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

double cndf(double x) {
	return 0.5 * (1. + erf(x / sqrt(2)));
}

void AmericanOption::Wdiff_und_X_erstellen() {
	RNG generator;
	//generator.setSeed(seed);
	WDiff = DoubleFeld(M, N, D);
	X = DoubleFeld(M, N, D);
	Y = DoubleFeld(M, N, D);
	reihe=IntFeld(M,N,D);

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
	for (int m = 0; m < M; ++m) {
		if (m % 10 == 0) {
			printf("Y erstellen ... %d %%            \r", m * 100 / M);
			cout.flush();
		}

		{
			for (int n = 0; n < N; ++n)
				for (int d = 0; d < D; ++d)
					Y[m][n][d] = log(X[m][n][d] / Strike)
							/ sqrt(T - (double) n * dt + 0.00001);
//				Pfadgenerieren(X[m], WDiff[m], 0, X0);
		}
	}

	printf("Y erstellen ... fertig\n");
	for (int m = 0; m < M; ++m) {
		if (m % 10 == 0) {
			printf("reihe erstellen ... %d %%            \r", m * 100 / M);
			cout.flush();
		}
			for (int n = 0; n < N; ++n)
				reihe[m][n]=BubbleSort(X[m][n],D);
	}

	printf("reihe erstellen ... fertig\n");
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

void AmericanOption::CG(double* alpha) {
	real_1d_array x;
	x.setlength(K);
	for (int k = 0; k < K; ++k)
		x[k] = alpha[k];

	printf("CG starten\n");
	double epsg = 0.000001;
	double epsf = 0;
	double epsx = 0;
	ae_int_t maxits = 10;
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

void AmericanOption::BFGS(double* alpha) {
	real_1d_array x;
	x.setlength(K);
	for (int k = 0; k < K; ++k)
		x[k] = alpha[k];

	printf("L-BFGS starten\n");
//	double epsg = 0.000001;
//	double epsf = 0;
//	double epsx = 0;
	ae_int_t maxits = 0;
//	mincgstate state;
//	mincgreport rep;
//	minlbfgscreate(x, state);
//	mincgsetcond(state, epsg, epsf, epsx, maxits);
//	alglib::minlbfgsoptimize(state, Delegate_static_objfs);
//	mincgresults(state, x, rep);
//	printf("Termination Type(expected 4): %d\n", int(rep.terminationtype));

	double epsg = 0.000001;
	double epsf = 0.000001;
	double epsx = 0.000001;
//    double diffstep = 1.0e-6;
//    ae_int_t maxits = 0;
	minlbfgsstate state;
	minlbfgsreport rep;

	minlbfgscreate(1, x, state);
	minlbfgssetcond(state, epsg, epsf, epsx, maxits);
	alglib::minlbfgsoptimize(state, Delegate_static_objfs);
	minlbfgsresults(state, x, rep);

	printf("Termination Type(expected 4): %d\n", int(rep.terminationtype));

	for (int k = 0; k < K; ++k)
		alpha[k] = x[k];
}

void AmericanOption::BFGS_aufrufen() {
	LongstaffSchwartz();
//	testen();
//	exit(0);
	Wdiff_und_X_erstellen();
	stoppzeiten_erstellen();
	StochInt_erstellen();

	zeiger = this;
	BFGS_Iterations = 20000;

	double alpha[K];

	for (int k = 0; k < K; ++k)
		alpha[k] = 0;

//	for (int k = 0; k < K; ++k)
//		alpha[k] = 1.0/(double)(K)*sigma[0];
//
//
//	double gamma[29] = { 0.270214964, 0.262294286, 0.002871799, -0.246069280,
//			-0.399297234, -0.170046342, -0.171978632, 0.172005816, 0.645239092,
//			0.236397054, -0.313426534, -0.162664364, -0.367309602, 1.350555822,
//			0.324501753, -0.602584800, 0.460306357, -1.600847590, -0.087618737,
//			-0.262300287, -0.768527750, 2.017575277, 0.460799091, -0.859946013,
//			0.023943856, -0.058172667, -0.148199704, 0.059818865, 0.091840144 };

//	double beta[55]={
//	0.277948720416063,
//	0.24423532258458316,
//	0.21498579753441713,
//	0.30548029127014187,
//	0.07596160010792914,
//	0.17717133764387039,
//	-0.0359596498149467,
//	0.019051028547944274,
//	-0.05694838135839331,
//	-0.02049939469416402,
//	-0.010284431548030496,
//	0.049726063784567696,
//	0.032027443847918224,
//	0.10354965626819945,
//	0.026277038025868713,
//	0.06107121048688391,
//	-0.008343681790585648,
//	-0.021872505829892808,
//	-0.025457205202658727,
//	-0.03802663161950771,
//	-0.00854632060651845,
//	0.02728098347057547,
//	0.01577893523858833,
//	0.07760747607627068,
//	0.016385688818893584,
//	0.03653242590788444,
//	-0.00393480076742058,
//	-0.04411998993559421,
//	-0.016285151685483782,
//	-0.05046490146171571,
//	-0.0054576661573987675,
//	0.030988578473800325,
//	0.011431923977195769,
//	0.08050450830322112,
//	0.01051058005751125,
//	0.012303428864031301,
//	-0.005327122203231325,
//	-0.07962719634229135,
//	-0.011703493579703182,
//	-0.04004836697642605,
//	4.7677198759575406E-4,
//	0.08626722255828474,
//	0.010819147033526062,
//	0.07060081345218915,
//	0.0019340291683042287,
//	-0.10252942638945799,
//	-0.009827666667450573,
//	-0.07503242902380206,
//	-0.0011635021521699568,
//	0.19304402802457313,
//	0.009745194009910483,
//	-0.12865494899791535,
//	-0.006490612666559615,
//	0.03330476034858769,
//	0.0014488325711449366
//	};
//
////	for (int k = 0; k < K; ++k)
////			alpha[k] = 1.0/(double)(K)*sigma[0];
//
//	printf("mit gegebenen koeff: %f\n", obj(gamma));
//exit(0);

//	printf("\n\nobjfs %f\n\n", obj(alpha));
//	exit(0);
	M /= 10;
	CG(alpha);
	M *= 10;
	BFGS(alpha);

	//Training Wert eintragen
	double ergETr = 0;

	for (int m = 0; m < M / 2; ++m) {
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][ex];
		}
		ergETr += Max(d, number_of_Exercise_Dates) / (double) (M / 2);
	}

	ErgebnisAnhaengen(ergETr, (char*) "BFGS_training.txt");

	//Testing Wert eintragen
	double ergETe = 0;
	for (int m = M / 2; m < M; ++m) {
		double d[number_of_Exercise_Dates];
		for (int ex = 0; ex < number_of_Exercise_Dates; ++ex) {
			d[ex] = payoff(X[m][Exercise_Dates[ex]], Exercise_Dates[ex]);
			for (int k = 0; k < K; ++k)
				d[ex] -= alpha[k] * StochIntegrals[k][m][ex];
		}
		ergETe += Max(d, number_of_Exercise_Dates) / (double) (M / 2);
	}

	ErgebnisAnhaengen(ergETe, (char*) "BFGS_testing.txt");

}
