#include "AmericanOption.h"
#include "RNG.h"

using namespace std;

AmericanOption::AmericanOption() {
	BFGS_Nesterov_Intervals = 1;
	Iterations_Nummer = 0;
 nesterov = false;
	bfgs =  longstaffschwarz = false;
	Daten();
	dt = T / (double) (N - 1);
	neueExerciseDates(Dates);
	K = BFGS_Nesterov_Intervals * (KpI); //Anzahl der Basisfunktionen

	Anfangszeit = time(NULL);
	//(double*)malloc(sizeof(double)*Threadanzahl);
}

AmericanOption::~AmericanOption() {

}

double AmericanOption::BoxMuller(double U1, double U2) {
	double R = -2 * log(U1);
	double V = 2 * 3.1415926535 * U2;
	return sqrt(R) * cos(V);
}

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

void AmericanOption::neueExerciseDates(int n) {
	Exercise_Dates = new int[n + 1];
	number_of_Exercise_Dates = n;

	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		Exercise_Dates[e] =  (N - 1) *  e / (number_of_Exercise_Dates - 1);

	printf("neue Exercise dates: ");

	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		printf("%f, ", (double) Exercise_Dates[e]);
	printf("\n");

	for (int e = 0; e < number_of_Exercise_Dates; ++e)
		printf("%f, ", (double) Exercise_Dates[e]*dt);
	printf("\n");
}

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
	}
		if (X[N - 1][0] <= 0)
			printf("Error0\n");
	}
}
