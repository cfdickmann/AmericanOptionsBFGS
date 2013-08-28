#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten() {
	int Example =3;

	X0 = (double*) malloc(sizeof(double) * 100);
	sigma = (double*) malloc(sizeof(double) * 100);

	if (Example == 1) {            //Rogers Put
		PfadModell = ITO;
		option = MIN_PUT;
		delta = 0;
		D = 1;
		for (int j = 0; j < D; ++j) {
			X0[j] = 100.;
			sigma[j] = 0.4;
		}
		Strike = 100.;
		r = 0.06;
		T = 0.5;
		Dates = 50;
		N = 150;
		KpI = 20;
		M = 5000;
	}

	if (Example == 3) {					//Glasserman Example MaxCall
		PfadModell = ITO;
		option = MAX_CALL;
		delta = 0.1;
		D = 2;
		for (int j = 0; j < D; ++j) {
			X0[j] = 90.;
			sigma[j] = 0.2;
		}
		Strike = 100.;
		r = 0.05;
		T = 3;
		N = 2*90;
		Dates = 10;
		KpI = 48;
		M = 10000;
	}

	BFGS_Nesterov_Intervals = 1;
	Threadanzahl = 10;
	NpI = N / BFGS_Nesterov_Intervals;
}


//	if (Example == 2) {					//Rogers MinPut
//		PfadModell = ITO;
//		option = MIN_PUT;
//		delta = 0;
//		D = 2;
//		for (int j = 0; j < D; ++j) {
//			X0[j] = 100.;
//			sigma[j] = 0.6;
//		}
//		Strike = 100.;
//		r = 0.06;
//		T = 0.5;
//		N = 400;
//		Dates = 20;
//		KpI = 1;
//		M = 50000;
//	}

