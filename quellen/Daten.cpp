#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten(){
	int Example=3;

	X0=(double*)malloc(sizeof(double)*100);
	sigma=(double*)malloc(sizeof(double)*100);

	if(Example==3){					//Glasserman Example MaxCall
		PfadModell=ITO;
		option=MAX_CALL;
		delta=0*0.1;
		D=1;
		for(int j=0;j<D;++j)
		{
			X0[j]=100.;
			sigma[j]=0.2;
		}
		Strike = 100.;
		r = 0.05;
		T = 3;
//		Testing_Dates=40;
//		Training_Dates=40;
		N = 8;
		KpI=6;
		M=1000;
	}

    BFGS_Nesterov_Intervals=4;
	Threadanzahl=10;


	NpI=N/BFGS_Nesterov_Intervals;
}
