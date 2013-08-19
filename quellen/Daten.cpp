#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten(){
	int Example=3;

	X0=(double*)malloc(sizeof(double)*100);
	sigma=(double*)malloc(sizeof(double)*100);

//	if(Example==1){					//
//			PfadModell=ITO;
//			option=MIN_PUT;
//			delta=0;
//			D=1;
//			for(int j=0;j<D;++j)
//			{
//				X0[j]=100.;
//				sigma[j]=0.4;
//			}
//			Strike = 100.;
//			r = 0.06;
//			T = 0.5;
//			Dates=10;
//			N = 100;
//			KpI=7;
//			M=1000;
//		}

	if(Example==3){					//Glasserman Example MaxCall
		PfadModell=ITO;
		option=MAX_CALL;
		delta=0.1;
		D=2;
		for(int j=0;j<D;++j)
		{
			X0[j]=90.;
			sigma[j]=0.2;
		}
		Strike = 100.;
		r = 0.05;
		T = 3;
		N=90;
		Dates = 9;
		KpI=20;
		M=10000;
	}

    BFGS_Nesterov_Intervals=3;
	Threadanzahl=10;

	NpI=N/BFGS_Nesterov_Intervals;
}
