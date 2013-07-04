#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten(){
	int Example=3;

	X0=(double*)malloc(sizeof(double)*100);
	sigma=(double*)malloc(sizeof(double)*100);


	if(Example==1){				//Rogers Example 1d
		PfadModell=ITO;
		//	PfadModell=EULER;
		option=MIN_PUT;
		delta=0;
		X0[0] = 40.;      		//Spot
		Strike = 100.; 			//Ausuebungspreis
		r = 0.06;   			//interest rate
		sigma[0] = 0.4;  		//Volatility
		T = 0.5; 			//Gesamtzeit
		Testing_Dates=10;
		Training_Dates=10;
		N =10;  			//time discretization
		D=1;
		K1=20;                          // Basisfunktionen fuer
		K2=0;                           // BFGS und Nesterov
		K3=0;
		K4=0;
		K5=0;
		M=15000;   		        //Number of replications fuer BFGS und Nesterov
	}

	if(Example==3){					//Glasserman Example MaxCall
		PfadModell=ITO;
		option=MAX_CALL;
		delta=0;
		D=1;
		for(int j=0;j<D;++j)
		{
			X0[j]=100.;
			sigma[j]=0.2;
		}
		Strike = 100.;
		r = 0;
		T = 3;
		Testing_Dates=10;
		Training_Dates=10;
		N = 101;
		K1=1;
		K2=1;
		K3=1;
		K4=1;
		K5=1;
		M=5000;
	}

        BFGS_Nesterov_Intervals=10;
	Threadanzahl=10;

	KpI=K1+K2+K3+K4+K5;
	NpI=(N-1)/BFGS_Nesterov_Intervals;
}
