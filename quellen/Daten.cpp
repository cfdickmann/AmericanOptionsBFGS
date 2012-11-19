#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten(){
	int Example=103;

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

	if(Example==101){				//Rogers Example 1d
			PfadModell=ITO;
			//	PfadModell=EULER;
			option=MIN_PUT;
			delta=0;
			X0[0] = 40.;      		//Spot
			Strike = 45.; 			//Ausuebungspreis
			r = 0.0488;   			//interest rate
			sigma[0] = 0.2;  		//Volatility
			T = 1./3.; 			//Gesamtzeit
			Testing_Dates=20;
			Training_Dates=20;
			N =20;  			//time discretization
			D=1;
			K1=20;                          // Basisfunktionen fuer
			K2=0;                           // BFGS und Nesterov
			K3=0;
			K4=0;
			K5=0;
			M=15000;   		        //Number of replications fuer BFGS und Nesterov
		}

	if(Example==2){					//Rogers Example MIN_PUT
		PfadModell=ITO;
		//	PfadModell=EULER;
		option=MIN_PUT;
		delta=0;
		X0[0] =120;  sigma[0]=0.4;
		X0[1]= 80; sigma[1]=0.8;
		Strike = 100.;
		r = 0.06;
		T = 0.5;
		Testing_Dates=100;
		Training_Dates=50;
		N = 401;
		D=2;
		K1=15; 
		K2=15;
		K3=15;
		K4=15;
		K5=0;
		M=20000;
	}

	if(Example==103){					//Glasserman Example MaxCall
			PfadModell=ITOrho;
			option=MAX_CALL;
			delta=0.1;
			D=3;
			rho=0.3;
			for(int j=0;j<D;++j)
			{
				X0[j]=110.;
				sigma[j]=0.2;
			}
			Strike = 100.;
			r = 0.05;
			T = 1;
			Testing_Dates=4;
			Training_Dates=4;
			N = 4;
			K1=15;
			K2=15;
			K3=15;
			K4=15;
			K5=0;
			M=5000;
		}

	if(Example==3){					//Glasserman Example MaxCall
		PfadModell=ITO;
		option=MAX_CALL;
		delta=0.1;
		D=10;
		for(int j=0;j<D;++j)
		{
			X0[j]=110.;
			sigma[j]=0.2;
		}
		Strike = 100.;
		r = 0.05;
		T = 3;
		Testing_Dates=10;
		Training_Dates=10;
		N = 10;
		K1=15;
		K2=15;
		K3=15;
		K4=15;
		K5=0;
		M=5000;
	}
//
//	if(Example==4){					//Cuffignals Example
//		PfadModell=ITO;
//		option=MIN_PUT;
//		delta=0;
//		X0[0] = 10;
//		Strike = X0[0];
//		r = 0.05;
//		sigma[0] = 0.2;
//		T = 1;
//		Testing_Dates=128;
//		Training_Dates=128;
//		N = 128;
//		D=1;
//		K1=15;
//		K2=15;
//		K3=15;
//		K4=15;
//		K5=0;
//		M=8000;
//	}
//
//	if(Example==6){					//CIR Example
//		PfadModell=CIR;
//		option=MAX_CALL;
//		delta=0;
//		kappa=0.5;
//		theta=1.5;
//		D=2;
//		for(int j=0;j<D;++j)
//		{
//			X0[j]= theta;
//			sigma[j]=0.2;
//		}
//		Strike = 1.;
//		r = 0.05;
//		T = 1;
//		Testing_Dates=10;
//		Training_Dates=10;
//		N = 201;
//		K1=15;
//		K2=0;
//		K3=0;
//		K4=0;
//		K5=0;
//		M=15000;
//	}
//
//	if(Example==7){					//firth
//		PfadModell=ITO;
//		eta=0;
//		option=MIN_PUT;
//		delta=0;
//		lambdaJump=1;
//		D=1;
//		for(int j=0;j<D;++j)
//		{
//			X0[j]= 36.;
//			sigma[j]=0.2;
//		}
//		Strike = 40.;
//		r = 0.06;
//		T = 2;
//		Testing_Dates=10;
//		Training_Dates=10;
//		N = 401;
//		//		K1=K2=K3=K4=15;
//		K1=10;
//		K2=10;
//		K3=10;
//		K4=10;
//		//		K5=J*3;
//		K5=0;
//		M=20000;
//	}
//
//	if(Example==8){					//Jump diffusion
//		PfadModell=JDI;
//		eta=0.1;
//		option=MAX_CALL;
//		delta=0.1;
//		lambdaJump=1;
//		D=10;
//		for(int j=0;j<D;++j)
//		{
//			X0[j]= 90.;
//			sigma[j]=0.2;
//		}
//		Strike = 100.;
//		r = 0.05;
//		T = 3;
//		Testing_Dates=10;
//		Training_Dates=10;
//		N = 201;
//		//		K1=K2=K3=K4=15;
//		K1=15;
//		K2=10;
//		K3=10;
//		K4=10;
//		//		K5=J*3;
//		K5=0;
//		M=8000;
//	}

        BFGS_Nesterov_Intervals=1;
	Threadanzahl=8;
}
