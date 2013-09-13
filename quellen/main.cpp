#include "AmericanOption.h"

using namespace std;


int main(int argc, char* args[]) {

	EuroBewerter EB;

//	printf("Margrabe : %f\n", EB.exchange_option(150, 100., 0., 3., 0.05, 0, 0.2));
//	printf("Margrabe : %f\n", EB.exchange_option(300, 200., 0., 3., 0.05, 0, 0.2));

//
//	printf("test: %f\n", EB.exchange_option_diff   (150, 100., 0., 3., 0.05, 0, 0.2,0));
//	printf("test2: %f\n", EB.exchange_option_diff2 (150, 100., 0., 3., 0.05, 0, 0.2,0));
//	printf("test: %f\n", EB.exchange_option_diff   (150, 100., 0., 3., 0.05, 0, 0.2,1));
//	printf("test2: %f\n", EB.exchange_option_diff2 (150, 100., 0., 3., 0.05, 0, 0.2,1));

//		printf("test: %f\n", EB.exchange_option_diff(110, 100., 0., 3., 0.05, 0, 0.2,0));
//
//	printf("teueres y hinten: %f\n", EB.exchange_option(90., 100., 0., 3., 0.05, 0, 0.2));
//	printf("teueres x vorne: %f\n", EB.exchange_option(100., 90., 0., 3., 0.05, 0, 0.2));
		printf("E0  : %f\n",EB.put(0.,1.,110,100,0.07,0.4,0.25));
		printf("E3  : %f\n",EB.put(3.,4.,110,100,0.07,0.4,0.25));
		printf("E3  : %f\n",EB.put_diff(3.,4.,110,100,0.07,0.4,0.25));

	int runden = 1;
	AmericanOption AMO;

//	double d[5]={5,4,2,3,1};
//	int* reihe=BubbleSort(d,5);
//
//	for(int i=0;i<5;++i)
//	printf("%d,", reihe[i]);
//		printf("\n");
//		delete[] reihe;

//	if (argc == 1) {
//		printInfo();
//		int auswahl = 6;
//		cin >> auswahl;
//		if (auswahl == 0)AMO.nesterov = AMO.speedup = AMO.verbose = true;
//		if (auswahl == 1)AMO.nesterov = AMO.speedup = AMO.parallelTest = true;
//		if (auswahl == 2)AMO.longstaffschwarz = true;
//		if (auswahl == 3)AMO.andersenbroadie = true;
//		if (auswahl == 4)AMO.semiinf = AMO.verbose = true;
//	}
	bool wieder = false;
	for (int i = 0; i < argc; ++i) {
		string arg = args[i];
		bool geaendert = false;
//		if(! arg.compare("-testing"))               {geaendert=true;AMO.testing=true;}
		if (!arg.compare("-w")) {
			geaendert = true;
			wieder = true;
		}
		if (!arg.compare("-10")) {
			geaendert = true;
			runden = 10;
		}
		if (!arg.compare("-20")) {
			geaendert = true;
			runden = 20;
		}
		if (!arg.compare("-50")) {
			geaendert = true;
			runden = 50;
		}
		if (!arg.compare("-100")) {
			geaendert = true;
			runden = 100;
		}
		if (!arg.compare("-1000")) {
			geaendert = true;
			runden = 1000;
		}
//		if(! arg.compare("-extremeTesting"))        {geaendert=true;AMO.extremTest=true;}
//		if(! arg.compare("-parallelTesting"))       {geaendert=true;AMO.parallelTest=true;}
		if (!arg.compare("-Nesterov")) {
			geaendert = true;
			AMO.nesterov = true;
		}
//		if(! arg.compare("-verbose"))               {geaendert=true;AMO.verbose=true;}
//		if(! arg.compare("-loadAlphas"))            {geaendert=true;AMO.loadAlphas=true;}
		if (!arg.compare("-zehnmal")) {
			geaendert = true;
			runden = 10;
		};
		if (!arg.compare("-fuenfzigmal")) {
			geaendert = true;
			runden = 50;
		};
		if (!arg.compare("-hundertmal")) {
			geaendert = true;
			runden = 100;
		};
//		if(! arg.compare("-speedup"))               {geaendert=true;AMO.speedup=true;}
		if (!arg.compare("-BFGS")) {
			geaendert = true;
			AMO.bfgs = true;
		}
		if ((!arg.compare("-LongstaffSchwartz")) || (!arg.compare("-LSM"))) {
			geaendert = true;
			AMO.longstaffschwarz = true;
		}
//		if(! arg.compare("-AndersenBroadie"))       {geaendert=true;AMO.andersenbroadie=true;}
//		if(! arg.compare("-antithetics"))           {geaendert=true;AMO.antithetics=true;}
//		if(! arg.compare("-semi"))                  {geaendert=true;AMO.semiinf=true;}
		if (i > 0 && !geaendert) {
			printf("Unverständliche Parameter!\n");
			return 0;
		}
	}
	if (wieder) {
		printf("Anzahl der Wiederholungen\n");
		cin >> runden;
	}

	for (int i = 0; i < runden; ++i) {
		if (AMO.longstaffschwarz)
			AMO.LongstaffSchwartz();
		if (AMO.bfgs)
			AMO.BFGS_aufrufen();
		if (AMO.nesterov)
			AMO.Nesterov_aufrufen();
	}
}
