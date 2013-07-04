#include <stdio.h>
#include "AmericanOption.h"
#include <cstring>
#include <string.h>
#include <iostream>

using namespace std;

void printInfo()
{
	printf("\nInfo - Folgende Argumente sind möglich: (Zahl eingeben)\n");
	printf("0: ./AmericanOption -Nesterov -speedup -verbose\n");
	printf("1: ./AmericanOption -Nesterov -speedup -parallelTesting\n");
	printf("2: ./AmericanOption -LongstaffSchwartz\n");
	printf("3: ./AmericanOption -AndersenBroadie\n");
	printf("4: ./AmericanOption -semi -verbose\n");
}

void Test(){
	double daten[5]={10,2,50,4,30};
	int* index=BubbleSort(daten,5);
	for(int i=0;i<5;++i)
		printf("%d\n",index[i]);
}

int main( int argc, char* args[]) {
	int runden=1;
	AmericanOption AMO;

	if (argc == 1) {
		printInfo();
		int auswahl = 6;
		cin >> auswahl;
		if (auswahl == 0)AMO.nesterov = AMO.speedup = AMO.verbose = true;
		if (auswahl == 1)AMO.nesterov = AMO.speedup = AMO.parallelTest = true;
		if (auswahl == 2)AMO.longstaffschwarz = true;
		if (auswahl == 3)AMO.andersenbroadie = true;
		if (auswahl == 4)AMO.semiinf = AMO.verbose = true;
	}
	bool wieder=false;
	for(int i=0;i<argc;++i)
	{
		string arg=args[i];
		bool geaendert=false;
		if(! arg.compare("-testing"))               {geaendert=true;AMO.testing=true;}
		if(! arg.compare("-w"))                     {geaendert=true;wieder=true;}
		if(! arg.compare("-10"))                    {geaendert=true;runden=10;}
		if(! arg.compare("-20"))                    {geaendert=true;runden=20;}
		if(! arg.compare("-50"))                    {geaendert=true;runden=50;}
		if(! arg.compare("-100"))                   {geaendert=true;runden=100;}
		if(! arg.compare("-1000"))                  {geaendert=true;runden=1000;}
		if(! arg.compare("-extremeTesting"))        {geaendert=true;AMO.extremTest=true;}
		if(! arg.compare("-parallelTesting"))       {geaendert=true;AMO.parallelTest=true;}
		if(! arg.compare("-Nesterov"))              {geaendert=true;AMO.nesterov=true;}
		if(! arg.compare("-verbose"))               {geaendert=true;AMO.verbose=true;}
		if(! arg.compare("-loadAlphas"))            {geaendert=true;AMO.loadAlphas=true;}
		if(! arg.compare("-zehnmal"))               {geaendert=true;runden=10;};
		if(! arg.compare("-fuenfzigmal"))           {geaendert=true;runden=50;};
		if(! arg.compare("-hundertmal"))            {geaendert=true;runden=100;};
		if(! arg.compare("-speedup"))               {geaendert=true;AMO.speedup=true;}
		if(! arg.compare("-BFGS"))                  {geaendert=true;AMO.bfgs=true;}
		if((! arg.compare("-LongstaffSchwartz"))
				|| (! arg.compare("-LSM")))         {geaendert=true;AMO.longstaffschwarz=true;}
		if(! arg.compare("-AndersenBroadie"))       {geaendert=true;AMO.andersenbroadie=true;}
		if(! arg.compare("-antithetics"))           {geaendert=true;AMO.antithetics=true;}
		if(! arg.compare("-semi"))                  {geaendert=true;AMO.semiinf=true;}
		if(i>0 && !geaendert){printf("Unverständliche Parameter!\n");return 0;}
	}
	if (wieder) {
		printf("Anzahl der Wiederholungen\n");
		cin >> runden;
	}

	for(int i=0;i<runden;++i){
		if (AMO.longstaffschwarz)AMO.LongstaffSchwartz();
		if (AMO.andersenbroadie)AMO.AndersenBroadie();
		if (AMO.bfgs)AMO.BFGS();
		if (AMO.nesterov)AMO.Nesterov();
		if (AMO.testing)AMO.BFGS_testing(50000.);
		if (AMO.parallelTest)AMO.BFGS_parallelTesting(50000.);
		if (AMO.extremTest)AMO.BFGS_extremeTesting(25, 50000.);
		if (AMO.semiinf)AMO.semi();
	}
}
