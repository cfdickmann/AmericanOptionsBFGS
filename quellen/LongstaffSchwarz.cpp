#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "MTRand.h"
#include "math.h"
#include "stdlib.h"
#include "AmericanOption.h"
using namespace std;

AmericanOption* zeiger;

void* DELEGATE_LSM_mittelwert(void* data) {
	zeiger->LSM_mittelwert(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void AmericanOption::LongstaffSchwartz() { //TODO
	LSM_Mtraining=1000000;
	LSM_Mtesting =10000000;

	LSM_setting();

	X=DoubleFeld(M,N,D);
	V=DoubleFeld(M,N);
	B =DoubleFeld(Mphi,Mphi,Threadanzahl);
	BV=DoubleFeld(Mphi,Threadanzahl);
	betas=DoubleFeld(N,Mphi);

	//i
	MT.seed(time(NULL)+getpid());
	srand(time(NULL));
	if(verbose)printf("Pfade erstellen\n");
	RNG generator;
	for (int m = 0; m < M; ++m)
		Pfadgenerieren(X[m],0,X0,&generator);

	//ii
	for(int m=0;m<M;++m)  // letzter zeitschritt
		V[m][N-1]=payoff(X[m][N-1],N-1);

	//iii
	for (int lauf = N - 2; lauf >= 1; --lauf) {
		LSlauf=lauf;
		printf("Schritt %d - ",lauf);


		double** P=DoubleFeld(Mphi,Mphi);
		double PV[Mphi];

		if(verbose)printf("Threads starten\n");
		pthread_t threads[Threadanzahl];
		for(int  i = 0; i < Threadanzahl; i++ ) {
			int* z=(int*)malloc(sizeof(int));
			z[0]=i;
			pthread_create( &threads[i], NULL, DELEGATE_LSM_mittelwert, (void*)z );
			//                       Argument für thread_function ---^
		}

		if(verbose)printf("auf Treads warten\n");
		for(int  i = 0; i < Threadanzahl; i++ )
			pthread_join(threads[i], NULL);

		if(verbose)printf("Ergebnisse von Threads zusammenfuegen\n");
		for (int r = 0; r < Mphi; ++r) {
			PV[r]=0;
			for(int t=0;t<Threadanzahl;++t)
				PV[r]+=BV[r][t]*M/(M);
		}

		if(verbose)printf("Ergebnisse von Threads zusammenfuegen\n");
		for (int r = 0; r < Mphi; ++r)
			for (int q = 0; q < Mphi; ++q) {
				P[r][q]=0;
				for(int t=0;t<Threadanzahl;++t)
					P[r][q]+=B[r][q][t]*M/(M);
			}

		if(verbose)printf("ausrechnen\n");
		betas[lauf]=LGSloesen((double**)P,(double*)PV,Mphi);

		for (int m = 0; m < M; ++m)
			V[m][LSlauf]=max(   (double)payoff(X[m][LSlauf],LSlauf)   ,   (double)LSM_C_estimated(X[m][LSlauf],LSlauf)   );

		if(verbose)
			for(int m=0;m<Mphi;++m)
				printf("%f ,\n",betas[lauf][m]);
	}
	printf("\n");

	//High-Estimation
	double h=0;
	for (int m = 0; m < M; ++m)
		h+=V[m][1]/double(M);
	printf("Regression: %f\n",h);

	// geschaetzte C_i speichern
	fstream f;
	f.open("LSMkoeff.dat", ios::out);   //Achtung, beim ersten und letzten Zeitpunkt Nullen
	for(int n=0;n<N;++n)
		for(int m=0;m<Mphi;++m){
			f<<betas[n][m];//f<<floor(betas[n][m]*100000.)/100000.;
			if(m<Mphi-1 || n<N-1)f<<endl;
		}
	f.close();
	printf("(Koeff der C_i in LSMkoeff.dat gespeichert)\n");

	ErgebnisAnhaengen(h,(char*)"ergebnisse_LSM_high.txt"); // high biased ergebnis ausgeben

	//Low-Estimation
	{int ergebnispipe[Threadanzahl][2];
	for(int z=0;z<Threadanzahl;++z)
		pipe(ergebnispipe[z]);

	int pid;
	for(int f=0;f<Threadanzahl;++f)
	{
		pid = fork();
		if (pid == 0)
		{
			srand(f);
			MT.seed(getpid()+f+time(NULL));
			double e=0;
			double** X=DoubleFeld(N,D);
			double** wdiff=DoubleFeld(N,D);
			int durchlaeufe=LSM_Mtesting/Threadanzahl; //10000000/Threadanzahl
			for (int k = 0; k < durchlaeufe; ++k){
				RNG generator;
				for(int n=0;n<N;++n)
					for(int j=0;j<D;++j)
						wdiff[n][j]=sqrt(dt)*generator.nextGaussian();
				Pfadgenerieren(X,wdiff);
				double erg=payoff(X[N - 1],N - 1 );
				for (int lau = 1; lau<N-1; ++lau)
					if (LSM_C_estimated(X[lau],lau)<= payoff( X[lau],lau)){
						erg= payoff(X[lau],lau);break;
					}
				e += erg / (double)durchlaeufe;
			}
			InPipeSchreiben(ergebnispipe[f],e);
			exit(0);
		}
	}

	double erg=0;
	for(int f=0;f<Threadanzahl;++f)
		erg+=AusPipeLesen(ergebnispipe[f])/(double)(Threadanzahl);

	printf("low bound: %f\n\n",erg);
	ErgebnisAnhaengen(erg,(char*)"ergebnisse_LSM_low.txt");
	}
}

double AmericanOption::LSM_C_estimated(double* x, int time) {
	double erg = 0;
	for (int k = 0; k < Mphi; ++k)
		erg = erg + betas[time][k] * LSM_phi(x, k, time);
	return erg;
}

double AmericanOption::LSM_phi(double* x, int j, int time) {
	if(D==1){
		if (j == 0)
			return 1;
		if (j == 1)
			return x[0];
		if (j == 2)
			return x[0] * x[0];
		if (j == 3)
			return x[0] * x[0] * x[0];
		if (j == 4)
			return x[0] * x[0] * x[0]* x[0];
		if (j == 5)
			return payoff(x,time);
	}

	if(D==2){
		if (j == 0)
			return 1;
		if (j == 1)
			return x[0];
		if (j == 2)
			return x[0] * x[1];
		if (j == 3)
			return x[1];
		if (j == 4)
			return x[0] * x[0];
		if (j == 5)
			return x[1] * x[1] ;
		if (j == 6)
			return payoff(x,time);
	}

	if(mlsm){
		if(D>2){
			int reihe[D];
			for(int jj=0;jj<D;++jj)
				reihe[jj]=jj;
			BubbleSort(x,reihe,D);

			if(j<LSM_K0)return 1;

			if(j<LSM_K0 + LSM_K1 && j>=LSM_K0)
				return pow(x[reihe[0]],j-LSM_K0+3);

			if(j<LSM_K0 + LSM_K1 + LSM_K2 && j>=LSM_K0 + LSM_K1)
			{
				int a=(j-LSM_K0-LSM_K1)%D;
				int b=((j-LSM_K0-LSM_K1)-a)/D;
				//if(verbose)printf("gemischt %d,%d\n",a,b);
				return pow(x[reihe[a]],b+1);
			}

			if(j<LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 && j>=LSM_K0 + LSM_K1 + LSM_K2)
				return x[reihe[j-LSM_K1-LSM_K2-LSM_K0]]*x[reihe[j-LSM_K1-LSM_K2-LSM_K0+1]];

			if(j<LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4 && j>=LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3)
			{
				double product=1;
				for(int jj=0;jj<D;++jj)
					product*=x[reihe[jj]];
				return product;
			}
		}
	}
	else{ // not mlsm
		if(D>2){
			int reihe[D];
			for(int jj=0;jj<D;++jj)
				reihe[jj]=jj;
			BubbleSort(x,reihe,D);

			if(j<LSM_K0)return 1;

			if(j<LSM_K0 + LSM_K1 && j>=LSM_K0)
			{
				if(j-LSM_K0==0)return x[reihe[0]];
				if(j-LSM_K0==1)return x[reihe[1]];
				if(j-LSM_K0==2)return pow(x[reihe[0]],2);
				if(j-LSM_K0==3)return pow(x[reihe[1]],2);
				if(j-LSM_K0==4)return pow(x[reihe[0]],3);
				if(j-LSM_K0==5)return pow(x[reihe[1]],3);
				if(j-LSM_K0==6)return pow(x[reihe[0]],4);
				if(j-LSM_K0==7)return pow(x[reihe[1]],4);
				if(j-LSM_K0==8)return pow(x[reihe[0]],5);
				if(j-LSM_K0==9)return pow(x[reihe[1]],5);
			}

			if(j<LSM_K0 + LSM_K1 + LSM_K2 && j>=LSM_K0 + LSM_K1)
			{
				if(j-LSM_K0-LSM_K1==0)return pow(x[reihe[0]],1)*pow(x[reihe[1]],1);
				if(j-LSM_K0-LSM_K1==1)return pow(x[reihe[0]],1)*pow(x[reihe[1]],2);
				if(j-LSM_K0-LSM_K1==2)return pow(x[reihe[0]],2)*pow(x[reihe[1]],1);
				if(j-LSM_K0-LSM_K1==3)return pow(x[reihe[0]],2)*pow(x[reihe[1]],2);
				if(j-LSM_K0-LSM_K1==4)return pow(x[reihe[0]],3)*pow(x[reihe[1]],1);
				if(j-LSM_K0-LSM_K1==5)return pow(x[reihe[0]],3)*pow(x[reihe[1]],2);
				if(j-LSM_K0-LSM_K1==6)return pow(x[reihe[0]],2)*pow(x[reihe[1]],3);
				if(j-LSM_K0-LSM_K1==7)return pow(x[reihe[0]],1)*pow(x[reihe[1]],3);
				if(j-LSM_K0-LSM_K1==8)return pow(x[reihe[0]],4)*pow(x[reihe[1]],1);
				if(j-LSM_K0-LSM_K1==9)return pow(x[reihe[0]],1)*pow(x[reihe[1]],4);
			}
		}
	}

	printf("ERROR56 %d\n",j);
	return 0;
}

double GBM(double drift, double sigma, double t, double normal) {
	return exp((drift - 0.5 * sigma * sigma) * t + sqrt(t) * sigma * normal);
}

void AmericanOption::LSM_mittelwert(int threadnummer)
{
	int mAnfang=M*threadnummer/Threadanzahl;
	int mEnde=M*(threadnummer+1)/Threadanzahl;

	if(verbose)printf("Thread %d: von %d bis %d\n",threadnummer, mAnfang,mEnde);

	for (int r = 0; r < Mphi; ++r) {
		for (int q = 0; q < Mphi; ++q) {
			double erg = 0;
			for (int m = mAnfang; m < mEnde; ++m)
				//				if(1 ==0 && LSlauf!=0 && LSlauf!=1 && LSlauf!=N-1 && LSlauf!=N-2 && Cestimated(X[m][LSlauf+1],LSlauf+1)<= payoff( X[m][LSlauf+1],LSlauf+1))
				//					ausgesetzt[threadnummer]++; else
				erg = erg + LSM_phi(X[m][LSlauf], r, LSlauf) * LSM_phi(X[m][LSlauf], q, LSlauf);
			B[r][q][threadnummer] = erg / (double) M;
		}

		double erg2 = 0;
		for (int m = mAnfang; m < mEnde; ++m)
			//			if(1==0 && LSlauf!=0 && LSlauf!=1 && LSlauf!=N-1 && LSlauf!=N-2 && Cestimated(X[m][LSlauf+1],LSlauf+1)<= payoff( X[m][LSlauf+1],LSlauf+1))
			//				ausgesetzt[threadnummer]++; else
			erg2 = erg2 + LSM_phi(X[m][LSlauf], r, LSlauf) *V[m][LSlauf+1];
		BV[r][threadnummer] = erg2 / (double) M;
	}
}

void AmericanOption::LSM_setting(){
	zeiger=this;  // damit threada auf diese objekt zugreifen koennen
	Daten(); //Problemdaten laden

	if(D==1)Mphi=6;
	if(D==2)Mphi=7;
	if(D>2){
		if(mlsm)   //Anzahl der Basisfunktionen initialisieren
		{
			LSM_K0=1;//Konstante
			LSM_K1=3;// 3,4,5 polynom des teuersten assets
			LSM_K2=D*2;//Polynome zweiter ordnung in allen einzelnen assets
			if(D>2)LSM_K3=(D-1);//Produkte von Verfolgern
			else LSM_K3=0;
			LSM_K4=1;//Produkte aller assets
		}else{
			LSM_K0=1;//Konstante
			LSM_K1=10;
			LSM_K2=10;//Polynome zweiter ordnung in allen einzelnen assets
			LSM_K3=0;
			LSM_K4=0;//Produkte aller assets
		}
		Mphi=LSM_K0+LSM_K1+LSM_K2+LSM_K3+LSM_K4;
	}

	//Basisfunktionen testen
	//	double tt[Mphi];
	//	for(int i=0;i<Mphi;++i)
	//		tt[i]=i%D+2;

	//	for(int t=0;t<Mphi;++t){
	//		if( t==K_0 || t==K_1+K_0 || t==K_1+K_2+K_0 || t==K_1+K_2+K_0+K_3 || t==K_1+K_2+K_0+K_3+K_4)printf("\n");
	//		printf("t=%d, %f\n",t,phi(tt,t,1));
	//	}

	N=Testing_Dates;  // fuer longstaff schwartz keine hoehere genauigkeit noetig, da explizite formel fuer pfade verwendet wurde
	dt=T/(double(N-1));
	M=LSM_Mtraining;
}

