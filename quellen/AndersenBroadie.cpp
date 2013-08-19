//#include <iostream>
//#include <math.h>
//#include <stdlib.h>
//#include "MTRand.h"
//#include "math.h"
//#include "stdlib.h"
//#include "AmericanOption.h"
//
//using namespace std;
//
//AmericanOption* zeiger4;
//
//void AmericanOption::AndersenBroadie()	//Andersen-Broadie Multilevel //TODO
//{
//	int L=4;   //Anzahl der Level
//	double n0=200; //Anzahl der Pfade
//	double k0=300; // Anzahl der subsimulaions
//	double faktor=2; // Das ist KAPPA !
//
//	LSM_setting();
//
//	betas=betasAusDateiLaden();
//
//	double Gesamtergebnis=0;
//
//	int k[L];
//	int n[L];
//	k[0]=k0;
//	n[0]=n0;
//
//	for(int l=0;l<L;++l) //Anfang des Levels
//	{
//		n[l]=(int)ceil(n0*pow(faktor,-(double)l)/Threadanzahl);  // neue n und k setzen
//		k[l]=(int)(k0*pow(faktor,(double)l));
//
//		// Pipes zur Uebergabe von ergebnissen an das Hauptprogramm definieren
//		int ergPipe[Threadanzahl][2];			//zur uebergabe der differenz
//		int mittelwertPipe[Threadanzahl][2];   	// zur uebergabe der mittelwerte
//		int quadratsummePipe[Threadanzahl][2]; 	// zur uebergabe von summe (x_i^2)
//
//		for(int z=0;z<Threadanzahl;++z){// pipes erzeugen
//			pipe(mittelwertPipe[z]);
//			pipe(quadratsummePipe[z]);
//			pipe(ergPipe[z]);
//		}
//
//		int pid;
//		for(int f=0;f<Threadanzahl;++f)
//		{
//			pid = fork();  // Prozess duplizieren
//			if (pid == 0)
//			{
//				srand(getpid()+f+time(NULL));  //seed von Zufallsgeneratoren setzen, damit nicht jeder instanz das gleiche
//				MT.seed(getpid()+f+time(NULL)); // Ergebnis berechnet
//
//				double ** x=DoubleFeld(N,D);
//
//				double ergErgebnis;
//				double ergMittelwert=0;
//				double ergQuadratsumme=0;
//				for (int nn = 0; nn < n[l]; ++nn)  // nn-mal simulieren
//				{
//				RNG generator;
//					Pfadgenerieren(x,0,X0,&generator);
//					double ABminus=0;
//					if(l>0)ABminus=AndersenBroadieEinzel(x,k[l-1]); //im ersten Schritt keine Subtraktion
//					double ABplus=AndersenBroadieEinzel(x,k[l]);
//					ergErgebnis+=  (ABplus-ABminus)/(double)(n[l]);
//					ergMittelwert+=(ABplus-ABminus)/(double)(n[l]);
//					ergQuadratsumme+=pow(ABplus-ABminus,2)/(double)(n[l]);
//				}
//				//printf("zwischen %d: ergErgebnis: %f\n",l,ergErgebnis);
//				InPipeSchreiben(ergPipe[f],ergErgebnis);                // Ergebnisse an Hauptinstanz des Programms uebermitteln
//				InPipeSchreiben(mittelwertPipe[f],ergMittelwert);
//				InPipeSchreiben(quadratsummePipe[f],ergQuadratsumme);
//				exit(0);  // Duplizierte Instanzen beenden
//			}
//		}
//
//		double ergMittelwerte=0;
//		double ergQuadratsummen=0;
//		for(int f=0;f<Threadanzahl;++f)  //Ergebnisse aus pipe auslesen
//		{
//			double e=AusPipeLesen(ergPipe[f]);
//                        printf("%f\n",e);
//                        Gesamtergebnis+=e/(double)(Threadanzahl);
//			ergMittelwerte+=e/(double)(Threadanzahl);
//			ergQuadratsummen+=AusPipeLesen(quadratsummePipe[f])/(double)(Threadanzahl);
//		}
//		double var=ergQuadratsummen-pow(ergMittelwerte,2);
//		printf("Level %d,  n_l=%d, k_l=%d: %.2lf (%.5lf)\n", l, Threadanzahl*n[l], k[l],ergMittelwerte,sqrt(var));
//	}
//	printf("AndersenBroadie");
//	if(L>1)printf(" (Multilevel)");
//	printf(": %f\n\n",Gesamtergebnis);
//	ErgebnisAnhaengen(Gesamtergebnis);// Ergebnisse werden in ergebnisse.dat ausgegeben
//}
//
//double** AmericanOption::betasAusDateiLaden() {
//    double** betas = DoubleFeld(N, Mphi);
//    //Koeff aus Datei laden
//    string line;
//    ifstream myfile("LSMkoeff.dat");
//    if (myfile.is_open()) {
//        int k = 0;
//        while ((!myfile.eof()) && k < N * Mphi) {
//            //printf("gelesen\n");
//            getline(myfile, line);
//            char * buffer = new char[line.length()];
//            strcpy(buffer, line.c_str());
//            //printf("%d,%d\n",(k-k%Mphi)/Mphi,k%Mphi);
//            betas[(k - k % Mphi) / Mphi][k % Mphi] = (double) (atof(buffer));
//            k++;
//        }
//        myfile.close();
//    }
//
//    printf("(Koeff aus LSMkoeff.dat geladen)\n");
//    if (verbose)
//        for (int n = 0; n < N; ++n)
//            for (int m = 0; m < Mphi; ++m)
//                printf("%f\n", betas[n][m]);
//    return betas;
//}
//
//double AmericanOption::AndersenBroadieEinzel(double ** x, int nsubpaths) {
//    return AndersenBroadieEinzel(x, nsubpaths, 0);
//}
//
//double AmericanOption::AndersenBroadieEinzel(double ** x, int nsubpaths, int startzeit) {
//    double Delta, Etau[N];
//    double M = 0;
//    double max = 0;
//    double ** xx = DoubleFeld(N, D);
//    for (int n = startzeit; n < N; ++n) {
//        if (verbose)printf("\nZeitschritt %d\n", n);
//        Etau[n] = 0;
//        if (n < N - 1) {
//            int seed = (int)(x[n + 1][0]*100000);
//           // MT.seed(seed);
//            srand(seed);
//        }
//        for (int laufsub = 0; laufsub < nsubpaths; ++laufsub) {
//            if (verbose)printf("Nested simulation no. %d: \n", laufsub);
//            Pfadgenerieren(xx, n, x[n]);
//            if (verbose) for (int i = 0; i < N; i++)printf("%.2lf, ", xx[i][0]);
//            int ex_time = N - 1;
//            for (int testen = N - 1; testen > n; --testen)
//                if (payoff(xx[testen], testen) >= LSM_C_estimated(xx[testen], testen))
//                    ex_time = testen;
//            Etau[n] += payoff(xx[ex_time], ex_time) / (double) (nsubpaths);
//            if (verbose)printf("best payoff %f at %d\n", payoff(xx[ex_time], ex_time), ex_time);
//        }
//        if (n > startzeit) {
//            if (payoff(x[n], n) < LSM_C_estimated(x[n], n))Delta = Etau[n] - Etau[n - 1];
//            else Delta = payoff(x[n], n) - Etau[n - 1];
//        } else
//            Delta = 0;
//        M += Delta;
//        double d = payoff(x[n], n) - M;
//        if (d > max)max = d;
//    }
//    return max;
//}
//
///*
//double AmericanOption::AndersenBroadieEinzelRand(double ** x, int nsubpaths, int startzeit) {
//    double Delta, Etau[N];
//    double M = 0;
//    double max = 0;
//    double ** xx = DoubleFeld(N, D);
//    for (int n = startzeit; n < N; ++n) {
//        if (verbose)printf("\nZeitschritt %d\n", n);
//        Etau[n] = 0;
//        if (n < N - 1) {
//            int seed = x[n + 1][0]*100000;
//            // MT.seed(seed);
//            srand(seed);
//        }
//        for (int laufsub = 0; laufsub < nsubpaths; ++laufsub) {
//            if (verbose)printf("Nested simulation no. %d: \n", laufsub);
//
//            //Pfadgenerieren(xx, n, x[n]);
//            //----------------------------
//            {
//                double** wdiff = DoubleFeld(N, D);
//                double** sprue = DoubleFeld(N, D);
//                for (int nn = 0; nn < N; ++nn)
//                    for (int j = 0; j < D; ++j) {
//                        wdiff[nn][j] = sqrt(dt) * nextGaussian();
//                        int NumberOfJumps = Poisson(lambdaJump * dt);
//                        sprue[nn][j] = 0;
//                        for (int jump = 0; jump < NumberOfJumps; ++jump)
//                            sprue[nn][j] += newSprung();
//                    }
//                Pfadgenerieren(xx, wdiff, sprue, n, x[n]);
//                for (int i = 0; i < N; ++i) {
//                    free(wdiff[i]);
//                    free(sprue[i]);
//                }
//            }
//            //-------------
//            if (verbose) for (int i = 0; i < N; i++)printf("%.2lf, ", xx[i][0]);
//            int ex_time = N - 1;
//            for (int testen = N - 1; testen > n; --testen)
//                if (payoff(xx[testen], testen) >= LSM_C_estimated(xx[testen], testen))
//                    ex_time = testen;
//            Etau[n] += payoff(xx[ex_time], ex_time) / (double) (nsubpaths);
//            if (verbose)printf("best payoff %f at %d\n", payoff(xx[ex_time], ex_time), ex_time);
//        }
//        if (n > startzeit) {
//            if (payoff(x[n], n) < LSM_C_estimated(x[n], n))Delta = Etau[n] - Etau[n - 1];
//            else Delta = payoff(x[n], n) - Etau[n - 1];
//        } else
//            Delta = 0;
//        M += Delta;
//        double d = payoff(x[n], n) - M;
//        if (d > max)max = d;
//    }
//    return max;
//}
//
//
//
//void* DELEGATE_AndersenBraodie5000(void* data) {
//    zeiger4->AndersenBroadie5000(((int*) data)[0]);
//    pthread_exit(NULL);
//    return NULL;
//}
//
//void AmericanOption::AndersenBroadie5000(int threadnummer) {
//    int seed = getpid() + threadnummer + time(NULL);
//    srand(seed);
//    MT.seed(seed);
//
//    printf("nummer: %d\n", threadnummer);
//
//    for (int l = 0; l < AB_Level; ++l) //Anfang des Levels
//    {
//        //    printf("level: %d\n",l);
//        double ** x = DoubleFeld(N, D);
//        double ** wdiff = DoubleFeld(N, D);
//        double ** sprue = DoubleFeld(N, D);
//
//        for (int nn = 0; nn < N; ++nn)
//            for (int d = 0; d < D; ++d)
//                sprue[nn][d] = 0; // keine Sprungprozesse
//
//        double ergMittelwert = 0;
//        double ergQuadratsumme = 0;
//        for (int nn = 0; nn < AB_n[l]; ++nn) // nn-mal simulieren
//        {
//            //  printf("hier: %d\n",nn);
//            for (int nnn = 0; nnn < N; ++nnn) // Fuer jeden Zeitschritt
//                for (int d = 0; d < D; ++d) // Fue jedes asset
//                    if (antithetics && nn % 2 == 1) // jedes 2te mal
//                        wdiff[nnn][d] *= -1;
//                    else //Antithetics!
//                        wdiff[nnn][d] = sqrt(dt) * BoxMuller((double) (random()) / (double) (RAND_MAX), (double) (random()) / (double) (RAND_MAX)); // neue Wdiff
//            Pfadgenerieren(x, wdiff, sprue);
//            double ABminus = 0;
//            if (l > 0)ABminus = AndersenBroadieEinzel(x, AB_k[l - 1]); //im ersten Schritt keine Subtraktion					double ABplus=AndersenBroadieEinzel(x,k[l]);
//            double ABplus = AndersenBroadieEinzel(x, AB_k[l]);
//
//            ergMittelwert += (ABplus - ABminus) / (double) (AB_n[l]);
//            ergQuadratsumme += pow(ABplus - ABminus, 2) / (double) (AB_n[l]);
//        }
//        AB_mittelwerte[testlaufactual][l][threadnummer] = ergMittelwert;
//        AB_quadratsummen[testlaufactual][l][threadnummer] = ergQuadratsumme;
//    }
//
//}
//
//void AmericanOption::AndersenBroadieTest() {
//    zeiger4 = this;
//    AB_Level = 4;
//    double AB_kappa = 2;
//    AB_k = (int*) malloc(sizeof (int) *AB_Level);
//    AB_n = (int*) malloc(sizeof (int) *AB_Level);
//    AB_k[0] = 10;
//    AB_n[0] = 10;
//    int AB_durchlaeufe = 10;
//    for (int l = 1; l < AB_Level; ++l) {
//        AB_n[l] = (int) ceil(AB_n[0] * pow(AB_kappa, -(double) l) / Threadanzahl); // neue n und k setzen
//        AB_k[l] = (int) (AB_k[0] * pow(AB_kappa, (double) l));
//    }
//    AB_mittelwerte = DoubleFeld(AB_durchlaeufe, AB_Level, Threadanzahl);
//    AB_quadratsummen = DoubleFeld(AB_durchlaeufe, AB_Level, Threadanzahl);
//
//    LSM_setting();
//
//    betas = betasAusDateiLaden();
//
//    printf("AndersenBroadie");
//    if (AB_Level > 1)printf(" (Multilevel)");
//    printf("\n");
//
//    for (int l = 0; l < AB_Level; ++l)
//        printf("Level %d,  n_l=%d, k_l=%d \n", l, Threadanzahl * AB_n[l], AB_k[l]);
//
//    for (int testlauf = 0; testlauf < AB_durchlaeufe; ++testlauf) {
//        testlaufactual = testlauf;
//        pthread_t threads[Threadanzahl];
//
//        for (int t = 0; t < Threadanzahl; t++)
//            pthread_create(&threads[t], NULL, DELEGATE_AndersenBraodie5000, (void*) array_machen(t));
//
//        for (int t = 0; t < Threadanzahl; t++)
//            pthread_join(threads[t], NULL);
//
//        for (int l = 0; l < AB_Level; ++l) {
//            double ergMittel = 0;
//            for (int t = 0; t < Threadanzahl; ++t)
//                for (int n = 0; n <= testlaufactual; ++n)
//                    ergMittel += AB_mittelwerte[n][l][t] / (double) (Threadanzahl * (testlaufactual + 1));
//            double ergQuadrat = 0;
//            for (int t = 0; t < Threadanzahl; ++t)
//                for (int n = 0; n <= testlaufactual; ++n)
//                    ergQuadrat += AB_quadratsummen[testlaufactual][l][t] / (double) (Threadanzahl * (testlaufactual + 1));
//
//            printf("L %d,  %.2lf (%.5lf)\n", l, ergMittel, sqrt(ergQuadrat - pow(ergMittel, 2)));
//
//        }
//    }
//
//    double Gesamtergebnis = 0;
//    for (int t = 0; t < Threadanzahl; ++t)
//        for (int n = 0; n < AB_durchlaeufe; ++n)
//            for (int l = 0; l < AB_Level; ++l)
//                Gesamtergebnis += AB_mittelwerte[n][l][t] / (double) (Threadanzahl * n);
//
//    printf(": %f\n\n", Gesamtergebnis);
//    ErgebnisAnhaengen(Gesamtergebnis);
//}
// *
// * */
