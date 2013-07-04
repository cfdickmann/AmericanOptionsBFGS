
#include "AmericanOption.h"

using namespace std; 

double AmericanOption::semi_Basisfunktionen(int zeit, int j, double* x) {
	if (D == 1)	return semi_Basisfunktionen1D(zeit,j,x);
	if (D == 2)	return semi_Basisfunktionen2D(zeit,j,x);
	if (D > 2) 	return semi_BasisfunktionenHigherD(zeit,j,x);
	printf("Error598");
	return -1;
}

double AmericanOption::semi_Basisfunktionen1D(int zeit, int j, double* x)
{
	if (j < 6) {
		if (j == 0)
			return 1;
		if (j == 1)
			return x[0];
		if (j == 2)
			return x[0] * x[0];
		if (j == 3)
			return x[0] * x[0] * x[0];
		if (j == 4)
			return x[0] * x[0] * x[0] * x[0]; //Polynom 5. Grades bringt nur Probleme
		if (j == 5)
			return payoff(x, zeit);
	}
	j -= 6;

	if (option == MIN_PUT)
		return max(Strike * 4. * (double) (j) / (double) (Mphi - 6) - x[0], 0);
	else
		return max(x[0] - Strike * (double) (j) / (double) (Mphi - 6), 0);
}

double AmericanOption::semi_Basisfunktionen2D(int zeit, int j, double* x) {
//	int reihe[D];
//	for(int jj=0;jj<D;++jj)
//		reihe[jj]=jj;
//	BubbleSort(x,reihe,D);
//
//	double y[D];
//	for (int d = 0; d < D; ++d) {
//		y[d] = x[d];
////		y[d] -= 100;
////		y[d] /= 100;
//	}
	if (j == 0)
		return 1;
	if (j == 1)
		return x[0];
	if (j == 2)
		return x[1];
	if (j == 3)
		return x[0] * x[1];
	if (j == 4)
		return x[0] * x[0];
	if (j == 5)
		return x[1] * x[1];
	if (j == 6)
		return payoff(x, zeit);
	j -= 7;

	if (j <1000){
		return max((x[0] + x[1]) -  j / 1000. *4* X0[0], 0);
	}
	j -= 1000;

if (j <1000){
		return max(fabs(x[0] - x[1])-  j / 1000. * X0[0], 0);//Achtung fabs(x[0] - x[1])
	}
	j -= 1000;

	if (j < 1000) {
		double xx[2];
		xx[0] = j / 1000. * 2. * x[0];
		xx[1] = j / 1000. * 2. * x[1];
		return payoff(xx, zeit);
	}
	j -= 1000;





	printf("Error 653 \n");return -1;
}

double AmericanOption::semi_BasisfunktionenHigherD(int zeit, int j, double* x) {
//	double y[D];
//	for (int d = 0; d < D; ++d) {
//		y[d] = x[d];
//		y[d] -= Strike;
//		y[d] /= Strike;
//	}

	int reihe[D];
	for(int jj=0;jj<D;++jj)
		reihe[jj]=jj;
	BubbleSort(x,reihe,D);

	if(j<1)return 1;
	j-=1;

	if(j<3)
		return pow(x[reihe[0]],j+3);
	j-=3;

	if(j<D*2 )
	{
		int a=j%D;
		int b=(j-a)/D;
		//if(verbose)printf("%d:asset auf platz %d hoch %d\n",j,a,b);
		return pow(x[reihe[a]],b+1);
	}
	j-=D*2;

	if(j<D-1)
		return x[reihe[j]]*x[reihe[j+1]];
	j-=D-1;

	if(j<1)
	{
		double product=1;
		for(int jj=0;jj<D;++jj)
			product*=x[reihe[jj]];
		return product;
	}
	j-=1;

	if(j<1)return x[reihe[0]]*x[reihe[0]]*x[reihe[1]];
	j-=1;

	if(j<1)return x[reihe[1]]*x[reihe[1]]*x[reihe[0]];
	j-=1;

	if (j < 1000) {
		double diff = x[reihe[0]]-x[reihe[1]];
		return max(diff - (double) (j) / 1000. * X0[0], 0);
	}
	j -= 1000;

	if (j < 1000) {
		double beide=x[reihe[0]]+x[reihe[1]];
		return max(beide - (double) (j) / 1000. * 20. * X0[0], 0);
	}
	j -= 1000;

	if (j < 1000) {
		double summe = 0;
		for (int d = 0; d < D; ++d)
			summe += x[d];
		return max(summe - (double) (j) / 1000. * 10. * (double) D * X0[0], 0);
	}
	j -= 1000;

	if (j < 1000) {
		double xx[D];
		for (int d = 0; d < D; ++d)
			xx[d] = (double) j / 1000. * 10. * D * x[d];
		return payoff(xx, zeit);
	}
	j -= 1000;

	if (j < 1000) {
		double xx[D];
		for (int d = 0; d < D; ++d){
			xx[d] = x[d];
			if(j%D==d)
				xx[d] *=(double) j / 1000. * 10.;
		}
		return payoff(xx, zeit);
	}
	j -= 1000;

	if (j < 1000) {
		double xx[D];
		for (int d = 0; d < D; ++d){
			xx[d] = x[d];
			if(j%D!=d)
				xx[d] *=(double) j / 1000. * 10.;
		}
		return payoff(xx, zeit);
	}
	j -= 1000;

	if (j < 1000) {  //NEU
		double xx[D];
		for (int d = 0; d < D; ++d){
			xx[d] = x[d];
			if(j%D==reihe[d])
				xx[d] *=(double) j / 1000. * 10.;
		}
		return payoff(xx, zeit);
	}
	j -= 1000;

	printf("Error 654 \n");return -1;
}

