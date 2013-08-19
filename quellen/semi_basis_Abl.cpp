//
//#include "AmericanOption.h"
//
//using namespace std;
//
//double* AmericanOption::semi_BasisfunktionenAbl(int zeit, int j, double* x) {
//	if (D == 1)	return semi_Basisfunktionen1DAbl(zeit,j,x);
//	if (D == 2)	return semi_Basisfunktionen2DAbl(zeit,j,x);
//	if (D > 2) 	return semi_BasisfunktionenHigherDAbl(zeit,j,x);
//	printf("Error596");
//	return NULL;
//}
//
//double* AmericanOption::semi_Basisfunktionen1DAbl(int zeit, int j, double* x)
//{
//	if (j == 5)
//		return payoffAbl(x, zeit);
//
//	double* grad=DoubleFeld(D);
//	if (j < 6) {
//		if (j == 0)
//			grad[0]=0;//1;
//		if (j == 1)
//			grad[0]=1;//x[0];
//		if (j == 2)
//			grad[0]=2.*x[0];//x[0] * x[0];
//		if (j == 3)
//			grad[0]=3.* x[0] * x[0];
//		if (j == 4)
//			grad[0]=4. * x[0] * x[0] * x[0]; //Polynom 5. Grades bringt nur Probleme
//		return grad;
//	}
//	j -= 6;
//
//	if (option == MIN_PUT){
//		grad[0]=max(Strike * 4. * (double) (j) / (double) (Mphi - 6) - x[0], 0)==0?0:-1.;
//	}else
//		grad[0]=max(x[0] - Strike * (double) (j) / (double) (Mphi - 6), 0)==0?0:1.;
//	return grad;
//}
//
//double* AmericanOption::semi_Basisfunktionen2DAbl(int zeit, int j, double* x) {
//	//	int reihe[D];
//	//	for(int jj=0;jj<D;++jj)
//	//		reihe[jj]=jj;
//	//	BubbleSort(x,reihe,D);
//	//
//	//	double y[D];
//	//	for (int d = 0; d < D; ++d) {
//	//		y[d] = x[d];
//	////		y[d] -= 100;
//	////		y[d] /= 100;
//	//	}
//
//	if (j == 6)
//		return payoffAbl(x, zeit);
//
//	if(j<7){
//		double* grad=DoubleFeld(D);
//if(j==0)
//	grad[0]=0;
//		if (j == 1)
//			grad[0]=1;
//		if (j == 2)
//			grad[1]=1;
//		if (j == 3)
//		{
//			grad[1]= x[0];
//			grad[0]=x[1];
//		}
//		if (j == 4)
//			grad[0]=2. * x[0];
//		if (j == 5)
//			grad[1]=2. * x[1];
//		return grad;
//	}
//	j -= 7;
//
//
//
//
//	if (j <1000)
//	{
//		double* grad=DoubleFeld(D);
//		grad[1]=max((x[0] + x[1]) -  j / 1000. *4* X0[0], 0)==0?0:1;
//		grad[0]=grad[1];
//		return grad;
//	}
//	j -= 1000;
//
//if (j <1000){
//		double* grad=DoubleFeld(D);
//		grad[0]=max(fabs(x[0] - x[1])-  j / 1000. * X0[0], 0)==0?0:1*pow((-1),(x[0]<x[1]));
//		grad[1]=max(fabs(x[0] - x[1])-  j / 1000. * X0[0], 0)==0?0:1*pow((-1),(x[0]>x[1]));
//		return grad;
//	}
//	j -= 1000;
//
//	if (j < 1000) {
//		double xx[2];
//		xx[0] = j / 1000. * 2. * x[0];
//		xx[1] = j / 1000. * 2. * x[1];
//		double* gr= payoffAbl(xx, zeit);
//		gr[0]*=j / 1000. * 2.;
//		gr[1]*=j / 1000. * 2.;
//				return gr;
//	}
//	j -= 1000;
//
//	printf("Error 651 \n");return NULL;
//}
//
//double* AmericanOption::semi_BasisfunktionenHigherDAbl(int zeit, int j, double* x) {
//double* grad=DoubleFeld(D);
//return grad;
///*
//int reihe[D];
//	for(int jj=0;jj<D;++jj)
//		reihe[jj]=jj;
//	BubbleSort(x,reihe,D);
//
//	if(j<1)return grad;
//	j-=1;
//
//	if(j<3){
//		grad[0]=(j+2.)*pow(x[reihe[0]],j+2);
//		return grad;
//	}
//	j-=3;
//
//	if(j<D*2 )
//	{
//		int a=j%D;
//		int b=(j-a)/D;
//		//if(verbose)printf("%d:asset auf platz %d hoch %d\n",j,a,b);
//		grad[reihe[a]]= (b+1)*pow(x[reihe[a]],b);
//		return grad;
//	}
//	j-=D*2;
//
//	if(j<D-1){
//		grad[reihe[j]]=x[reihe[j+1]];
//		grad[reihe[j+1]]=x[reihe[1]];
//		return grad;
//	}
//	j-=D-1;
//
//	if(j<1)
//	{
//		double product=1;
//		for(int jj=0;jj<D;++jj)
//			product*=x[reihe[jj]];
//		for(int d=0;d<D;++d)
//			grad[reihe[d]]=product/x[reihe[d]];
//		return grad;
//	}
//	j-=1;
//
//	if(j<1){
//		grad[reihe[0]]=2.*x[reihe[0]]*x[reihe[1]];
//		grad[reihe[1]]=x[reihe[0]]*x[reihe[0]];return grad;
//	}
//	j-=1;
//
//	if(j<1){
//		grad[reihe[1]]=2.*x[reihe[1]]*x[reihe[0]];
//		grad[reihe[0]]=x[reihe[1]]*x[reihe[1]];
//		return grad;
//	}
//	j-=1;
//
//	if (j < 1000) {
//		double diff = x[reihe[0]]-x[reihe[1]];
//		return max(diff - (double) (j) / 1000. * X0[0], 0);
//	}
//	j -= 1000;
//
//	if (j < 1000) {
//		double beide=x[reihe[0]]+x[reihe[1]];
//		return max(beide - (double) (j) / 1000. * 20. * X0[0], 0);
//	}
//	j -= 1000;
//
//	if (j < 1000) {
//		double summe = 0;
//		for (int d = 0; d < D; ++d)
//			summe += x[d];
//		return max(summe - (double) (j) / 1000. * 10. * (double) D * X0[0], 0);
//	}
//	j -= 1000;
//
//	if (j < 1000) {
//		double xx[D];
//		for (int d = 0; d < D; ++d)
//			xx[d] = (double) j / 1000. * 10. * D * x[d];
//		return payoffAbl(xx, zeit);
//	}
//	j -= 1000;
//
//	if (j < 1000) {
//		double xx[D];
//		for (int d = 0; d < D; ++d){
//			xx[d] = x[d];
//			if(j%D==d)
//				xx[d] *=(double) j / 1000. * 10.;
//		}
//		return payoffAbl(xx, zeit);
//	}
//	j -= 1000;
//
//	if (j < 1000) {
//		double xx[D];
//		for (int d = 0; d < D; ++d){
//			xx[d] = x[d];
//			if(j%D!=d)
//				xx[d] *=(double) j / 1000. * 10.;
//		}
//		return payoff(xx, zeit);
//	}
//	j -= 1000;
//
//	if (j < 1000) {  //NEU
//		double xx[D];
//		for (int d = 0; d < D; ++d){
//			xx[d] = x[d];
//			if(j%D==reihe[d])
//				xx[d] *=(double) j / 1000. * 10.;
//		}
//		return payoff(xx, zeit);
//	}
//	j -= 1000;
//*/
//	printf("Error 652 \n");return NULL;
//}
//
