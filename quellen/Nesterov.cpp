#include "AmericanOption.h"
using namespace std;

#include "Hilfsmittel.h"

int argMinlauf;
AmericanOption* zeiger2;

//static void D_argMinfunction1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *ptr)
//{
//	zeiger2->argMinfunction1_grad(y,func,grad,NULL);
//}
//
//void AmericanOption::argMinfunction1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *ptr)
//{
//	double erg=0;
//	for(int k=0;k<K;++k)
//		for(int l=0;l<=argMinlauf;++l)
//			erg+=((double)l+1.)/2.*( Dfx[l][k]*(y[k]-x[l][k]));
//
//	for(int l=0;l<=argMinlauf;++l)
//		erg+=((double)l+1.)/2.*fx[l];
//
//	double s=0;
//	for(int k=0;k<K;++k)s+=0.5*L*pow(y[k],2);
//
//	func = erg;
//
//	//Gradient
//	for(int k=0;k<K;++k)
//	{
//		grad[k]=2*L*y[k];
//		for(int l=0;l<=argMinlauf;++l)
//			grad[k]+=((double)l+1.)/2.*Dfx[l][k];
//	}
//}

double* AmericanOption::argminProblem(int lauf){
	//bool restringiert=false;
	//if(restringiert)
	//{
//		argMinlauf=lauf;
//		real_1d_array u;
//		u.setlength(K);for(int k=0;k<K;++k)u[k]=0;
//
//		real_1d_array bndl;
//		bndl.setlength(K);for(int k=0;k<K;++k)bndl[k]=-0.3;
//
//		real_1d_array bndu;
//		bndu.setlength(K);for(int k=0;k<K;++k)bndu[k]=0.3;
//		minbleicstate state;
//		minbleicreport rep;
//
//		double epsg = 0.001;
//		double epsf = 0;
//		double epsx = 0;
//		double epso = 0.001;
//		double epsi = 0.001;
//
//		minbleiccreate(u, state);
//		minbleicsetbc(state, bndl, bndu);
//		minbleicsetinnercond(state, epsg, epsf, epsx);
//		minbleicsetoutercond(state, epso, epsi);
//		alglib::minbleicoptimize(state, D_argMinfunction1_grad);
//		minbleicresults(state, u, rep);
//
//		//printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
//		//printf("%s\n", u.tostring(2).c_str()); // EXPECTED: [-1,1]
//		double* e;
//		e=(double*)(malloc(sizeof(double)*K));
//		for(int k=0;k<K;++k)e[k]=u[k];
//		return e;
	//}else
	{//unrestringiert
		double* erg;
		erg=(double*)malloc(sizeof(double)*K);
		for(int k=0;k<K;++k)
		{
			double S=0;
			for(int l=0;l<=lauf;++l)
				S+=((double)l+1.)/2.*Dfx[l][k];
			erg[k]=-S/Nesterov_L;
		}
		return erg;
	}
}

int TQlauf;

//void AmericanOption::TQfunction1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *ptr)
//{
//	double erg=0;
//	for(int k=0;k<K;++k)erg+=Dfx[TQlauf][k]*(y[k]-x[TQlauf][k]);
//	double s=0;
//	for(int k=0;k<K;++k)s+=pow(x[TQlauf][k]-y[k],2);
//	func = erg+0.5*L*s;
//	for(int k=0;k<K;++k)
//		grad[k]=Dfx[TQlauf][k]+L*(y[k]-x[TQlauf][k]);
//}
//
//static void D_TQfunction1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *ptr)
//{
//	zeiger2->TQfunction1_grad(y,func,grad,NULL);
//}
//
double* AmericanOption::TQx(int lauf)
{
//	bool restringiert=false;
//	if(restringiert)
//	{
//		TQlauf=lauf;
//		real_1d_array u;
//		u.setlength(K);for(int k=0;k<K;++k)u[k]=0.;
//
//		real_1d_array bndl;
//		bndl.setlength(K);for(int k=0;k<K;++k)bndl[k]=-0.1;
//
//		real_1d_array bndu;
//		bndu.setlength(K);for(int k=0;k<K;++k)bndu[k]=0.1;
//		minbleicstate state;
//		minbleicreport rep;
//
//		double epsg = 0.001;
//		double epsf = 0;
//		double epsx = 0;
//		double epso = 0.001;
//		double epsi = 0.001;
//
//		minbleiccreate(u, state);
//		minbleicsetbc(state, bndl, bndu);
//		minbleicsetinnercond(state, epsg, epsf, epsx);
//		minbleicsetoutercond(state, epso, epsi);
//		alglib::minbleicoptimize(state, D_TQfunction1_grad);
//		minbleicresults(state, u, rep);
//
//		//printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
//		//printf("%s\n", u.tostring(2).c_str()); // EXPECTED: [-1,1]
//		double* e;
//		e=(double*)(malloc(sizeof(double)*K));
//		for(int k=0;k<K;++k)e[k]=u[k];
//		return e;
//	}else
	{//unrestringiert
		double* erg;
		erg=(double*)malloc(sizeof(double)*K);
		for(int k=0;k<K;++k)
			erg[k]=x[lauf][k]-Dfx[lauf][k]/Nesterov_L;
		return erg;
	}
}

void AmericanOption::Nesterov()
{
	//Nesterov Parameter
    if(testing || parallelTest)Nesterov_Iterations=400;else
	Nesterov_Iterations=4000;
	Nesterov_L=3000;

	zeiger2=this;
	BFGS_setting();
	Daten();

	if(testing  || extremTest)Nesterov_Iterations=BFGS_Iterations-2;

	fx=(double*)malloc(sizeof(double)*max(Nesterov_Iterations,705));
	x=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
	for(int i=0;i<max(Nesterov_Iterations,705);++i)
		x[i]=(double*)malloc(sizeof(double)*K);
	y=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
	for(int i=0;i<max(Nesterov_Iterations,705);++i)
		y[i]=(double*)malloc(sizeof(double)*K);
	z=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
	for(int i=0;i<max(Nesterov_Iterations,705);++i)
		z[i]=(double*)malloc(sizeof(double)*K);
	Dfx=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,705));
	for(int i=0;i<max(Nesterov_Iterations,705);++i)
		Dfx[i]=(double*)malloc(sizeof(double)*K);

	for(int k=0;k<K;++k)
		x[0][k]=0.000;

	if(verbose)printf("Alphas aus Datei laden\n");
	if(loadAlphas){
		double* alpha=alphasLaden(K);
		for (int k = 0; k < K; ++k) x[k][0]=alpha[k];
	}

	int M_store=M;
	if(verbose)printf("Nesterov\n");

	for(int lauf=0;(lauf<Nesterov_Iterations-1)||(speedup);++lauf)
	{
		printf("\n\nNesterov-Iteration: %d: \n",lauf);
		//Schritt 1

		if(speedup)
		{if(lauf<300)M=max(M_store/250,1000);
		if(lauf<200)M=M_store/100;
		if(lauf<100)M=M_store/5;
		if(lauf>=300){
			speedup=false;
			lauf=0;
		}
		}else M=M_store;
		TQlauf=argMinlauf=lauf;
		objfs_aufrufen(x[lauf],fx[lauf],Dfx[lauf]);


		//	printf("fx_%d = %.4lf (Minimum: %.4lf, diff: %.4lf)\n",lauf,fx[lauf],min, min-fx[lauf]);//printf("fx_%d= %f  (=%f)\n",lauf,fx[lauf],tan(fx[lauf]));
		if(verbose){printf("Df_%d =  ",lauf);ausgeben(Dfx[lauf],K);}

		//Schritt 2
		//		if(lauf<10)
		//		{
		y[lauf]=TQx(lauf);
		//		}else{double* y_; //TODO hier ist ein fehler drin
		//		y_=TQx(lauf);
		//		double m3;objfs(y_,m3,Dfx[lauf-2],MM);
		//		double m2;objfs(x[lauf],m2,Dfx[lauf-2],MM);
		//		double m1;objfs(y[lauf-1],m1,Dfx[lauf-2],MM);
		//		if(m1<m2)m1<m3?y[lauf]=y[lauf-1]:y[lauf]=y_;if(m1<m3)m1<m2?y[lauf]=y[lauf-1]:y[lauf]=x[lauf];
		//		//	printf("fy: %f\n",objfs(y[lauf]));
		//		//							if(fx[lauf]<min)
		//		//							{min=fx[lauf];argmin[k]=x[lauf]}}
		//		}
		//		printf("y_%d =   ",lauf);ausgeben(y[lauf],K);

		//Schritt 3
		z[lauf]=argminProblem(lauf);
		if(verbose){printf("z_%d =   ",lauf);ausgeben(z[lauf],K);}

		//Schritt 4
		for(int k=0;k<K;++k)
			x[lauf+1][k]=2./((double)lauf+3.)*z[lauf][k]+((double)lauf+1.)/((double)lauf+3.)*y[lauf][k];
		if(verbose){printf("x_%d+1 = ",lauf);ausgeben(x[lauf+1],K);}
                for(int k=0;k<K;++k)
                    alpha[k]=y[lauf][k];
	}
        alphasSchreiben(alpha,K);
}

//void AmericanOption::Lpruefen()
//{
//	double max=0;
//	for(int o=0;o<10;++o)
//	{
//		double x[K];
//		double y[K];
//		for(int k=0;k<K;++k)
//		{
//			x[k]=fabs(nextGaussian())*0.1;
//			y[k]=fabs(nextGaussian())*0.1;
//		}
//		//		if(          fabs(objfs(x[0])-objfs(y[0]))>100*d(x[0],y[0])				)
//		//printf("Fail-------------------------------------\n");
//		double d=0;
//		for(int k=0;k<K;++k)
//			d+=(x[k]-y[k])*(x[k]-y[k]);
//		d=sqrt(d);
//		objfsNesterov(x);
//		for(int k=0;k<K;++k)
//			x[k]=GRAD[k];
//		objfsNesterov(y);
//		for(int k=0;k<K;++k)
//			y[k]=GRAD[k];
//		double c=0;
//		for(int k=0;k<K;++k)
//			c+=pow(x[k]-y[k],2);
//		c=sqrt(c);
//		printf("max L: %f , max x:%f, max y:%f\n",max,Max(x,K),Max(y,K));
//		if(c/d>max)max=c/d;
//	}
//}

//void AmericanOption::objfsNesterov(double* alphas,double& func, double* gradient, int MM)
//{
//	//	real_1d_array x;
//	//	x.setlength(K);
//	//	real_1d_array grad;
//	//	grad.setlength(K);
//	//
//	//	for(int k=0;k<K;++k)
//	//	{
//	//		grad[k]=gradient[k];
//	//		x[k]=alphas[k];
//	//	}
//
//	objfs(alphas,func,gradient,MM);
//	//	for(int k=0;k<K;++k)
//	//	{
//	//		gradient[k]=grad[k];
//	//		alphas[k]=x[k];
//	//	}
//
//	//for(int k=0;k<K;++k)
//	//	GRAD[k]=grad[k];
//	//return func;
//
//	//	double gradient[K][M];
//	//	double sup_glatt[M];
//	//	double p=2;
//	//	time_t t;
//	//	t = time(NULL);
//	//	//	if(verbose)printf("Bisherige Rechenzeit: %ld Sekunden\n", t-ttt);
//	//	double mean_glatt = 0;
//	//	double mean_unglatt=0;
//	//	for(int k=0;k<K;++k)
//	//		GRAD[k]=0;
//	//
//	//	for (int m = 0; m < M; ++m) { //Summieren  ueber alle replications
//	//		double d[number_of_Exercise_Dates];
//	//		for(int ex=0;ex<number_of_Exercise_Dates;++ex)
//	//		{
//	//			d[ex]=payoff(Exercise_Dates[ex], X[m][Exercise_Dates[ex]]);
//	//			for (int k = 0; k < K; ++k)
//	//				d[ex] -= alphas[k] * StochIntegrals[k][m][Exercise_Dates[ex]];
//	//		}
//	//		double s=0;
//	//		for(int ex=0;ex<number_of_Exercise_Dates;++ex)
//	//			s+=exp((p * d[ex]));
//	//
//	//		for (int k = 0; k < K; ++k) {
//	//			double S =0;
//	//			for(int ex=0;ex<number_of_Exercise_Dates;++ex)
//	//				S+=StochIntegrals[k][m][Exercise_Dates[ex]] * exp(p * d[ex]);
//	//			gradient[k][m] = -S / s;
//	//		}
//	//
//	//		sup_glatt[m]= log(s) / p;
//	//		mean_glatt +=sup_glatt[m] /(double)(M);
//	//		mean_unglatt += Max(d,number_of_Exercise_Dates);
//	//	}
//	//	if(verbose)printf("Zwischenwert Ohne glaettung: %f\n",mean_unglatt/(double)(M));
//	//
//	//	//penelization...
//	//	//...fuer f
//	//	double lambda=1.;
//	//	double W=0;
//	//	for(int m=0;m<M;++m)
//	//		W+=pow(sup_glatt[m]-mean_glatt,2);
//	//	double erg=mean_glatt+lambda/sqrt(M)*sqrt(W);
//	//	if(verbose)printf("penelization term: %f\n", -mean_glatt+erg);
//	//
//	//	//...fuer den gradient
//	//	double mean_gradient[K];
//	//	for(int k=0;k<K;++k)
//	//	{
//	//		mean_gradient[k]=0;
//	//		for(int m=0;m<M;++m)
//	//			mean_gradient[k]+=gradient[k][m]/(double)M;
//	//	}
//	//	for(int k=0;k<K;++k)
//	//	{
//	//		GRAD[k]=0;
//	//		for(int m=0;m<M;++m)
//	//			GRAD[k]+=gradient[k][m]/(double)(M)+lambda/sqrt(M)*(sup_glatt[m]-mean_glatt)
//	//			*(gradient[k][m]-mean_gradient[k])/sqrt(W);
//	//	}
//	//
//	//	//arcus tangens
//	//	double q=10.;
//	//	bool atanAnwenden=false;
//	//	if(atanAnwenden)
//	//	{
//	//		double q=10.;
//	//		for (int k = 0; k < K; ++k)
//	//			GRAD[k] *=1./(1.+erg*erg/q/q)/q;
//	//		erg=atan(erg/q);
//	//	}
//	//
//	//	//for (int k = 0; k < K; ++k)
//	//	//	GRAD[k] *=1./(1.+erg*erg/q/q)/q;
//	//	//return atan(erg/q);
//	//	return erg;
//}

//	fx=(double*)malloc(sizeof(double)*max(Nesterov_Iterations,205));
//	x=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,205));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		x[i]=(double*)malloc(sizeof(double)*K);
//	y=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,205));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		y[i]=(double*)malloc(sizeof(double)*K);
//	z=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,205));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		z[i]=(double*)malloc(sizeof(double)*K);
//	Dfx=(double **)malloc(sizeof(double*)*max(Nesterov_Iterations,205));
//	for(int i=0;i<max(Nesterov_Iterations,705);++i)
//		Dfx[i]=(double*)malloc(sizeof(double)*K);

