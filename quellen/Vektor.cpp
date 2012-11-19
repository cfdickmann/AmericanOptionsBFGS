/*
 * Vektor.cpp
 *
 *  Created on: Oct 29, 2012
 *      Author: cfdickmann
 */

#include "Vektor.h"
#include <stdlib.h>

Vektor::Vektor() {
	eintraege=NULL;
}
//
//Vektor::Vektor(int laenge):Vektor() {
//	eintraege=NULL;
//	setlaenge(laenge);
//}

void Vektor::setlaenge(int l){
	if( eintraege!=NULL)delete[] eintraege;
	eintraege=new double[l];
	for(int i=0;i<laenge;++i)
		eintraege[l]=0;
	laenge=l;
}

Vektor::~Vektor() {
	delete[] eintraege;
}

double Vektor::get(int i){
	return eintraege[i];
}

double Vektor::sum(){
double s=0;
for(int i=0;i<laenge;++i)
	s+=eintraege[i];
return s;
}

double Vektor::mean(){
	return sum()/(double)(laenge);
}

void Vektor::set(int i, double wert){
	eintraege[i]=wert;
}

void Vektor::add(int i, double wert){
	eintraege[i]+=wert;
}

int Vektor::getLaenge(){
	return laenge;
}
