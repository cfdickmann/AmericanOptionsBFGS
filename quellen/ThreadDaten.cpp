/*
 * NummerUndZeiger.cpp
 *
 *  Created on: Oct 29, 2012
 *      Author: cfdickmann
 */

#include "ThreadDaten.h"

ThreadDaten::ThreadDaten() {
	// TODO Auto-generated constructor stub
}

ThreadDaten::ThreadDaten(int t, double * erg) {
	nummer=t;
	ergebnis=erg;
}

ThreadDaten::~ThreadDaten() {
	// TODO Auto-generated destructor stub
}

void ThreadDaten::setErgebnis(double D){
	ergebnis[0]=D;
}

double ThreadDaten::getErgebnis(){
	return ergebnis[0];
}
