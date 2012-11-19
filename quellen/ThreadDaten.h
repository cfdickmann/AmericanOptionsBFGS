/*
 * NummerUndZeiger.h
 *
 *  Created on: Oct 29, 2012
 *      Author: cfdickmann
 */

#ifndef THREADDATEN_H_
#define THREADDATEN_H_

class ThreadDaten {
public:
	double* ergebnis;
	int nummer;
	void setErgebnis(double D);
	double getErgebnis();

	ThreadDaten();
	ThreadDaten(int t, double* erg);
	virtual ~ThreadDaten();
};

#endif /* THREADDATEN_H_ */
