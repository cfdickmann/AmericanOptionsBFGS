/*
 * EuroBewerter.h
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#ifndef EUROBEWERTER_H_
#define EUROBEWERTER_H_

namespace std {

class EuroBewerter {
public:
	EuroBewerter();
	virtual ~EuroBewerter();

	double call(double t, double T, double X0, double Strike, double r, double delta, double sigma);
	double call_diff(double t, double T, double X0, double Strike, double r, double delta, double sigma);
	double put(double t, double T, double X0, double Strike, double r, double delta, double sigma);
	double max_call(double t, double T, double* X0, int D, double Strike, double r, double delta, double sigma);
};

} /* namespace std */
#endif /* EUROBEWERTER_H_ */
