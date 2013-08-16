/*
 * Polynom.h
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#ifndef POLYNOM_H_
#define POLYNOM_H_

namespace std {
class Polynom;

class Polynom {
public:
	Polynom();
	virtual ~Polynom();
	int length;
	double* koeff;
	void ausgeben();
	void add(Polynom* b);
	void malnehmen(double c);
	double auswerten(double u);
	void multiply_with(Polynom* p);
	void lerweitern(int l);
	void verschieben(double v);
	void set(double *a, int l);
};

} /* namespace std */
#endif /* POLYNOM_H_ */
