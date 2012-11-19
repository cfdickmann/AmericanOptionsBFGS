/* 
 * File:   LCG_Generator.h
 * Author: dickmann
 *
 * Created on August 6, 2012, 1:51 PM
 */

#ifndef LCG_GENERATOR_H
#define	LCG_GENERATOR_H

#include "mt.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "math.h"
#include "stdlib.h"
#include "Hilfsmittel.h"

class RNG {
public:
    RNG();
    RNG(const RNG& orig);
    virtual ~RNG();
    MersenneTwister mt;
    void setSeed(int seed);
    double nextUnif();
      double BoxMuller(double U1, double U2);
    double nextGaussian();

    void generate_numbers();
private:
   long int x,a,c,m;

//    int* MT;
// int index;


};

#endif	/* LCG_GENERATOR_H */

