/* 
 * File:   LCG_Generator.cpp
 * Author: dickmann
 * 
 * Created on August 6, 2012, 1:51 PM
 */

#include "RNG.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

RNG::RNG() {
    int r = rand();
    setSeed(getpid() + r + time(NULL));
}

RNG::RNG(const RNG& orig) {
}

RNG::~RNG() {
}

void RNG::setSeed(int seed) {
    mt.init_genrand(seed);
    // Ein Mersennetwister braucht bei schlechten seeds bis zu
    // 600 Durchlaeufe um auf Temperatur zu kommen.
    for (int lauf = 0; lauf < 600; ++lauf)
        nextGaussian();
}

double RNG::BoxMuller(double U1, double U2) {
    if (U1 == 0)U1 = 0.000001;
    double R = -2 * log(U1);
    double V = 2. * 3.1415926535 * U2;
    return sqrt(R) * cos(V);
}

double RNG::nextUnif() {
    return mt.random();
}

double RNG::nextGaussian() {
    return BoxMuller(nextUnif(), nextUnif());
}

//x=seed;

//   Xx ^= Xx << 16;
//   Xx ^= Xx >> 5;
//   Xx ^= Xx << 1;
//   unsigned long t;
//   t = Xx;
//   Xx = Xy;
//   Xy = Xz;
//   Xz = t ^ Xx ^ Xy;
//
//  return (double)(Xz%100000)/100000.;



    //      MT[0] = seed;
    //     for i from 1 to 623 { // loop over each other element
    //         MT[i] := last 32 bits of(1812433253 * (MT[i-1] xor (right shift by 30 bits(MT[i-1]))) + i) // 0x6c078965
    //     }


//Bauanleitung auf Wikipedia
//void LCG_Generator::generate_numbers() {
    //     for(int i=0;i<624++i) {
    //         int y := 32nd bit of(MT[i]) + last 31 bits of(MT[(i+1) mod 624])
    //         MT[i] := MT[(i + 397) % 624] ^ (y>>1))
    //         if (y % 2 != 0) { // y is odd
    //             MT[i] := MT[i] ^ (2567483615) // 0x9908b0df
    //         }
//}


    //
    //   x=(a*x+c)>>16;
    //    if(x>=0)
    //   return (double)(x%100000)/100000.;
    //    else
    //        return (double)(-x%1000000)/1000000.;


    //      if index == 0 {
    //         generate_numbers()
    //     }
    //
    //     int y := MT[index]
    //     y := y ^ (y>>11)
    //     y := y ^ ((y<<7) & 2636928640) // 0x9d2c5680
    //     y := y ^ ((y<<15)& 4022730752) // 0xefc60000
    //     y := y ^ (y>>18)
    //
    //     index := (index + 1) % 624
    //     return y




    //   x=1;
    //   c=1013904223;
    //   a=1664525;
    //    MT=new int[624];
    //    index=0;
