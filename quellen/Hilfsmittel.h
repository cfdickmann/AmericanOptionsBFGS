/*
 * Hilfsmittel.h
 *
 *  Created on: Feb 6, 2012
 *      Author: dickmann
 */

#include "math.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef HILFSMITTEL_H_
#define HILFSMITTEL_H_

int argMin(double* v, int l);
int argZweiter(double* v, int l);
int argDritter(double* v, int l);
int argMax(double* v, int l);

void MatrixAusgeben(double**  a, int D);
double** MatrixMultiplizieren(double**  a,double**  b,int D);
double betrag(double x);

void BubbleSort(double* werte, int* index, int l);
int* BubbleSort(double* werte, int l);

int* array_machen(int z);
//int hargMax(double* v, int l, int ll);
//int hMax(double* v, int l, int ll);
double summe(double* vector,int laenge);
double mean(double* vector,int laenge);

double Max(double* v, int l);
double* LGSloesen(double** A, double* b, int Mphi);
double Min(double* v, int l);
void ErgebnisAnhaengen(double d);
void ErgebnisAnhaengenML(double d);
double* alphasLaden(int K);
void alphasSchreiben(double* alpha,int K);
void werteSchreiben(double* w,int K, int N);
void ErgebnisAnhaengen(double d, char* filename);

//double qnorm(double p);

int * IntFeld(int m);
int ** IntFeld(int m,int n);
int *** IntFeld(int m,int n, int o);
double * DoubleFeld(int m);
double ** DoubleFeld(int m,int n);
double *** DoubleFeld(int m, int n, int o);
double **** DoubleFeld(int m, int n, int o, int p);
double ***** DoubleFeld(int m, int n, int o, int p, int q);

void deleteDoubleFeld(double * D, int m);
void deleteDoubleFeld(double ** D  ,int m,int n);
void deleteDoubleFeld(double *** D  ,int m, int n, int o);
void deleteDoubleFeld(double **** D , int m, int n, int o, int p);
void deleteDoubleFeld(double ***** D, int m, int n, int o, int p, int q);

void deleteIntFeld(int * D, int m);
void deleteIntFeld(int ** D  ,int m,int n);
void deleteIntFeld(int *** D  ,int m,int n,int o);

double** faurepts(int n0, int npts, int d, int b);
void ausgeben(double* x, int j);

void InPipeSchreiben(int* pipe, double wert );
double AusPipeLesen(int* pipe);

double max(double x, double y);

double CumulativeNormalDistribution(double x);


#endif /* HILFSMITTEL_H_ */
