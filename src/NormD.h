#ifndef NORMD_H
#define NORMD_H
#include "math.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;
#define PI 3.1415926535897932384626433832795

class NormD{
	public:
		NormD();
		~NormD();
		NormD(long double newCoeff, long double newMean, long double newVar);
		NormD(NormD* N1, NormD* N2);
		NormD(NormD* N1);
		
		long double mean;
		long double var;
		long double coeff;
        double sig;

		void SetSig(long double newCoeff, long double newVar);
		void Set(long double newCoeff, double long newMean, long double newVar);
		void Set(NormD* N1, NormD* N2);
		void Set(NormD* N1);
		void Print();
		void PlotPrint();		
};

#endif
