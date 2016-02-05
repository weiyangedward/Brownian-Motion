#include "NormD.h"

NormD::NormD()
{
	coeff = 1.0;
	mean = 0.0;
	var = 1.0;
}
NormD::~NormD()
{
}
NormD::NormD(NormD* newN)
{
	coeff = newN->coeff;
	mean = newN->mean;
	var = newN->var;
}

NormD::NormD(long double newCoeff, long double newMean, long double newVar)
{
	coeff = newCoeff;
	mean = newMean;
	var = newVar;
}

NormD::NormD(NormD* N1, NormD* N2) // the same as NormD::Set
{
	coeff = N1->coeff * N2->coeff * 1.0 / (sqrt((N1->var + N2->var) * 2.0 * PI)) * exp(-1.0 * pow(N1->mean - N2->mean,2.0) / (2.0*(N1->var + N2->var)));
	var = 1.0/((1.0/N1->var)+(1.0/N2->var));
	mean = var*(N1->mean/N1->var + N2->mean/N2->var);  
//    mean = (N1->mean/N1->var + N2->mean/N2->var);
}

void NormD::Set(long double newCoeff, long double newMean, long double newVar)
{
	coeff = newCoeff;
	mean = newMean;
	var = newVar;
}

void NormD::SetSig(long double newCoeff, long double newVar)
{
	coeff = newCoeff;
	var = newVar;
}
void NormD::Set(NormD* N1, NormD* N2) // calculate the scalers of two NormD
{
	//if (N1->var > 0 && N2->var > 0)
    //{
        //coeff = N1->coeff * N2->coeff * 1.0 / ((N1->var + N2->var) * sqrt(2 * PI)) * exp(-1.0 * pow(N1->mean - N2->mean,2.0) / (2.0*pow(N1->var + N2->var, 2.0)));
        coeff = N1->coeff * N2->coeff * 1.0 / (sqrt((N1->var + N2->var) * 2.0 * PI)) * exp(-1.0 * pow(N1->mean - N2->mean,2.0) / (2.0*(N1->var + N2->var)));
    //}
    //if (N1->var > 0 && N2->var < 0)
    //{
    //   coeff = N1->coeff * N2->coeff * 1.0 / ((N1->var - N2->var) * sqrt(2 * PI)) * exp(-1.0 * pow(N1->mean - N2->mean,2.0) / (2.0*pow(N1->var + N2->var, 2.0))); 
    //}
    //if (N1->var < 0 && N2->var > 0)
    //{
    //    coeff = N1->coeff * N2->coeff * 1.0 / ((-N1->var + N2->var) * sqrt(2 * PI)) * exp(-1.0 * pow(N1->mean - N2->mean,2.0) / (2.0*pow(N1->var + N2->var, 2.0)));
    //}
    //if (N1->var < 0 && N2->var < 0)
    //{
    //    coeff = N1->coeff * N2->coeff * 1.0 / ((-N1->var - N2->var) * sqrt(2 * PI)) * exp(-1.0 * pow(N1->mean - N2->mean,2.0) / (2.0*pow(N1->var + N2->var, 2.0)));
    //}
	var = 1.0/((1.0/N1->var)+(1.0/N2->var));
	mean = var*(N1->mean/N1->var + N2->mean/N2->var);
//    mean = (N1->mean/N1->var + N2->mean/N2->var);
}
void NormD::Set(NormD* newN)
{
	coeff = newN->coeff;
	mean = newN->mean;
	var = newN->var;
}
void NormD::Print()
{
	cout << coeff<<" "<<mean<< " "<<var<<"\n";
}
void NormD::PlotPrint() // print out normal distribution
{
	cout << "plot(x, " << coeff << "*dnorm(x, " << mean << ", "<<sqrt(var)<< "))\n";
}


