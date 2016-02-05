/***************************************************************************
 *   Questions please direct to weiyang4@illinois.edu                      *
 *   Author: Wei Yang                                                      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>

#include "MSTree.h"
#include "MSTMain.h"

#define MAX_N_SP 20 // max number of species
#define NCHRS 1 // number of chromosoms
string chrs[1]={"2L"};

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

using namespace std;

double my_f (const gsl_vector * v, void * params); // gsl my_function
double integration_f (double* x, size_t dim, void* params); //monte carlo integration function

int main(int argc, char **argv)
{
	if (argc < 4) 
    {
        printf("usage: %s <tree> <motif tab> <num of spe> <num of vectors>\n", argv[0]);
		exit(1);
	}

    int nSp = atoi(argv[3]);
    int NumParam = atoi(argv[4]);

    ifstream fs;
    fs.open(argv[2]);

    int line_total = 2 + nSp;
    int spe_start = 2;

    if (fs.is_open())
    {
        int line_num = 1;
        string line;
        string* spe_name;
        spe_name = new string [nSp];
        while (getline (fs, line))
        {
            double* sigma;
            sigma = new double [NumParam];
            string* arr;
            arr = new string [line_total];
            double* scores;
            scores = new double [nSp];
            stringstream ssin(line);
            int i = 0;
            if (line_num == 1)
            {
                while (ssin.good() && i < line_total)
                {
                    ssin >> arr[i];
                    cout << arr[i] << "\t";
                    i++;
                }
                for (int i=1; i<NumParam+1; i++)
                    cout << "estimated_sigma" << i << "\t";

                cout << "LL\t";
                cout << endl;

                for (int i = spe_start, j = 0; i < line_total; i++, j++)
                    spe_name[j] = arr[i];
            }
            else
            {
               double initial_sigma [] = {0.00001, 0.0001, 0.001, 0.01, 0.1}; // for bee10 motif scores
              //  double initial_sigma [] = {0.1, 0.5, 2.0, 4.0, 6.0}; // for hourglass 6 fly motif score and gene expression
                long double coeff = 0.0;
                double* final_sigma;
                final_sigma = new double[NumParam];
                double integration;
                for (int k = 0; k < 5; k++)
                {
                    while (ssin.good() && i < line_total)
                    {
                        ssin >> arr[i];
                        cout << arr[i] << "\t";
                        i++;
                    }
                    for (int j = 0; j < NumParam; j++)
                        sigma[j] = initial_sigma[k];

                    for (int i = spe_start, j=0; i < line_total; i++, j++)
                        scores[j] = atof(arr[i].c_str());

                    MSTree tree1(argv[1], NumParam, spe_name, scores, nSp, sigma);

                    /*======== start of gsl Multidimensional Minimizer =============*/
                    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
                    gsl_multimin_fminimizer *s ;
                    gsl_multimin_function ex4_fn;
                    gsl_vector *ss;

                    ss = gsl_vector_alloc (NumParam); /* Initial vertex size vector */
                    gsl_vector_set_all (ss, 1.0);

                    gsl_vector *xvec = gsl_vector_alloc (NumParam);
                    double* VecParam;
                    VecParam = new double[NumParam];
                    for (int i=0; i<NumParam; i++)
                    {
                        gsl_vector_set (xvec, i, sigma[i]);
                        VecParam[i] = gsl_vector_get(xvec, i);
                    }

                    ex4_fn.f = &my_f;
                    ex4_fn.n = NumParam;
                    ex4_fn.params = &tree1;
                    s = gsl_multimin_fminimizer_alloc (T, NumParam);
                    gsl_multimin_fminimizer_set (s, &ex4_fn, xvec, ss);

                    int iter =0;
                    int status;
                    double size;

                    do
                    {
                        iter++;
                        status = gsl_multimin_fminimizer_iterate(s); // performs one iteration to update the state of the minimizer, and the vector corresponding to the highest function value.

                        if (status) // If the iteration encounters an unexpected problem then an error code will be returned
                        {
                            cout << "encounters an unexpected problem: ";
                            cout << GSL_ENOPROG << endl;
                            break;
                        }

                        size = gsl_multimin_fminimizer_size (s);
                        status = gsl_multimin_test_size (size, 1e-5); // The minimizer-specific characteristic size is calculated as the average distance from the geometrical center of the simplex to all its vertices. This size can be used as a stopping criteria, as the simplex contracts itself near the minimum. RETURN 0 if converged (size <= 0.01).

                        if (status == GSL_SUCCESS) //GSL_SUCCESS = 0
                            break;
                    }
                    while (status == GSL_CONTINUE && iter < 10000); // GSL_CONTINUE = -2

                    if (coeff < (- s->fval))
                    {
                        double result, error;
                        coeff = - s->fval;
                        for (int i=0; i<NumParam; i++)
                            final_sigma[i] = gsl_vector_get (s->x, i);
                    }
                    gsl_vector_free(xvec);
                    gsl_vector_free(ss);
                    gsl_multimin_fminimizer_free (s);
                    delete [] VecParam;
                }
                
                for (int i=0; i<NumParam; i++)
                {
                    printf ("%10.3e\t", final_sigma[i]);
                }
                printf ("%7.5f\t", log(coeff));
                cout << endl;
                delete [] final_sigma;
            }
            line_num++;
            delete [] arr;
            delete [] scores;
            delete [] sigma;
        }
        delete [] spe_name;
        fs.close();
    }
	return EXIT_SUCCESS;
    
}

double my_f (const gsl_vector * v, void * params)
{
    double* VecParam_sub;
    MSTree * t = (MSTree *)params;
    int NP;
    NP = t->ParamNum();
    VecParam_sub = new double[NP];
    for (int i=0; i<NP ; i++)
        VecParam_sub[i] = gsl_vector_get(v, i);
    return (-t->TreeScore(VecParam_sub));
}

double integration_f (double* x, size_t dim, void* params)
{
    double* vector = x;
    MSTree * t = (MSTree *)params;
    return (t->Integrat(vector));
}
