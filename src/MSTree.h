
#ifndef MSTREE_H
#define MSTREE_H
#include "MSTreeNode.h"
#include "math.h"
#include <list>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cfloat>


class MSTree{
	public:
		MSTree();
		~MSTree();
		MSTree(const char *filename, int NumParam, string* names, double* scores, int arr_size, double* VP);
		
		void SetSigmaSq(double value);
		void SetLeafScores(string* names, double* scores, int arr_size, double* VP);
		double TreeScore(double* VP);
        double Integrat(double* VP);
		void ClearAllScores(double* VP, char* NumParam);
		double ComputeTreeScore();
    	void CalcAllValues(MSTreeNode* rt);
        double CalcRootValues();
		void PrintTree();
		void ResetAllSigmaSq(double* VP);
        int ParamNum();
        long double GammaPDF(double* VP);
        long double Gamma_integrat(double* VP);
		
		MSTreeNode* root;
        string* speName;
        int speNum;
        long double gamma;
		double sigmaSq;
		double tScore;
		double TBL; //sum of tree BL
        int NP;

	private:
		
		void ResetSigmaSq(MSTreeNode* rt, double* VP);
		void ClearScores(MSTreeNode* rt, double* VP);
		void PrintTree(MSTreeNode* rt);
	    	void SetLeafScore(string name, double score, MSTreeNode* rt, double* VP);
		MSTreeNode* BuildTreePreorder(FILE* f);
		void CalcUp(MSTreeNode* rt);
		double CalcDown(MSTreeNode* rt, MSTreeNode* parent, bool L);
};

#endif

