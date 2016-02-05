#ifndef MSTREENODE_H
#define MSTREENODE_H
#include "NormD.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#define INF 1000000
using namespace std;

class MSTreeNode{
	public:
		MSTreeNode();
		~MSTreeNode();
		MSTreeNode(double bl=1.0, MSTreeNode* LC=NULL, MSTreeNode* RC=NULL, char* tname="node", double tmark=1.0);
		void InitializeSig(double newSig);
		void Initialize(double sc=0.0, double sigmaSq=1.0);
		void PrintNode();
		
		MSTreeNode* LeftChild;
		MSTreeNode* RightChild;
        double Tmark;
        double Tvector;
		NormD* B;
		NormD* Bparent;
		NormD* A;
		NormD* Final;
		char name[100];
		double BL; // Branch Length
};

#endif

