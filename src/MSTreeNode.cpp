#include "MSTreeNode.h"

MSTreeNode::MSTreeNode()
{
	BL=1.0;
	LeftChild = NULL;
	RightChild = NULL;
}

MSTreeNode::~MSTreeNode()
{
}

MSTreeNode::MSTreeNode(double bl, MSTreeNode* LC, MSTreeNode* RC, char* tname, double tmark)
{
	BL=bl;
	strcpy(name,tname);
	LeftChild=LC;
	RightChild=RC;

    Tmark = tmark;
    Tvector = tmark;
	
	B = new NormD();
	Bparent = new NormD();
	A = new NormD();
	Final = new NormD();
}


void MSTreeNode::Initialize(double sc, double sigmaSq)
{
	B->Set(1.0, sc, 0.0);
	Bparent->Set(1.0, sc, sigmaSq*BL);
	A->Set(1.0, 0.0, 1.0);
	Final->Set(1.0, 0.0, 1.0);
    Tvector = sigmaSq;
}

void MSTreeNode::InitializeSig(double newSig)
{
	B->SetSig(1.0, 0.0);
	Bparent->SetSig(1.0, newSig*BL);
	A->Set(1.0, 0.0, 1.0);
	Final->Set(1.0, 0.0, 1.0);
    Tvector = newSig;
}

void MSTreeNode::PrintNode()
{
	cout << "# " << name << " "<< BL << "\n";
	B->PlotPrint();
	Bparent->PlotPrint();
	A->PlotPrint();
	Final->PlotPrint();
	cout << "\n\n";
}

