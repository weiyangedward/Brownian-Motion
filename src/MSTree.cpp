#include "MSTree.h"

MSTree::MSTree(){}

MSTree::~MSTree()
{
	// I should remove 5(free) the tree here recursively
}

MSTree::MSTree(const char *filename, int NumParam, string* names, double* scores, int arr_size, double* VP)
{
	FILE *f=fopen(filename,"r");
	TBL=0.0;
	sigmaSq = 1.0;
	tScore = 0.0;
    NP = NumParam;
	
	root=new MSTreeNode(0, NULL, NULL, "root", 0.1); // initialize the value of the root (1, 0, 1), BL=1

	if(!feof(f)){
		root->LeftChild=BuildTreePreorder(f); // leftchild of root
		root->RightChild=BuildTreePreorder(f); // rightchild of root
	}
	fclose(f);
    speName = names;
    speNum = arr_size;

    SetLeafScores(names, scores, arr_size, VP);
}

MSTreeNode* MSTree::BuildTreePreorder(FILE* f) 
{
	if(feof(f)) 
        return NULL;
	
	char blstr[100], *spename, *mark_char;
    char name[100];
	double bl, mark;

	fgets (name , 100 , f); // node name
    spename = strtok (name," ");
    mark_char = strtok (NULL," \n");
    mark = atof(mark_char);

	fgets (blstr , 100 , f); // branch length
	bl=atof(blstr); // bl = branch length

	if(bl==0) 
        return NULL; // return NULL to the children of leaf

    TBL+=bl;
    MSTreeNode* rt=new MSTreeNode(bl, NULL, NULL, spename, mark);
	rt->LeftChild=BuildTreePreorder(f); // leftchild of node
	rt->RightChild=BuildTreePreorder(f); // rightchild of node
	return rt;
}

void MSTree::SetSigmaSq(double value)
{
	sigmaSq = value;
}

void MSTree::SetLeafScores(string* names, double* scores, int arr_size, double* VP)
{
    for (int i=0; i<arr_size; i++)
		SetLeafScore(names[i], scores[i],root, VP); // set leaf score for each leaf
}

void MSTree::SetLeafScore(string name, double score, MSTreeNode* rt, double* VP)
{
	if (rt==NULL) 
        return;
    
    if (rt->name == name)
    {
		for (int i=0; i<NP; i++)
        {
            if (i == rt->Tmark)
            {
                rt->Initialize(score, VP[i]);
                rt->Tvector = VP[i];
            }
        }
	}
	SetLeafScore(name, score, rt->LeftChild, VP);
	SetLeafScore(name, score, rt->RightChild, VP);
}

void MSTree::ClearAllScores(double* VP, char* NumParam)
{
	NP = atof(NumParam);
    //cout << NP << "..." << endl;
    ClearScores(root, VP);
	tScore = 0.0;
}


void MSTree::ClearScores(MSTreeNode* rt, double* VP) // initialize the score
{
	if (rt==NULL) return;
    for (int i=0; i<NP; i++)
    {
	    //cout << "i=" << i <<" Tmark=" << rt->Tmark << " VP[Tmark]=" << VP[i]  << "..." << endl;
        if (i == rt->Tmark)
        {
            rt->Initialize(0.0, VP[i]);
            rt->Tvector = VP[i];
            //cout << "Tmark=" << *(rt->Tmark) << " VP[Tmark]=" << VP[i]  << "..." << endl;
        }
    }       // initialize the var = sigma*bl
	ClearScores(rt->LeftChild, VP);
	ClearScores(rt->RightChild, VP);
}

void MSTree::ResetAllSigmaSq(double* VP)
{
	MSTreeNode* reset_root = root;
    reset_root->Initialize(0.0, 1.0);
    ResetSigmaSq(root, VP);
	tScore = 0.0;
}


void MSTree::ResetSigmaSq(MSTreeNode* rt, double* VP)
{
	int k = 0;
    if (rt==NULL) 
        return;

    for (int i=0; i<NP; i++)
    {
	    if (i == rt->Tmark)
        {
            for (int j=0; j<speNum; j++)
            {
                if (speName[j] == rt->name)
                {
                    rt->InitializeSig(VP[i]);
                    k = 1;
                }
            }

            if (k == 0)
                rt->Initialize(-0.0000000001, VP[i]);

            rt->Tvector = VP[i];
        }
    }
	ResetSigmaSq(rt->LeftChild, VP);
	ResetSigmaSq(rt->RightChild, VP);
}

void MSTree::PrintTree()
{
	PrintTree(root); // call MSTree::PrintTree(MSTreeNode* rt)
	cout<< "sigmaSq " << sigmaSq << " Score "<<tScore<<" Estimate " <<tScore/TBL<<"\n\n\n";
}

void MSTree::PrintTree(MSTreeNode* rt)
{
	if (rt==NULL) return;
	rt->PrintNode(); // call PrintNode, print out the NormD
	PrintTree(rt->LeftChild);
	PrintTree(rt->RightChild);
}

double MSTree::CalcRootValues()
{
    double constant = 1.0/20.0;    
    root->Final->Set(constant*root->B->coeff, root->B->mean, root->B->var);
    return root->Final->coeff;
}

void MSTree::CalcAllValues(MSTreeNode* rt)
{
	double constant = 1.0/20.0;
	
    if(rt->RightChild == NULL || (rt->RightChild->LeftChild == NULL && rt->RightChild->RightChild == NULL)) // only calculate the nodes' final
        return;
    else{
        rt->RightChild->A->Set(constant*rt->LeftChild->Bparent->coeff, rt->LeftChild->Bparent->mean, rt->LeftChild->Bparent->var + rt->RightChild->Tvector * rt->RightChild->BL);
        NormD* temp = new NormD(rt->RightChild->A, rt->RightChild->B);
        rt->RightChild->Final->Set(temp);
        tScore= tScore+(rt->Final->mean + rt->RightChild->Final->mean)/2.0*rt->RightChild->BL;
    }

	if(rt->LeftChild == NULL || (rt->LeftChild->LeftChild == NULL && rt->LeftChild->RightChild == NULL)) // only calculate the nodes' final
		return;
	else{
		rt->LeftChild->A->Set(constant*rt->RightChild->Bparent->coeff, rt->RightChild->Bparent->mean, rt->RightChild->Bparent->var + rt->LeftChild->Tvector * rt->LeftChild->BL);
		NormD* temp = new NormD(rt->LeftChild->A, rt->LeftChild->B); 
		rt->LeftChild->Final->Set(temp);
		tScore= tScore+(rt->Final->mean + rt->LeftChild->Final->mean)/2.0*rt->LeftChild->BL;
	}
	CalcAllValues(rt->LeftChild);
    CalcAllValues(rt->RightChild);

	tScore = tScore + CalcDown(rt->LeftChild->LeftChild,rt->LeftChild, 1);
	tScore = tScore + CalcDown(rt->LeftChild->RightChild,rt->LeftChild, 0);
	tScore = tScore + CalcDown(rt->RightChild->LeftChild,rt->RightChild, 1);
	tScore = tScore + CalcDown(rt->RightChild->RightChild,rt->RightChild, 0);
}

void MSTree::CalcUp(MSTreeNode* rt) // return nothing !1
{
	if(rt->LeftChild==NULL && rt->RightChild==NULL) return;

    if(rt->LeftChild->B->mean != -DBL_MIN && rt->RightChild->B->mean != -DBL_MIN && rt->B->var == 0.0)
    {
        rt->B->Set(rt->LeftChild->Bparent, rt->RightChild->Bparent);
        rt->Bparent->Set(rt->B->coeff, rt->B->mean, rt->Tvector * rt->BL + rt->B->var);
    }
    
    CalcUp(rt->LeftChild);
    CalcUp(rt->RightChild);
}

double MSTree::CalcDown(MSTreeNode* rt, MSTreeNode* parent, bool L)
{
	if(rt==NULL) return 0;

	if(L == 1)
	{
		NormD* temp = new NormD(parent->A, parent->RightChild->Bparent); 
		rt->A->Set(temp->coeff,temp->mean,temp->var+sigmaSq*rt->BL);
	}
	else
	{
		NormD* temp = new NormD(parent->A, parent->LeftChild->Bparent); 
		rt->A->Set(temp->coeff,temp->mean,temp->var+sigmaSq*rt->BL);
	}

	if(rt->LeftChild == NULL && rt->RightChild==NULL)
		rt->Final->Set(rt->B);
	else{
		NormD* temp2 = new NormD(rt->A, rt->B); 
		rt->Final->Set(temp2);
	}

	double LScore = CalcDown(rt->LeftChild, rt, 1);
	double RScore = CalcDown(rt->RightChild, rt, 0);

	return (parent->Final->mean + rt->Final->mean)/2.0*rt->BL + LScore + RScore;
}

long double MSTree::GammaPDF(double* VP)
{
    long double totalGamma = 1.0;
    long double tmpGamma = 0.0;
    for (int i=0; i<NP; i++)
    {
        tmpGamma = (40000.0 * exp(-200.0 * VP[i]) * VP[i]); // alpha=2, beta=0.005
        if (VP[i] > 0.0)
        {
            totalGamma *= tmpGamma;
        }
        else
        {
            totalGamma *= 0.0;
        }

    }
    return totalGamma;
}


long double MSTree::Gamma_integrat(double* VP)
{
    long double totalGamma = 1.0;
    long double tmpGamma = 0.0;
    for (int i=0; i<NP; i++)
    {
        tmpGamma = (40000.0 * exp(-200.0 * VP[i]) * VP[i]); // alpha=2, beta=0.005
        totalGamma *= tmpGamma;
    }
    return totalGamma;
}

double MSTree::TreeScore(double* VP) // call CalcAllValues(); input sigma!!
{
	ResetAllSigmaSq(VP);
    tScore = 0.0;
    gamma = 0.0;
    gamma = GammaPDF(VP);
    while (root->B->mean == 0)
        CalcUp(root); // call CalcUp

    tScore = CalcRootValues();
    int key=0;
    int tScore_negative=0;
    for (int i =0; i< NP; i++)
    {
        if (VP[i] < 0)
            key=1;
    }

    if (key==0)
        return (tScore);
    else
        return (tScore_negative);
}

double MSTree::Integrat(double* VP)
{
    ResetAllSigmaSq(VP);
    tScore = 0.0;
    gamma = 0.0;
    gamma = Gamma_integrat(VP);
    while (root->B->mean == 0.0)
        CalcUp(root); // call CalcUp

    tScore = CalcRootValues();
    return (tScore * gamma);
}


int MSTree::ParamNum()
{
    return NP;
}

