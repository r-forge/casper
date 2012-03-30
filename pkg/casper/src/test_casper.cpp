#include "seppelsmart.h"
#include <cstdlib>
#include <stdio.h>
//#include <tchar.h>
#include <time.h>

using namespace std;

double fragsta_cumu(double x)
{
	return x;
}

Casper* example()
{	
	double* fraglens = new double[1];
        int* lenvals= new int[1];
        fraglens[0]= 1;
        lenvals[0]= 200;
	DiscreteDF* fraglen_dist = new DiscreteDF(fraglens, lenvals, 1);

	DataFrame* f = new DataFrame(fraglen_dist, fragsta_cumu);

	Exon* e1 = new Exon(1, 300);
	Exon* e2 = new Exon(2, 100);
	Exon* e3 = new Exon(3, 500);
	f->addExon(e1);
	f->addExon(e2);
	f->addExon(e3);

	Gene* g1 = new Gene(1);
	g1->addExon(e1);
	g1->addExon(e2);
	g1->addExon(e3);
	f->addGene(g1);

	Fragment* f0 = new Fragment(1, 1, 1824);
	f0->left[0] = 1;
	f0->right[0] = 1;
	f->addData(f0);
	Fragment* f1 = new Fragment(1, 1, 4068);
	f1->left[0] = 3;
	f1->right[0] = 3;
	f->addData(f1);
	Fragment* f2 = new Fragment(2, 1, 649);
	f2->left[0] = 2;
	f2->left[1] = 3;
	f2->right[0] = 3;
	f->addData(f2);
	Fragment* f3 = new Fragment(1, 2, 339);
	f3->left[0] = 1;
	f3->right[0] = 1;
	f3->right[1] = 3;
	f->addData(f3);
	Fragment* f4 = new Fragment(1, 2, 1013);
	f4->left[0] = 1;
	f4->right[0] = 1;
	f4->right[1] = 2;
	f->addData(f4);
	Fragment* f5 = new Fragment(1, 1, 368);
	f5->left[0] = 1;
	f5->right[0] = 2;
	f->addData(f5);
	Fragment* f6 = new Fragment(1, 1, 268);
	f6->left[0] = 1;
	f6->right[0] = 3;
	f->addData(f6);
	Fragment* f7 = new Fragment(1, 2, 202);
	f7->left[0] = 1;
	f7->right[0] = 2;
	f7->right[1] = 3;
	f->addData(f7);
	Fragment* f8 = new Fragment(2, 2, 430);
	f8->left[0] = 1;
	f8->left[1] = 2;
	f8->right[0] = 2;
	f8->right[1] = 3;
	f->addData(f8);
	Fragment* f9 = new Fragment(2, 1, 368);
	f9->left[0] = 1;
	f9->left[1] = 3;
	f9->right[0] = 3;
	f->addData(f9);
	Fragment* f10 = new Fragment(1, 1, 240);
	f10->left[0] = 2;
	f10->right[0] = 3;
	f->addData(f10);
	Fragment* f11 = new Fragment(2, 1, 231);
	f11->left[0] = 1;
	f11->left[1] = 2;
	f11->right[0] = 3;
	f->addData(f11);

	vector<Exon*>* v1v = new vector<Exon*>();
	v1v->push_back(e1);
	v1v->push_back(e2);
	v1v->push_back(e3);
	vector<Exon*>* v2v = new vector<Exon*>();
	v2v->push_back(e1);
	v2v->push_back(e3);
	vector<Exon*>* v3v = new vector<Exon*>();
	v3v->push_back(e1);
	v3v->push_back(e2);

	Variant* v1 = new Variant(g1, v1v);
	v1->id = 1;
	Variant* v2 = new Variant(g1, v2v);
	v2->id = 2;
	Variant* v3 = new Variant(g1, v3v);
	v3->id = 3;

	vector<Variant*>* varis = new vector<Variant*>();
	varis->push_back(v1);
	varis->push_back(v2);
	//varis->push_back(v3);

	Model* model = new Model(varis);
	Casper* casp = new Casper(model, f);

	return casp;
}

int main() {
//int _tmain(int argc, _TCHAR* argv[])
	srand((unsigned)time( NULL ));

	Casper* c = example();
	double* pi = c->calculateMode();
	//double inte = c->calculateIntegral();

	int n = c->model->count();

	for (int i = 0; i < n; i++)
	{
		printf("%f\n", pi[i]);
	}
	printf("\n");

	//vector<Model*>* m1 = c->frame->allModels(c->frame->genes[1]);
	//vector<Model*>* m2 = c->frame->allModels(c->frame->genes[1]);
	//Model* n1 = m1->at(0);
	//Model* n2 = m2->at(0);

	SeppelPrior* sep2 = new SeppelPrior(c->frame, c->frame->genes[1]);
	map<Model*, double, ModelCmp> res2 = sep2->calculate();

	printf("PRIOR\n");

	SeppelSmart* sep3 = new SeppelSmart(c->frame, c->frame->genes[1]);
	sep3->calculate(c->model);
	map<Model*, double, ModelCmp> res3 = sep3->resProbs;

	printf("SMART\n");

	SeppelExact* sep1 = new SeppelExact(c->frame, c->frame->genes[1]);
	map<Model*, double, ModelCmp> res1;
	sep1->calculate();
	res1 = sep1->resProbs;

	printf("EXACT\n");
	
	map<Model*, double, ModelCmp>::const_iterator mi;
	for (mi = res1.begin(); mi != res1.end(); mi++)
	{
		Model* m = mi->first;
		printf("%f\t%f\t%f\n", res1[m], res2[m], res3[m]);
	}

	getchar();

	return 0;
}
