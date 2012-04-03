#include "seppel.h"
#include <cstdlib>
#include <stdio.h>
//#include <tchar.h>
#include <time.h>

using namespace std;

double fragsta_cumu(double x)
{
	return x;
}
void debugdf(DataFrame* df)
{		
	// Exons
	printf("Exons:\n");
	vector<Exon*>::const_iterator ei;
	for (ei = df->exons.begin(); ei != df->exons.end(); ei++)
	{
		Exon* e = *ei;
		printf("%i\t%i\n", e->id, e->length);
	}
	printf("\n");

	// Fragments
	printf("Fragments:\n");
	list<Fragment*>::const_iterator fi;
	for (fi = df->data.begin(); fi != df->data.end(); fi++)
	{
		Fragment* f = *fi;
		printf("%i\t%i\t%i\n", f->leftc, f->rightc, f->count);
		for (int l = 0; l < f->leftc; l++)
		{
			printf("%i\n", f->left[l]);
		}
		for (int r = 0; r < f->rightc; r++)
		{
			printf("%i\n", f->right[r]);
		}
	}
	printf("\n");
}
void debugmodel(Model* model)
{
	// Transcripts
	printf("Model:\n");
	vector<Variant*>::const_iterator vi;
	for (vi = model->items.begin(); vi != model->items.end(); vi++)
	{
		Variant* v = *vi;
		printf("%i\t%i\n", v->id, v->exonCount);
		for (int e = 0; e < v->exonCount; e++)
		{
			printf("%i\n", v->exons[e]->id);
		}
	}
	printf("\n");
}

Variant* path2Variant(DataFrame *df, Fragment* f) 
{
	int eid; Exon *ex;
	vector<Exon*>::iterator itexon;
	map<int, Exon*> id2exon;
	for (itexon= df->exons.begin(); itexon != df->exons.end(); itexon++) {
		ex= (*itexon);
		id2exon[ex->id] = ex;
	}
	vector<Exon*>* el = new vector<Exon*>();
	for (itexon= df->exons.begin(); (*itexon)->id != f->left[0]; itexon++) {
		ex= (*itexon);
		el->push_back(ex);
	}
	for (int i=0; i< f->leftc; i++) {
		eid = f->left[i];
		ex = id2exon[eid];
		el->push_back(ex);
	}
	if (eid != f->right[0]) {
		eid = f->right[0];
		ex = id2exon[eid];
		el->push_back(ex);
	}
	for (int i=1; i< f->rightc; i++) {
		eid = f->right[i];
		ex = df->exons[eid];
		el->push_back(ex);
	}
	while ((*itexon)->id != eid) { itexon++; }
	itexon++;
	while (itexon != df->exons.end()) {
		ex= (*itexon);
		el->push_back(ex);
		itexon++;
	}
	Variant* v = new Variant(el);
	return v;
}
int fixUnexplFrags(set<Variant*, VariantCmp>* initvars, DataFrame* df)
{
	// copy all fragments
	set<Fragment*>* queue = new set<Fragment*>(df->data.begin(), df->data.end());

	// remove the fragments from the queue we can explain with our variants
	set<Variant*, VariantCmp>::iterator vi;
	for (vi = initvars->begin(); vi != initvars->end(); vi++) 
	{
		// remove the fragments that this variant can explain from our queue
		map<Fragment*, double> probs = df->probabilities(*vi);
		map<Fragment*, double>::iterator si;
		for (si = probs.begin(); si != probs.end(); si++) 
		{
			set<Fragment*>::iterator ri = queue->find(si->first);
			if (ri != queue->end())
			{
				queue->erase(ri);
			}
		}
	}

	// remove variants from model that cant be explained
	int discarded = 0;

	// while we still have unexplained fragments
	while (queue->size() > 0)
	{
		// pop the first fragment
		Fragment* frag = *queue->begin();
		queue->erase(queue->begin());

		Variant* nv = path2Variant(df, frag);

		// check if the new variant can explain the fragment
		map<Fragment*, double> probs = df->probabilities(nv);
		if (probs.count(frag) > 0)
		{
			initvars->insert(nv);

			// delete all fragments that this variant can explain
			map<Fragment*, double>::iterator si;
			for (si = probs.begin(); si != probs.end(); si++) 
			{
				set<Fragment*>::iterator ri = queue->find(si->first);
				if (ri != queue->end())
				{
					queue->erase(ri);
				}
			}
		}
		else
		{
			// this fragment cant be explained
			discarded++;
			df->data.remove(frag);
		}
	}

	return discarded;
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

	Variant* v1 = new Variant(v1v);
	v1->id = 1;
	Variant* v2 = new Variant(v2v);
	v2->id = 2;
	Variant* v3 = new Variant(v3v);
	v3->id = 3;

	vector<Variant*>* varis = new vector<Variant*>();
	varis->push_back(v1);
	varis->push_back(v2);
	//varis->push_back(v3);

	Model* model = new Model(varis);
	Casper* casp = new Casper(model, f);

	return casp;
}
Casper* david1()
{	
	double* fraglens = new double[1];
        int* lenvals= new int[1];
        fraglens[0]= 1;
        lenvals[0]= 200;
	DiscreteDF* fraglen_dist = new DiscreteDF(fraglens, lenvals, 1);

	DataFrame* f = new DataFrame(fraglen_dist, fragsta_cumu);

	Exon* e1 = new Exon(2181, 73);
	Exon* e2 = new Exon(2182, 144);
	Exon* e3 = new Exon(2183, 4868);
	Exon* e4 = new Exon(2184, 462);
	f->addExon(e1);
	f->addExon(e2);
	f->addExon(e3);
	f->addExon(e4);

	Fragment* f0 = new Fragment(1, 1, 500);
	f0->left[0] = 2181;
	f0->right[0] = 2181;
	f->addData(f0);
	Fragment* f1 = new Fragment(1, 1, 100);
	f1->left[0] = 2181;
	f1->right[0] = 2182;
	f->addData(f1);
	Fragment* f2 = new Fragment(2, 1, 20);
	f2->left[0] = 2181;
	f2->left[1] = 2182;
	f2->right[0] = 2183;
	f->addData(f2);
	Fragment* f3 = new Fragment(1, 1, 200);
	f3->left[0] = 2182;
	f3->right[0] = 2182;
	f->addData(f3);
	Fragment* f4 = new Fragment(1, 1, 20);
	f4->left[0] = 2182;
	f4->right[0] = 2183;
	f->addData(f4);
	Fragment* f5 = new Fragment(1, 1, 1000);
	f5->left[0] = 2183;
	f5->right[0] = 2183;
	f->addData(f5);
	Fragment* f6 = new Fragment(1, 1, 50);
	f6->left[0] = 2181;
	f6->right[0] = 2183;
	f->addData(f6);
	Fragment* f7 = new Fragment(1, 1, 100);
	f7->left[0] = 2181;
	f7->right[0] = 2184;
	f->addData(f7);

	vector<Exon*>* v1v = new vector<Exon*>();
	v1v->push_back(e1);
	v1v->push_back(e2);
	v1v->push_back(e3);
	v1v->push_back(e4);

	Variant* v1 = new Variant(v1v);
	v1->id = 1;

	set<Variant*, VariantCmp>* varis = new set<Variant*, VariantCmp>();
	varis->insert(v1);
	//varis->push_back(v3);

	debugdf(f);
	debugmodel(new Model(varis));

	fixUnexplFrags(varis, f);
	debugmodel(new Model(varis));

	Model* model = new Model(varis);
	Casper* casp = new Casper(model, f);

	return casp;
}
Casper* david2()
{	
	double* fraglens = new double[1];
        int* lenvals= new int[1];
        fraglens[0]= 1;
        lenvals[0]= 200;
	DiscreteDF* fraglen_dist = new DiscreteDF(fraglens, lenvals, 1);

	DataFrame* f = new DataFrame(fraglen_dist, fragsta_cumu);

	Exon* e1 = new Exon(3185, 157);
	Exon* e2 = new Exon(3186, 35706);
	Exon* e3 = new Exon(3187, 87);
	Exon* e4 = new Exon(3188, 149822);
	f->addExon(e1);
	f->addExon(e2);
	f->addExon(e3);
	f->addExon(e4);

	Fragment* f0 = new Fragment(1, 1, 50);
	f0->left[0] = 3185;
	f0->right[0] = 3185;
	f->addData(f0);
	Fragment* f1 = new Fragment(1, 1, 10);
	f1->left[0] = 3185;
	f1->right[0] = 3186;
	f->addData(f1);
	Fragment* f2 = new Fragment(2, 1, 20);
	f2->left[0] = 3185;
	f2->left[1] = 3186;
	f2->right[0] = 3187;
	f->addData(f2);
	Fragment* f3 = new Fragment(1, 1, 20);
	f3->left[0] = 3186;
	f3->right[0] = 3186;
	f->addData(f3);
	Fragment* f4 = new Fragment(1, 1, 2);
	f4->left[0] = 3186;
	f4->right[0] = 3187;
	f->addData(f4);
	Fragment* f5 = new Fragment(1, 1, 100);
	f5->left[0] = 3187;
	f5->right[0] = 3187;
	f->addData(f5);
	Fragment* f6 = new Fragment(1, 1, 5);
	f6->left[0] = 3185;
	f6->right[0] = 3187;
	f->addData(f6);
	Fragment* f7 = new Fragment(1, 1, 10);
	f7->left[0] = 3185;
	f7->right[0] = 3188;
	f->addData(f7);

	vector<Exon*>* v1v = new vector<Exon*>();
	v1v->push_back(e1);
	v1v->push_back(e2);
	v1v->push_back(e3);

	vector<Exon*>* v2v = new vector<Exon*>();
	v2v->push_back(e1);
	v2v->push_back(e2);
	v2v->push_back(e3);
	v2v->push_back(e4);

	Variant* v1 = new Variant(v1v);
	v1->id = 1;

	Variant* v2 = new Variant(v2v);
	v2->id = 2;

	set<Variant*, VariantCmp>* varis = new set<Variant*, VariantCmp>();
	varis->insert(v1);
	varis->insert(v2);
	//varis->push_back(v3);

	debugdf(f);
	debugmodel(new Model(varis));

	fixUnexplFrags(varis, f);

	Model* model = new Model(varis);
	Casper* casp = new Casper(model, f);

	return casp;
}

int main() {
//int _tmain(int argc, _TCHAR* argv[])
	srand((unsigned)time( NULL ));

	Casper* c = david2();
	double* pi = c->calculateMode();
	//double inte = c->calculateIntegral();

	int n = c->model->count();

	for (int i = 0; i < n; i++)
	{
		printf("%f\n", pi[i]);
	}
	printf("\n");
	
	vector<Variant*>* allvariants = c->frame->allVariants();
	vector<Model*>* allmodels = c->frame->allModels();
	
	/*Model* center = allmodels->at(20);
	printf("%s, %i\n", center->getCodeStr(allvariants), center->count());
	map<Model*, double, ModelCmp> counts;
	SmartModelDist* dis = new SmartModelDist(new Seppel(c->frame), c->frame, center, 0.8);
	double nums = 0;
	for (int i = 0; i < 1000; i++)
	{
		Model* a = dis->sample();
		//printf("%i\n", a->count());
		if (a != 0)
		{
			double count = 0;
			if (counts.count(a) > 0)
			{
				count = counts[a];
			}
			counts[a] = count + 1;
			nums++;
		}
		else
		{
			double count = 0;
			if (counts.count(center) > 0)
			{
				count = counts[center];
			}
			counts[center] = count + 1;
			nums++;
		}
	}

	map<Model*, double, ModelCmp>::const_iterator xi;
	for (xi = counts.begin(); xi != counts.end(); xi++)
	{
		printf("%s@%i\t%f\t%f\n", xi->first->getCodeStr(allvariants), xi->first->count(), xi->second / nums, exp(dis->densityLn(xi->first)));
	}*/

	Seppel* sep1 = new Seppel(c->frame);
	sep1->exploreExact();
	map<Model*, double, ModelCmp> res1 = sep1->resultPPIntegral();

	printf("EXACT\n");

	Seppel* sep2 = new Seppel(c->frame);
	sep2->explorePrior(100000);
	map<Model*, double, ModelCmp> res2 = sep2->resultPPMCMC();

	printf("PRIOR\n");

	Seppel* sep3 = new Seppel(c->frame);
	sep3->exploreSmart(c->model, 100000);
	map<Model*, double, ModelCmp> res3 = sep3->resultPPMCMC();

	printf("SMART\n");
	
	map<Model*, double, ModelCmp>::const_iterator mi;
	for (mi = res1.begin(); mi != res1.end(); mi++)
	{
		Model* m = mi->first;
		if (res1[m] > 0.0001)
		{
			const char* code = m->getCodeStr(allvariants);
			printf("%s\t%f\t%f\t%f\n", code, res1[m], res2[m], res3[m]);
		}
	}

	return 0;
}
