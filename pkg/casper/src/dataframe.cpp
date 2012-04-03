#include "dataframe.h"

DataFrame::DataFrame(DiscreteDF* fraglen_dist, double (*fragsta_cumu)(double))
{
	this->fraglen_dist = fraglen_dist;
	this->fragsta_cumu = fragsta_cumu;
	this->frag_readlen = 75;

        fraglen_minx= fraglen_dist->value(0);
        fraglen_maxx= fraglen_dist->value((fraglen_dist->size)-1);
	//	for (fraglen_minx = 0; fraglen_dist->probability(fraglen_minx) == 0; fraglen_minx++) ;
	//	for (fraglen_maxx = fraglen_dist->size - 1; fraglen_dist->probability(fraglen_maxx) == 0; fraglen_maxx--) ;
}

void DataFrame::addData(Fragment* f)
{
	this->data.push_back(f);
}
void DataFrame::addExon(Exon* e)
{
	e->num = this->exons.size();
	this->exons.push_back(e);
}

map<Fragment*, double> DataFrame::probabilities(Variant* v)
{
	if (this->cache.count(v) > 0)
	{
		return this->cache[v];
	}

	list<Fragment*>::const_iterator fi;
	for (fi = data.begin(); fi != data.end(); fi++)
	{
		Fragment* f = *fi;
		double p = probability(v, f);
		if (p > 0.0)
        {
			cache[v][f] = p;
		}
    }
	return cache[v];
}
double DataFrame::probability(Variant* v, Fragment* f)
{
	double p = 0.0;

	if (v->contains(f))
	{
		int fs = v->indexOf(f->left[0]);
		int fe = v->indexOf(f->left[f->leftc - 1]);
		int bs = v->indexOf(f->right[0]);
		int be = v->indexOf(f->right[f->rightc - 1]);

		p = prob(fs, fe, bs, be, v->positions, v->length);
	}
	return p;
}

double DataFrame::prob(int fs, int fe, int bs, int be, int* pos, double T)
{
	// lower bound for start of left transcript
	double a1 = max(pos[fs], pos[fe] - frag_readlen + 1);
	// upper bound for start of left transcript
	double b1 = min(pos[fs + 1] - 1, pos[fe + 1] - frag_readlen);
	// lower bound for start of right transcript
	double a2 = max(pos[bs] + frag_readlen, pos[be] + 1);
	// upper bound for start of right transcript
	double b2 = min(pos[bs + 1] + frag_readlen - 1, pos[be + 1]);

	double psum = 0;

        for (int i=0; i< fraglen_dist->size; i++) { //stop before T
	  double l= fraglen_dist->value(i);
	  double mb = 1.0 - l / T;
	  double rb = min(min(b1, b2 - l) / T, mb);
	  double lb = min((max(a1, a2 - l) - 1.0) / T, mb);

	  if (lb >= rb) { continue; }

	  double punc = (fragsta_cumu(rb) - fragsta_cumu(lb)) / fragsta_cumu(mb);

	  double factor = 0;
	  if (l <= T && punc > 0) {
	    factor = fraglen_dist->probability(i);
	    //			factor = fraglen_dist->probability((int)l);
	    if (T < fraglen_maxx) {
	      factor /= fraglen_dist->cumulativeProbability((int)(T-fraglen_minx));
	      //				factor /= fraglen_dist->cumulativeProbability((int)T);
	    }
	  }

	  psum += punc * factor;
	}
	return psum;
}

Variant* DataFrame::path2Variant(Fragment* f) 
{
	int eid; Exon *ex;
	vector<Exon*>::iterator itexon;
	map<int, Exon*> id2exon;
	for (itexon= exons.begin(); itexon != exons.end(); itexon++) {
		ex= (*itexon);
		id2exon[ex->id] = ex;
	}
	vector<Exon*>* el = new vector<Exon*>();
	for (itexon= exons.begin(); (*itexon)->id != f->left[0]; itexon++) {
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
		ex = exons[eid];
		el->push_back(ex);
	}
	while ((*itexon)->id != eid) { itexon++; }
	itexon++;
	while (itexon != exons.end()) {
		ex= (*itexon);
		el->push_back(ex);
		itexon++;
	}
	Variant* v = new Variant(el);
	return v;
}


int DataFrame::fixUnexplFrags(set<Variant*, VariantCmp>* initvars)
{
	// copy all fragments
	set<Fragment*>* queue = new set<Fragment*>(data.begin(), data.end());

	// remove the fragments from the queue we can explain with our variants
	set<Variant*, VariantCmp>::iterator vi;
	for (vi = initvars->begin(); vi != initvars->end(); vi++) 
	{
		// remove the fragments that this variant can explain from our queue
		map<Fragment*, double> probs = probabilities(*vi);
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

	int discarded = 0;

	// while we still have unexplained fragments
	while (queue->size() > 0)
	{
		// pop the first fragment
		Fragment* frag = *queue->begin();
		queue->erase(queue->begin());

		Variant* nv = path2Variant(frag);

		// check if the new variant can explain the fragment
		map<Fragment*, double> probs = probabilities(nv);
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
			data.remove(frag);
		}
	}

	return discarded;
}

void DataFrame::allVariantsRec(vector<Exon*>* stack, unsigned int level, vector<Variant*>* varis)
{
	if (exons.size() == level)
	{
		if (stack->size() > 0)
		{
			Variant* v = new Variant(stack);
			varis->push_back(v);
		}
		return;
	}
	
	stack->push_back(exons.at(level));
	allVariantsRec(stack, level + 1, varis);
	stack->pop_back();
	allVariantsRec(stack, level + 1, varis);
}
void DataFrame::allModelsRec(vector<Variant*>* stack, unsigned int level, vector<Variant*>* vars, vector<Model*>* models)
{
	if (vars->size() == level)
	{
		if (stack->size() > 0)
		{
			Model* m = new Model(stack);
			models->push_back(m);
		}
		return;
	}

	stack->push_back(vars->at(level));
	allModelsRec(stack, level + 1, vars, models);
	stack->pop_back();	
	allModelsRec(stack, level + 1, vars, models);
}
vector<Model*>* DataFrame::allModels()
{
	vector<Variant*>* varis = new vector<Variant*>();
	vector<Exon*>* estack = new vector<Exon*>();
	allVariantsRec(estack, 0, varis);
	
	vector<Variant*>* vstack = new vector<Variant*>();
	vector<Model*>* models = new vector<Model*>();
	allModelsRec(vstack, 0, varis, models);
	return models;
}
vector<Variant*>* DataFrame::allVariants()
{
	vector<Variant*>* varis = new vector<Variant*>();
	vector<Exon*>* estack = new vector<Exon*>();
	allVariantsRec(estack, 0, varis);
	return varis;
}

void DataFrame::debugprint()
{		
	// Exons
	printf("Exons:\n");
	vector<Exon*>::const_iterator ei;
	for (ei = exons.begin(); ei != exons.end(); ei++)
	{
		Exon* e = *ei;
		printf("%i\t%i\n", e->id, e->length);
	}
	printf("\n");

	// Fragments
	printf("Fragments:\n");
	list<Fragment*>::const_iterator fi;
	for (fi = data.begin(); fi != data.end(); fi++)
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