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
	this->exons[e->id] = e;
}
void DataFrame::addGene(Gene* g)
{
	this->genes[g->id] = g;
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

void DataFrame::allVariantsRec(vector<Exon*>* stack, unsigned int level, Gene* gene, vector<Variant*>* varis)
{
	if (gene->exons.size() == level)
	{
		if (stack->size() > 0)
		{
			Variant* v = new Variant(gene, stack);
			if (probabilities(v).size() > 0)
			{
				varis->push_back(v);
			}
		}
		return;
	}
	
	stack->push_back(gene->exons.at(level));
	allVariantsRec(stack, level + 1, gene, varis);
	stack->pop_back();
	allVariantsRec(stack, level + 1, gene, varis);
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
vector<Model*>* DataFrame::allModels(Gene* gene)
{
	vector<Variant*>* varis = new vector<Variant*>();
	vector<Exon*>* estack = new vector<Exon*>();
	allVariantsRec(estack, 0, gene, varis);
	
	vector<Variant*>* vstack = new vector<Variant*>();
	vector<Model*>* models = new vector<Model*>();
	allModelsRec(vstack, 0, varis, models);
	return models;
}
vector<Variant*>* DataFrame::allVariants(Gene* gene)
{
	vector<Variant*>* varis = new vector<Variant*>();
	vector<Exon*>* estack = new vector<Exon*>();
	allVariantsRec(estack, 0, gene, varis);
	return varis;
}
