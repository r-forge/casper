#include "seppel.h"

const double SmartModelDist::exon_weight = 1.0;
const double SmartModelDist::create_prob = 0.5;

SmartModelDist::SmartModelDist(Seppel* seppel, Model* center, double exp_exons)
{
	this->seppel = seppel;
	this->center = center;
	this->gene = seppel->gene;
	this->exp_exons = exp_exons;

	vector<Variant*>::const_iterator vi;
	for (vi = center->items.begin(); vi != center->items.end(); vi++)
	{
		Variant* v = *vi;
		if (v->gene != gene)
		{
			continue;
		}
		varis.push_back(v);
	}

	updatepks();
	buildrmtable();

	pnull = 0;
	for (vi = varis.begin(); vi != varis.end(); vi++)
	{
		Variant* v = *vi;
		pnull += prob(v);
	}
	pnull += prob(0);

	if (varis.size() == 1)
	{
		pcreate = 1;
	}
	else if (gene->exons.size() == (int)(log((center->count() + 1.0) / log(2.0)) + 0.5))
	{
		pcreate = 0;
	}
	else
	{
		pcreate = create_prob;
	}
}

void SmartModelDist::updatepks()
{
	int maxexused = 0;
	int sumexused = 0;
	exon_used = new int[gene->exons.size()];
	for (int u = 0; u < gene->exons.size(); u++)
	{
		exon_used[u] = 0;
	}
	
	vector<Variant*>::const_iterator vi;

	for (vi = varis.begin(); vi != varis.end(); vi++)
	{
		Variant* v = *vi;

		for (int ei = 0; ei < v->exonCount; ei++)
		{
			Exon* e = gene->exons[ei];
			int i = gene->indexOf(e);
			exon_used[i]++;

			if (exon_used[i] > maxexused)
			{
				maxexused = exon_used[i];
			}
		}
		sumexused += v->exonCount;
	}

	double gex = gene->exons.size();
	double explen = exp_exons;
	if (explen > gex * 0.9)
	{
		explen = gex * 0.9;
	}
	double esti = (explen * maxexused - sumexused) / (gex - explen);
	double fact = max(esti + 0.001, exon_weight);

	exon_prob = new double[gene->exons.size()];
	for (int i = 0; i < gene->exons.size(); i++)
	{
		double k = exon_used[i];
		exon_prob[i] = explen * ((double)k + fact) / (sumexused + fact * gex);
	}
}

void SmartModelDist::buildrmtable()
{
	Model** possible = new Model*[varis.size()];
	double* integrals = new double[varis.size()];
	int n = 0;

	list<Variant*>* copy = new list<Variant*>(varis.begin(), varis.end());
	for (int i = 0; i < varis.size(); i++)
	{
		Variant* v = copy->front();
		copy->pop_front();
		Model* m = new Model(copy);
		copy->push_back(v);

		double like = seppel->calcIntegral(m);
		if (like != 2)
		{
			possible[n] = m;
			integrals[n] = like;
			n++;
		}
	}

	if (n == 0)
	{
		return;
	}

	//double* probs = Seppel::normalizeIntegrals(integrals, n);
	double* probs = new double[n];
	for (int i = 0; i < n; i++)
	{
		probs[i] = (double)1 / (double)n;
		Model* m = possible[i];
		removeprobs[m] = probs[i];
	}
}

Variant* SmartModelDist::makevar()
{
	int rv = runifdisc(0, center->items.size() - 1);
	Variant* var = center->items[rv];
	//bool dir = var->strand;
	Gene* gene = var->gene;

	vector<Exon*>* nex = new vector<Exon*>();
	for (int i = 0; i < gene->exons.size(); i++)
	{
		double pk = exon_prob[i];

		double r = runif();
		if (r < pk)
		{
			nex->push_back(gene->exons[i]);
		}
	}

	Variant* nva = new Variant(gene, nex);
	return nva;
}
double SmartModelDist::prob(Variant* v)
{
	double p = 1;
	for (int i = 0; i < gene->exons.size(); i++)
	{
		double pk = exon_prob[i];
		if (v != NULL && v->contains(gene->exons[i]))
		{
			p *= pk;
		}
		else
		{
			p *= 1.0 - pk;
		}
	}
	return p;
}

Model* SmartModelDist::sample()
{
	vector<Variant*>* newm = new vector<Variant*>();
	vector<Variant*>::const_iterator vi;

	double x = runif();
	if (x < pcreate)
	{
		Variant* v;
		do
		{
			v = makevar();
		}
		while (v->exonCount == 0 || center->contains(v));
		
		for (vi = center->items.begin(); vi != center->items.end(); vi++)
		{
			newm->push_back(*vi);
		}
		newm->push_back(v);
		return new Model(newm);
	}
	else
	{
		double r = runif();
		double psum = 0;
		map<Model*, double, ModelCmp>::const_iterator mi;
		Model* last = NULL;
		for (mi = removeprobs.begin(); mi != removeprobs.end(); mi++)
		{
			psum = mi->second + psum;
			if (r <= psum)
			{
				return mi->first;
			}
			last = mi->first;
		}
		return last;

		/*int ri = runifdisc(0, varis.size() - 1);
		Variant* r = varis[ri];

		for (vi = center->items.begin(); vi != center->items.end(); vi++)
		{
			Variant* v = *vi;
			if (r->compare(v) != 0)
			{
				newm->push_back(v);
			}
		}
		return new Model(newm);*/
	}
}
double SmartModelDist::densityLn(Model* model)
{
	if (model->count() < center->count())
	{
		//return log(1.0 - pcreate) - log((double)center->count());
		return log(1.0 - pcreate) + log(removeprobs[model]);
	}
	/*else if (center->count() == model->count())
	{
		return log(pcreate) + log(
	}*/
	else
	{
		Variant* var = NULL;
		vector<Variant*>::const_iterator vi;
		for (vi = model->items.begin(); vi != model->items.end(); vi++)
		{
			Variant* v = *vi;
			if (center->contains(v) == false)
			{
				var = v;
				break;
			}
		}

		return log(pcreate) + log(prob(var)) - log(1.0 - pnull);
	}
	return DBL_MIN;
}
