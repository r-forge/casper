#include "smartmodeldist.h"

const double SmartModelDist::exon_weight = 1.0;
const double SmartModelDist::create_prob = 0.5;

SmartModelDist::SmartModelDist(Model* center, Gene* gene, double exp_exons)
{
	this->center = center;
	this->gene = gene;
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

	pnull = 0;
	for (vi = varis.begin(); vi != varis.end(); vi++)
	{
		Variant* v = *vi;
		pnull += prob(v);
	}
	pnull *= pcreate;

	if (varis.size() == 1)
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
	for (unsigned int u = 0; u < gene->exons.size(); u++)
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
	for (unsigned int i = 0; i < gene->exons.size(); i++)
	{
		double k = exon_used[i];
		exon_prob[i] = explen * ((double)k + fact) / (sumexused + fact * gex);
	}
}

Variant* SmartModelDist::makevar()
{
	int rv = Casper::randi(center->items.size());
	Variant* var = center->items[rv];
	//bool dir = var->strand;
	Gene* gene = var->gene;

	vector<Exon*>* nex = new vector<Exon*>();
	for (unsigned int i = 0; i < gene->exons.size(); i++)
	{
		double pk = exon_prob[i];

		double r = Casper::randd();
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
	for (unsigned int i = 0; i < gene->exons.size(); i++)
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

Model* SmartModelDist::Sample()
{
	vector<Variant*>* newm = new vector<Variant*>();
	vector<Variant*>::const_iterator vi;

	double x = Casper::randd();
	if (x < pcreate)
	{
		Variant* v = makevar();
		if (v->exonCount == 0)
		{
			return center;
		}
		if (center->contains(v))
		{
			return NULL;
		}
		
		for (vi = center->items.begin(); vi != center->items.end(); vi++)
		{
			newm->push_back(*vi);
		}
		newm->push_back(v);
	}
	else
	{
		int ri = Casper::randi(varis.size());
		Variant* r = varis[ri];

		for (vi = center->items.begin(); vi != center->items.end(); vi++)
		{
			Variant* v = *vi;
			if (r->compare(v) != 0)
			{
				newm->push_back(v);
			}
		}
	}
	return new Model(newm);
}
double SmartModelDist::DensityLn(Model* model)
{
	double p = -log(1.0 - pnull);
	if (center->count() == model->count())
	{
		p += log(1.0 - pcreate) - log((double)varis.size());
	}
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

		p += log(pcreate * prob(var));
	}
	return p;
}
