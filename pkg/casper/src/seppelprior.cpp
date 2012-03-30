#include "seppelprior.h"
#include <time.h>

SeppelPrior::SeppelPrior(DataFrame* frame, Gene* gene)
{
	this->frame = frame;
	this->gene = gene;

	this->runs = 10000;
}

void SeppelPrior::calculate()
{
	map<Model*, double, ModelCmp> result;

	vector<Model*>* allmodels = frame->allModels(gene);
	vector<Model*>* models = new vector<Model*>();
	vector<Model*>::const_iterator ami;
	for (ami = allmodels->begin(); ami != allmodels->end(); ami++)
	{
		Casper* ncasp = new Casper(*ami, frame);
		if (ncasp->isValid())
		{
			models->push_back(ncasp->model);
			result[ncasp->model] = 0;
		}
	}

	printf("OK HERE\n");

	int onum = Casper::randi(models->size());
	printf("picked %i out of %i:\n", onum, models->size());
	Model* omodl = models->at(onum);
	Casper* ocasp = new Casper(omodl, frame);
	double olike = ocasp->calculateIntegral();

	map<Model*, double*, ModelCmp> modes;

	double rescounts = 0;
	int accepted = 0;
	
	for (int r = 0; r < runs; r++)
	{
		int nnum = Casper::randi(models->size());
		printf("picked %i out of %i:\n", nnum, models->size());
		Model* nmodl = models->at(nnum);
		Casper* ncasp = new Casper(nmodl, frame);

		double * nmode = ncasp->calculateMode();
		modes[nmodl] = nmode;
		double nlike = ncasp->calculateIntegral(nmode, nmodl->count());
		
		double l = nlike - olike;
		double p = exp(l);
		double x = Casper::randd();
		if (x <= p)
		{
			olike = nlike;
			omodl = nmodl;
			accepted++;
		}

		result[omodl]++;
		rescounts++;
	}
	printf("ALL OK");

	map<Model*, double, ModelCmp>::const_iterator mvi;
	for (mvi = result.begin(); mvi != result.end(); mvi++)
	{
		if (mvi->second > 0)
		{
			resProbs[mvi->first] = mvi->second / rescounts;
			resModes[mvi->first] = modes[mvi->first];
		}
	}
}
