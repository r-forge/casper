#include "seppelsmart.h"

SeppelSmart::SeppelSmart(DataFrame* frame, Gene* gene)
{
	this->frame = frame;
	this->gene = gene;

	this->runs = 10000;
	this->burning = 0;
	this->thinning = 1;
}

hash_map<Model*, double, ModelCmp> SeppelSmart::calculate(Model* center)
{
	Model* omodl = center;
	Casper* ocasp = new Casper(omodl, frame);
	double olike = ocasp->calculateIntegral();
	SmartModelDist* odist = new SmartModelDist(omodl, gene, 0.8);
	
	hash_map<Model*, double, ModelCmp> result;
	double rescounts = 0;
	int accepted = 0;

	for (int r = 0; r < runs; r++)
	{
		Model* nmodl = odist->Sample();
		double nlike = 0;
		if (nmodl != NULL)
		{
			Casper* ncasp = new Casper(nmodl, frame);
			if (ncasp->isValid())
			{
				nlike = ncasp->calculateIntegral();
			}
		}

		if (nlike != 0)
		{
			SmartModelDist* ndist = new SmartModelDist(nmodl, gene, 0.8);

			double nprob = odist->DensityLn(nmodl);
			double oprob = ndist->DensityLn(omodl);

			double l = nlike - olike + oprob - nprob;
			double lp = exp(l);
			double x = Casper::randd();
			if (x < lp)
			{
				omodl = nmodl;
				odist = ndist;
				olike = nlike;
				accepted++;
			}
		}

		if (r >= burning && r % thinning == 0)
		{
			double count = 0;
			if (result.count(omodl) > 0)
			{
				count = result[omodl];
			}
			result[omodl] = count + 1;
			rescounts++;
		}
	}
	
	hash_map<Model*, double, ModelCmp> normalized;
	
	hash_map<Model*, double, ModelCmp>::const_iterator mvi;
	for (mvi = result.begin(); mvi != result.end(); mvi++)
	{
		if (mvi->second > 0)
		{
			normalized[mvi->first] = mvi->second / rescounts;
		}
	}

	return normalized;
}