#include "seppelexact.h"

SeppelExact::SeppelExact(DataFrame* frame, Gene* gene)
{
	this->frame = frame;
	this->gene = gene;
}

unordered_map<Model*, double, ModelCmp> SeppelExact::calculate()
{
	unordered_map<Model*, double, ModelCmp> result;

	vector<Model*>* models = frame->allModels(gene);

	vector<Model*>::const_iterator mi;
	double imax = -DBL_MAX;
	for (mi = models->begin(); mi != models->end(); mi++)
	{
		Model* model = *mi;

		Casper* casp = new Casper(model, frame);
		if (casp->isValid())
		{
			double inte = casp->calculateIntegral();
			result[model] = inte;
			imax = max(inte, imax);
		}
	}

	double asum = 0;
	unordered_map<Model*, double, ModelCmp>::const_iterator mvi;
	for (mvi = result.begin(); mvi != result.end(); mvi++)
	{
		asum += exp(mvi->second - imax);
	}
	double lsum = imax + log(asum);

	unordered_map<Model*, double, ModelCmp> normalized;
	for (mvi = result.begin(); mvi != result.end(); mvi++)
	{
		normalized[mvi->first] = exp(mvi->second - lsum);
	}

	return normalized;
}