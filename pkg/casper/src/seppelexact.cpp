#include "seppelexact.h"
#define DBL_MAX 1.79769e+308;

SeppelExact::SeppelExact(DataFrame* frame, Gene* gene) {
	this->frame = frame;
	this->gene = gene;
	this->models= frame->allModels(gene);
}

//map<Model*, double, ModelCmp> SeppelExact::calculate()
void SeppelExact::calculate()
{
	map<Model*, double, ModelCmp> resProbsUnNorm;

	vector<Model*>::const_iterator mi;
	double imax = -DBL_MAX;
	for (mi = models->begin(); mi != models->end(); mi++)
	{
		Model* model = *mi;

		Casper* casp = new Casper(model, frame);
		if (casp->isValid())
		{
			double* mode= casp->calculateMode();
			double inte= casp->calculateIntegral(mode, model->count());

			//double inte = casp->calculateIntegral();
			resProbsUnNorm[model] = inte;
			resModes[model] = mode;
			imax = max(inte, imax);
		}
	}

	double asum = 0;
	map<Model*, double, ModelCmp>::const_iterator mvi;
	for (mvi = resProbsUnNorm.begin(); mvi != resProbsUnNorm.end(); mvi++)
	{
		asum += exp(mvi->second - imax);
	}
	double lsum = imax + log(asum);

	double ppbest=0;
	for (mvi = resProbsUnNorm.begin(); mvi != resProbsUnNorm.end(); mvi++) 
	{
		double pp= exp(mvi->second - lsum);
		resProbs[mvi->first] = pp;
	}
}