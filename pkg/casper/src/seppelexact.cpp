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
  map<Model*, double, ModelCmp> posprobUnNorm;

  vector<Model*>::const_iterator mi;
  double imax = -DBL_MAX;
	for (mi = models->begin(); mi != models->end(); mi++)
	{
		Model* model = *mi;

		Casper* casp = new Casper(model, frame);
		if (casp->isValid())
		{
		  mode[model]= casp->calculateMode();
		  double inte= casp->calculateIntegral(mode[model], model->count());
		  //double inte = casp->calculateIntegral();
		  posprobUnNorm[model] = inte;
		  imax = max(inte, imax);
		}
	}

	double asum = 0;
	map<Model*, double, ModelCmp>::const_iterator mvi;
	for (mvi = posprobUnNorm.begin(); mvi != posprobUnNorm.end(); mvi++)
	{
		asum += exp(mvi->second - imax);
	}
	double lsum = imax + log(asum);

	double ppbest=0;
	for (mvi = posprobUnNorm.begin(); mvi != posprobUnNorm.end(); mvi++) {
	  double pp= exp(mvi->second - lsum);
	  posprob[mvi->first] = pp;
	  if (pp>ppbest) {
	    bestModel= mvi->first;
	    ppbest=pp;
	  }
	}
}

void SeppelExact::rmModels(double thre) {
  map<Model*, double, ModelCmp>::iterator mvi;
  for (mvi = posprob.begin(); mvi != posprob.end(); mvi++) {
    double p= posprob[mvi->first];
    if (!(p>=thre)) {
      posprob.erase(mvi);
      mode.erase(mvi->first);
    }
  }
}
