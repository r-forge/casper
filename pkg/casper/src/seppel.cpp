#include "seppel.h"
#include <iostream>

Seppel::Seppel(DataFrame* frame, Gene* gene)
{
	this->frame = frame;
	this->gene = gene;
}

double Seppel::calcIntegral(Model* model)
{
	if (model == NULL)
	{
		return 1;
	}
	if (integrals.count(model) > 0)
	{
		return integrals[model];
	}

	double like = 1;

	Casper* casp = new Casper(model, frame);
	if (casp->isValid())
	{
		double* mode = casp->calculateMode();
		modes[model] = mode;
		like = casp->calculateIntegral(mode, model->count());
	}
	integrals[model] = like;

	return like;
}

const char* getmodelcode2(vector<Variant*>* allvariants, Model* model)
{
	int n = allvariants->size();
	char* str = new char[n + 1];
	str[n] = '\0';

	for (int i = 0; i < (int)allvariants->size(); i++)
	{
		if (model->contains(allvariants->at(i)))
		{
			str[i] = '1';
		}
		else
		{
			str[i] = '0';
		}
	}

	return str;
}

void Seppel::exploreExact()
{
	vector<Model*>* models = frame->allModels(gene);

	vector<Model*>::const_iterator mi;
	for (mi = models->begin(); mi != models->end(); mi++)
	{
		Model* model = *mi;
		calcIntegral(model);
	}
}
void Seppel::explorePrior(int runs)
{
	vector<Model*>* allmodels = frame->allModels(gene);
	vector<Model*>* models = new vector<Model*>();
	vector<Model*>::const_iterator ami;
	for (ami = allmodels->begin(); ami != allmodels->end(); ami++)
	{
		Casper* ncasp = new Casper(*ami, frame);
		if (ncasp->isValid())
		{
			models->push_back(ncasp->model);
			counts[ncasp->model] = 0;
		}
	}

	int onum = Casper::randi(models->size());
	Model* omodl = models->at(onum);
	double olike = calcIntegral(omodl);

	int accepted = 0;
	
	for (int r = 0; r < runs; r++)
	{
		int nnum = Casper::randi(models->size());
		Model* nmodl = models->at(nnum);

		double nlike = calcIntegral(nmodl);
		double l = nlike - olike;
		double p = exp(l);
		double x = Casper::randd();
		if (x <= p)
		{
			olike = nlike;
			omodl = nmodl;
			accepted++;
		}

		counts[omodl]++;
	}
}
void Seppel::exploreSmart(Model* startmodel, int runs)
{
	vector<Variant*>* allvars = frame->allVariants(gene);
    //FILE * pFile = fopen("C:\\prop_file.txt", "w");
    //FILE * vFile = fopen("C:\\visi_file.txt", "w");

	Model* omodl = startmodel;
	double olike = calcIntegral(omodl);
	SmartModelDist* odist = new SmartModelDist(this, omodl, 0.8);
	
	int accepted = 0;

	for (int r = 0; r < runs; r++)
	{
		Model* nmodl = odist->sample();
		double nlike = calcIntegral(nmodl);

		//fprintf(vFile, "%s\n", getmodelcode2(allvars, omodl));
		//fprintf(pFile, "%s\n", getmodelcode2(allvars, nmodl));

		if (nlike != 1)
		{
			SmartModelDist* ndist = new SmartModelDist(this, nmodl, 0.8);

			double nprob = odist->densityLn(nmodl);
			double oprob = ndist->densityLn(omodl);

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

		double count = 0;
		if (counts.count(omodl) > 0)
		{
			count = counts[omodl];
		}
		counts[omodl] = count + 1;	
	}

	//fclose(pFile);
	//fclose(vFile);
}

map<Model*, double*, ModelCmp> Seppel::resultModes()
{
	return modes;
}
map<Model*, double, ModelCmp> Seppel::resultPPIntegral()
{
	map<Model*, double, ModelCmp> probs;

	double imax = -DBL_MAX;
	
	map<Model*, double, ModelCmp>::const_iterator mi;
	for (mi = integrals.begin(); mi != integrals.end(); mi++)
	{
		if (mi->second == 1)
		{
			continue;
		}
		imax = max(mi->second, imax);
	}

	double asum = 0;
	for (mi = integrals.begin(); mi != integrals.end(); mi++)
	{
		if (mi->second == 1)
		{
			continue;
		}
		asum += exp(mi->second - imax);
	}
	double lsum = imax + log(asum);

	for (mi = integrals.begin(); mi != integrals.end(); mi++) 
	{
		if (mi->second == 1)
		{
			continue;
		}
		probs[mi->first] = exp(mi->second - lsum);
	}

	return probs;
}
map<Model*, double, ModelCmp> Seppel::resultPPMCMC()
{
	map<Model*, double, ModelCmp> probs;
	map<Model*, double, ModelCmp>::const_iterator mi;
	double total = 0;
	for (mi = counts.begin(); mi != counts.end(); mi++)
	{
		total = total + mi->second;
	}

	for (mi = counts.begin(); mi != counts.end(); mi++)
	{
		if (mi->second > 0)
		{
			probs[mi->first] = mi->second / total;
		}
	}
	return probs;
}

double* Seppel::normalizeIntegrals(double* values, int n)
{
	double imax = -DBL_MAX;
	
	for (int i = 0; i < n; i++)
	{
		imax = max(values[i], imax);
	}
	
	double asum = 0;
	for (int i = 0; i < n; i++)
	{
		asum += exp(values[i] - imax);
	}

	double lsum = imax + log(asum);
	
	double* probs = new double[n];
	double psum = 0;
	for (int i = 0; i < n; i++)
	{
		probs[i] = exp(values[i] - lsum);
		psum += probs[i];
	}
	for (int i = 0; i < n; i++)
	{
		probs[i] = probs[i] / psum;
	}
	return probs;
}