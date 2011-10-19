#include "stdafx.h"
#include "casper.h"
#include <algorithm>
#include <time.h>
//#include <R.h>
//#include <Rinternals.h>

using namespace std;

const int Casper::frag_readlen = 75;
const double Casper::em_maxprec = 0.0000000000000001;
const int Casper::em_maxruns = 1000;

Casper::Casper()
{
	srand((int)time(NULL));
}
void Casper::addExon(Exon* e)
{
	this->exons[e->id] = e;
}
Exon* Casper::getExon(int exonid)
{
	return exons[exonid];
}
void Casper::addFragment(Fragment* f)
{
	this->fragments.push_back(f);
}
void Casper::addVariant(Variant* v)
{
	v->num = variants.size();
	this->variants.push_back(v);

	list<Fragment*>::const_iterator fi;
	for (fi = fragments.begin(); fi != fragments.end(); fi++)
	{
		Fragment* f = *fi;

		double p = 0.0;

		if (v->contains(f))
		{
			int fs = v->indexOf(f->left[0]);
			int fe = v->indexOf(f->left[f->leftc - 1]);
			int bs = v->indexOf(f->right[0]);
			int be = v->indexOf(f->right[f->rightc - 1]);

			p = prob(fs, fe, bs, be, v->positions, v->length);
			p = p;
		}

		if (p > 0.0)
		{
			mempprobs[f][v] = p;
			memvprobs[v][f] = p;
		}
	}
}

int Casper::varcount()
{
	return this->variants.size();
}

double Casper::psi(double x)
{
	//return fragsta_cumu(x);
	// return fragsta_cumu.Interpolate(x);
	return x;
}
double Casper::prob(int fs, int fe, int bs, int be, int* pos, double T)
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

	for (double l = fraglen_prob->MinimumX; l <= fraglen_prob->MaximumX; l++)
	{
		double mb = 1.0 - l / T;
		double rb = min(min(b1, b2 - l) / T, mb);
		double lb = min((max(a1, a2 - l) - 1.0) / T, mb);

		if (lb >= rb)
		{
			continue;
		}

		double punc = (psi(rb) - psi(lb)) / psi(mb);

		double factor = 0;
		if (l <= T && punc > 0)
		{
			factor = fraglen_prob->Get((int)l);
			if (T < fraglen_cumu->MaximumX)
			{
				factor /= fraglen_cumu->Get((int)T);
			}
		}

		psum += punc * factor;
	}
	return psum;
}
double* Casper::randomPi()
{
	int l = variants.size();

	double* pi = new double[l];
	double pisum = 0;
	for (int i = 0; i < l; i++)
	{
		pi[i] = (double)rand() / (double)(RAND_MAX + 1);
		pisum += pi[i];
	}
	for (int i = 0; i < l; i++)
	{
		pi[i] /= pisum;
	}
	return pi;
}
double* Casper::calculateEM(double* pi)
{
	int l = variants.size();

	for (int r = 0; r < em_maxruns; r++)
	{
		map<Fragment*, double> mem;
		double readcount = 0;
		
		map<Fragment*, map<Variant*, double> >::const_iterator fi;
		for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)
		{
			double sum = 0;

			map<Variant*, double>::const_iterator vi;
			for (vi = fi->second.begin(); vi != fi->second.end(); vi++)
			{
				int i = vi->first->num;
				sum += pi[i] * vi->second;
			}
			mem[fi->first] = sum;
			readcount += fi->first->count;
		}

		double* npi = new double[l];
		for (int d = 0; d < l; d++)
		{
			Variant* dv = variants[d];

			double nsum = 0;
			if (memvprobs.count(dv) > 0)
			{
				map<Fragment*, double> fragtable = memvprobs[dv];
				map<Fragment*, double>::const_iterator fi;
				for (fi = fragtable.begin(); fi != fragtable.end(); fi++)
				{
					nsum += (double)fi->first->count * fi->second / mem[fi->first];
				}
			}

			npi[d] = nsum / readcount * pi[d];
		}

		pi = npi;
	}

	return pi;
}
void Casper::setDistributions(DiscreteDF* fraglendist, double (*fragstacumu)(double))
{
	this->fraglen_prob = fraglendist;
	this->fraglen_cumu = fraglendist->GetCumulative();
	this->fragsta_cumu = fragstacumu;
}
