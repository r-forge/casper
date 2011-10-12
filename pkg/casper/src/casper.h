#include <map>
#include <list>
#include <vector>
#include "discretedf.h"
#include "fragment.h"
#include "variant.h"

using namespace std;

class Casper
{
public:
	Casper();

	void addVariant(Variant* v);
	void addExon(Exon* e);
	Exon* addExon(int id, int start, int end);
	Exon* getExon(int id);
	void addFragment(Fragment* f);
	void removeVariant(Variant* v);
	void setDistributions(DiscreteDF* fraglendist, double(*fragstadist)(double));

	int varcount();
	
	double* randomPi();
	double* calculateEM(double* pi);
	double* calcualteGibs(double* pi);
	double* calculateMH(double* pi);
	double* calcualteLogLikelihood(double* pi);
	
	map<int, Exon*> exons;
	list<Fragment*> fragments;
	vector<Variant*> variants;

	map<Variant*, map<Fragment*, double> > memvprobs;
private:

	map<Fragment*, map<Variant*, double> > mempprobs;
	
	DiscreteDF* fraglen_prob;
	DiscreteDF* fraglen_cumu;
	double (*fragsta_cumu)(double x);

	static const int frag_readlen;
	static const int em_maxruns;
	static const double em_maxprec;

	double psi(double x);
	double prob(int fs, int fe, int bs, int be, int* pos, double T);
	double randd();
};
