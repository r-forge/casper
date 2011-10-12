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
	void addFragment(Fragment* f);
	Fragment* addFragment(int leftc, int rightc, int count);
	void removeVariant(Variant* v);
	void setDistributions(DiscreteDF* fraglendist);

	int varcount();
	
	double* randomPi();
	double* calculateEM(double* pi);
	double* calcualteGibs(double* pi);
	double* calculateMH(double* pi);
	double* calcualteLogLikelihood(double* pi);
private:
	map<int, Exon*> exons;
	list<Fragment*> fragments;
	vector<Variant*> variants;

	map<Fragment*, map<Variant*, double> > mempprobs;
	map<Variant*, map<Fragment*, double> > memvprobs;

	DiscreteDF* fraglen_prob;
	DiscreteDF* fraglen_cumu;

	static const int frag_readlen;
	static const int em_maxruns;
	static const double em_maxprec;

	double psi(double x, double T);
	double prob(int fs, int fe, int bs, int be, int* pos, double T);
	double randd();
};
