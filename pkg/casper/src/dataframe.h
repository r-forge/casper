#include "model_cmp.h"
#include "discretedf.h"
#include <list>
#include <vector>
using namespace std;

class DataFrame
{
public:
	// all exons mapped by their id
	map<int, Exon*> exons;
	// all genes mapped by their id
	map<int, Gene*> genes;
	// all fragments
	list<Fragment*> data;

	DataFrame(DiscreteDF* fraglen_dist, double (*fragsta_cumu)(double));

	void addData(Fragment* f);
	void addExon(Exon* e);
	void addGene(Gene* g);

	// probabilities of all the fragments given a variant
	map<Fragment*, double> probabilities(Variant* v);
	double probability(Variant* v, Fragment* f);

	// returns a list of all possible models that could explain this gene (explicit calculation)
	vector<Model*>* allModels(Gene* gene);
private:
	static const int frag_readlen;

	int fraglen_minx;
    int fraglen_maxx;
	DiscreteDF* fraglen_dist;
	double (*fragsta_cumu)(double x);
    map<Variant*, map<Fragment*, double>, VariantCmp > cache;

	double prob(int fs, int fe, int bs, int be, int* pos, double T);
	void allVariantsRec(vector<Exon*>* stack, int level, Gene* gene, vector<Variant*>* vars);
	void allModelsRec(vector<Variant*>* stack, int level, vector<Variant*>* vars, vector<Model*>* models);
};
