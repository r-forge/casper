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
	map<Fragment*, double> probabilities(Variant* v);  //uses cache if available, otherwise fills cache and returns prob
	double probability(Variant* v, Fragment* f);  //does not use cache

	// returns a list of all possible models that could explain this gene (explicit calculation)
	vector<Model*>* allModels(Gene* gene);
	vector<Variant*>* allVariants(Gene* gene);
	int frag_readlen;
private:

	int fraglen_minx;
        int fraglen_maxx;
	DiscreteDF* fraglen_dist;
	double (*fragsta_cumu)(double x);
        map<Variant*, map<Fragment*, double>, VariantCmp > cache;

	double prob(int fs, int fe, int bs, int be, int* pos, double T);
	void allVariantsRec(vector<Exon*>* stack, unsigned int level, Gene* gene, vector<Variant*>* vars);
	void allModelsRec(vector<Variant*>* stack, unsigned int level, vector<Variant*>* vars, vector<Model*>* models);
};
