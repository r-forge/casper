#include "seppelprior.h"

class SmartModelDist
{
public:
	SmartModelDist(Model* center, Gene* gene, double exp_exons);
	
	// sample a proposal
	Model* Sample();
	// densitiy of the proposal
	double DensityLn(Model* model);

private:
	// extra weight of an exon if its already used
	static const double exon_weight;
	// probability to create, delete is 1-pcreate
    static const double create_prob;

    Model* center;
    vector<Variant*> varis;
    Gene* gene;

    // expected number of exons for a gene of this length
    double exp_exons;

    // number of times each exon of the gene was used in a variant
    int* exon_used;

    // probability of an exon with a given used count to be in a variant
    double* exon_prob;

    double pnull;
    double pcreate;

	void updatepks();
	Variant* makevar();
	double prob(Variant* v);
};
