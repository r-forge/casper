#include "casper.h"

class Seppel
{
public:
	Seppel(DataFrame* frame, Gene* gene);

	double calcIntegral(Model* model);

	// calculate posterior probability of all possible models
	void exploreExact();
	//
	void explorePrior(int runs);
	// does metropolis hastings given a start model, uses SeppelSmartDist as proposal
	void exploreSmart(Model* startmodel, int runs);

	map<Model*, double*, ModelCmp> resultModes();
	map<Model*, double, ModelCmp> resultPPIntegral();
	map<Model*, double, ModelCmp> resultPPMCMC();

	static double* normalizeIntegrals(double* values, int n);

	Gene* gene;

private:
	DataFrame* frame;

	map<Model*, double, ModelCmp> counts;
	map<Model*, double, ModelCmp> integrals;
	map<Model*, double*, ModelCmp> modes;
};

// moved header of SmartModelDist to here because Seppel and SmartModelDist now use eachother
class SmartModelDist
{
public:
	SmartModelDist(Seppel* seppel, Model* center, double exp_exons);
	
	// sample a proposal
	Model* sample();
	// densitiy of the proposal
	double densityLn(Model* model);

private:
	// extra weight of an exon if its already used
	static const double exon_weight;
	// probability to create, delete is 1-pcreate
    static const double create_prob;
	
	Seppel* seppel;
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
	
	map<Model*, double, ModelCmp> removeprobs;

	void updatepks();
	void buildrmtable();
	Variant* makevar();
	double prob(Variant* v);
};
