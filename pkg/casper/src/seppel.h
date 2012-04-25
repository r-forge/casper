#include "casper.h"

class Seppel
{
public:
	Seppel(DataFrame* frame);

	double calcIntegral(Model* model);

	void exploreExact(); // exhaustive enumeration of all possible models
	void exploreUnif(int runs); //Metropolis-Hastings MCMC with independent proposals (uniform)
	void exploreSmart(Model* startmodel, int runs); // Metropolis-Hastings MCMC with random walk (uses SeppelSmartDist as proposal)

	map<Model*, double*, ModelCmp> resultModes();
	map<Model*, double, ModelCmp> resultPPIntegral();  //compute post prob using marginal likelihoods. Sets integralSum & integralMax
	map<Model*, double, ModelCmp> resultPPMCMC();

	static double* normalizeIntegrals(double* values, int n);
	double integralSum; //sum integrals/exp(integralMax)
	double integralMax; //maximum log(integrals), i.e. log(marginal likelihood) + log(prior)

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
	SmartModelDist(Seppel* seppel, DataFrame* frame, Model* center, double exp_exons);
	
	// sample a proposal
	Model* sample();
	// densitiy of the proposal
	double densityLn(Model* model);

private:
	// extra weight of an exon if its already used
	static const double exon_weight;
	// probability to create, delete is 1-pcreate
    static const double create_prob;
	
	DataFrame* frame;
	Seppel* seppel;
    Model* center;

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
