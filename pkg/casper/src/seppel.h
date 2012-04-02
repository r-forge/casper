#include "smartmodeldist.h"

class Seppel
{
public:
	Seppel(DataFrame* frame, Gene* gene);

	// calculate posterior probability of all possible models
	void exploreExact();
	//
	void explorePrior(int runs);
	// does metropolis hastings given a start model, uses SeppelSmartDist as proposal
	void exploreSmart(Model* startmodel, int runs);

	map<Model*, double*, ModelCmp> resultModes();
	map<Model*, double, ModelCmp> resultPPIntegral();
	map<Model*, double, ModelCmp> resultPPMCMC();

private:
	DataFrame* frame;
	Gene* gene;
	int runs;

	map<Model*, double, ModelCmp> counts;
	map<Model*, double, ModelCmp> integrals;
	map<Model*, double*, ModelCmp> modes;

	double calcIntegral(Model* model);
};
