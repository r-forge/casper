#include "smartmodeldist.h"

class SeppelSmart
{
public:
	SeppelSmart(DataFrame* frame, Gene* gene, Model* startmodel);

	// does metropolis hastings given a start model, uses SeppelSmartDist as proposal
	void calculate(); 	// calculate posterior probability of all possible models

	map<Model*, double, ModelCmp> resProbs; 	// stores posterior probabilities
	map<Model*, double*, ModelCmp> resModes; 	// stores modes (estimated expression)

	int runs;
	int burning;
	int thinning;

private:
	DataFrame* frame;
	Gene* gene;
	Model* startmodel;
};
