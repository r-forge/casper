#include "casper.h"

class SeppelExact
{
public:
	SeppelExact(DataFrame* frame, Gene* gene);

	vector<Model*>* models;  // stores all possible models

	void calculate(); 	// calculate posterior probability of all possible models

	map<Model*, double, ModelCmp> resProbs; 	// stores posterior probabilities
	map<Model*, double*, ModelCmp> resModes; 	// stores modes (estimated expression)

private:
	DataFrame* frame;
	Gene* gene;
};
