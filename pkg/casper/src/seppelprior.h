#include "seppelexact.h"

class SeppelPrior
{
public:
	SeppelPrior(DataFrame* frame, Gene* gene);
	void calculate();
	
	map<Model*, double, ModelCmp> resProbs; 	// stores posterior probabilities
	map<Model*, double*, ModelCmp> resModes; 	// stores modes (estimated expression)

	int runs;

private:
	DataFrame* frame;
	Gene* gene;
};
