#include "seppelexact.h"

class SeppelPrior
{
public:
	SeppelPrior(DataFrame* frame, Gene* gene);
	hash_map<Model*, double, ModelCmp> calculate();

	int runs;

private:
	DataFrame* frame;
	Gene* gene;
};