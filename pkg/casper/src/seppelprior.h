#include "seppelexact.h"

class SeppelPrior
{
public:
	SeppelPrior(DataFrame* frame, Gene* gene);
	unordered_map<Model*, double, ModelCmp, ModelCmp> calculate();

	int runs;

private:
	DataFrame* frame;
	Gene* gene;
};