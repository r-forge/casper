#include "seppelexact.h"

class SeppelPrior
{
public:
	SeppelPrior(DataFrame* frame, Gene* gene);
	map<Model*, double, ModelCmp> calculate();

	int runs;

private:
	DataFrame* frame;
	Gene* gene;
};
