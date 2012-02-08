#include "casper.h"

class SeppelExact
{
public:
	SeppelExact(DataFrame* frame, Gene* gene);
	// calculates the probabilities of all possible explicit models through their integrals
	unordered_map<Model*, double, ModelCmp> calculate();
private:
	DataFrame* frame;
	Gene* gene;
};