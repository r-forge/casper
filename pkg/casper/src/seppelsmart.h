#include "smartmodeldist.h"

class SeppelSmart
{
public:
	SeppelSmart(DataFrame* frame, Gene* gene);
	// does metropolis hastings given a start model, uses SeppelSmartDist as proposal
	hash_map<Model*, double, ModelCmp> SeppelSmart::calculate(Model* center);

	int runs;
	int burning;
	int thinning;
private:
	DataFrame* frame;
	Gene* gene;
};