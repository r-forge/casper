#include "casper.h"

class SeppelExact
{
public:
	SeppelExact(DataFrame* frame, Gene* gene);
        // stores all possible models
	vector<Model*>* models;
	// calculates the posterior probability of all possible models
        void calculate();
	//map<Model*, double, ModelCmp> calculate();
	// stores posterior probabilities
	map<Model*, double, ModelCmp> posprob;
	// stores modes (estimated expression)
        map<Model*, double*, ModelCmp> mode;
private:
	DataFrame* frame;
	Gene* gene;
};
