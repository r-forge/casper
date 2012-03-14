#include "casper.h"

class SeppelExact
{
public:
	SeppelExact(DataFrame* frame, Gene* gene);

	vector<Model*>* models;  // stores all possible models

        void calculate(); 	// calculate posterior probability of all possible models
        void rmModels(double thre); // eliminate models with posprob<=thre from posprob and mode

	map<Model*, double, ModelCmp> posprob; 	// stores posterior probabilities
        Model* bestModel; //model with highest posterior prob

        map<Model*, double*, ModelCmp> mode; 	// stores modes (estimated expression)
        double* modeBest; //mode for best model
private:
	DataFrame* frame;
	Gene* gene;
};
