#include "discretedf.h"
#include "cppmemory.h"

DiscreteDF::DiscreteDF(double* data, int* vals, int size)
{
	this->size = size;
        this->values= new int[size];
	this->prob = new double[size];
	this->cumu = new double[size];

	double sum = 0;
	for (int i = 0; i < size; i++)
	{
	  values[i]= vals[i];
	  sum = sum + data[i];
	  prob[i] = data[i];
	  cumu[i] = sum;
	}
}

DiscreteDF::~DiscreteDF() {
  zaparray(values);   //delete [] values;
  zaparray(prob);   //delete [] prob;
  zaparray(cumu);   //delete [] cumu;
}

int DiscreteDF::value(int i) { return values[i]; }
double DiscreteDF::probability(int i) { return prob[i]; }
double DiscreteDF::cumulativeProbability(int i) { return cumu[i]; }
