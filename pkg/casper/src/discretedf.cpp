#include "discretedf.h"

DiscreteDF::DiscreteDF(double* data, int size)
{
	this->size = size;
	this->prob = new double[size];
	this->cumu = new double[size];

	double sum = 0;
	for (int i = 0; i < size; i++)
	{
		sum = sum + data[i];
		prob[i] = data[i];
		cumu[i] = sum;
	}
}
double DiscreteDF::probability(int i)
{
	return prob[i];
}
double DiscreteDF::cumulativeProbability(int i)
{
	return cumu[i];
}