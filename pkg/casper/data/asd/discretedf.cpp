#include "stdafx.h"
#include "discretedf.h"

DiscreteDF::DiscreteDF(int size)
{
	this->size = size;
	this->data = new double[size];
	this->MaximumX = -1;
	this->MinimumX = size;
}
double DiscreteDF::Get(int i)
{
	return data[i];
}
void DiscreteDF::Set(int i, double value)
{
	data[i] = value;

	if (i < MinimumX)
	{
		MinimumX = i;
	}
	if (i > MaximumX)
	{
		MaximumX = i;
	}
}
DiscreteDF* DiscreteDF::GetCumulative()
{
	DiscreteDF* cumu = new DiscreteDF(this->size);

	double sum = 0;
	for (int i = 0; i < this->size; i++)
	{
		sum = sum + data[i];
		cumu->Set(i, sum);
	}
	cumu->Set(this->size - 1, 1.0);

	return cumu;
}