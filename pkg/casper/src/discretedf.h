
// discrete density function for values 0 to size
class DiscreteDF
{
public:
	DiscreteDF(double* data, int size);

	int size;
	double probability(int i);
	double cumulativeProbability(int i);

private:
	double* prob;
	double* cumu;
};
