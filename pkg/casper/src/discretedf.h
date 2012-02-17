class DiscreteDF
{
public:
	int MinimumX;
	int MaximumX;

	DiscreteDF(int size);

	double Get(int i);
	void Set(int i, double value);

	DiscreteDF* GetCumulative();

private:
	double* data;
	int size;
};
