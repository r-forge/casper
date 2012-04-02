#include "variant_cmp.h"
#include <map>
#include <set>
#include <vector>

using namespace std;

class Model
{
public:
	Model(vector<Variant*>* variants);
	Model(set<Variant*, VariantCmp>* variants);

	// amount of variants in this model
	int count();
	// does this model contain the variant v
	bool contains(Variant* v);
	// index of the variant in the variant list
	int indexOf(Variant* v);
	// get the ith variant
	Variant* get(int i);

	// list of variants
	vector<Variant*> items;

	// hashcode for model
	int hashcode;

	char* toString();

	// compares two models. returns 0 if equal. -1 and +1 used for sorting in map
	int compare(Model* other);
	// hash for this model
	int gethash();

	map<Variant*, int, VariantCmp> idmap;
private:
};
