#include "variant_cmp.h"
#include <map>
#include <vector>

using namespace std;

class Model
{
public:
	Model(vector<Variant*>* list);

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
	// hashcode
	int hashcode;

	// compares two models. returns 0 if equal. -1 and +1 used for sorting in map
	int compare(Model* other);
	// hash for this model
	int gethash();

private:
	map<Variant*, int, VariantCmp> idmap;
};