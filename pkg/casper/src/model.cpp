#include "model.h"

using namespace std;

Model::Model(vector<Variant*>* variants)
{
	vector<Variant*>::const_iterator vi;
	for (vi = variants->begin(); vi != variants->end(); vi++)
	{
		Variant* v = *vi;
		int i = this->items.size();
		this->items.push_back(v);
		this->idmap[v] = i;
	}

	hashcode = gethash();
}
Model::Model(set<Variant*, VariantCmp>* variants)
{	
	set<Variant*, VariantCmp>::const_iterator vi;
	for (vi = variants->begin(); vi != variants->end(); vi++)
	{
		Variant* v = *vi;
		int i = this->items.size();
		this->items.push_back(v);
		this->idmap[v] = i;
	}

	hashcode = gethash();
}

bool Model::contains(Variant* v)
{
	return (this->idmap.count(v) > 0);
}
int Model::count()
{
	return this->items.size();
}
Variant* Model::get(int i)
{
	return this->items.at(i);
}
int Model::indexOf(Variant* v)
{
	return this->idmap[v];
}

int Model::compare(Model* other)
{
	if (this->count() < other->count())
    {
		return -1;
	}
	if (this->count() > other->count())
    {
		return +1;
	}

	map<Variant*, int, VariantCmp>::const_iterator ti = this->idmap.begin();
	map<Variant*, int, VariantCmp>::const_iterator oi = other->idmap.begin();
	while (ti != idmap.end())
	{
		int c = ti->first->compare(oi->first);
		if (c != 0)
		{
			return c;
		}

		ti++;
		oi++;
	}

	return 0;
}
int Model::gethash()
{
	int h = 0;
	map<Variant*, int, VariantCmp>::const_iterator ti;
	for (ti = idmap.begin(); ti != idmap.end(); ti++)
	{
		h ^= ti->first->hashcode;
	}
	return h;
}
