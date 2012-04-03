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
Model::Model(list<Variant*>* variants)
{	
	list<Variant*>::const_iterator vi;
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

char* Model::toString()
{
	char* str = new char[3000];
	str[0] = '\0';

	map<Variant*, int, VariantCmp>::const_iterator vi;
	for (vi = idmap.begin(); vi != idmap.end(); vi++)
	{
		Variant* v = vi->first;
		sprintf(str, "%s {%s}", str, v->toString());
	}

	return str;
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
		h += h * 7 + ti->first->hashcode;
	}
	return h;
}

const char* Model::getCodeStr(vector<Variant*>* allvariants)
{
	int n = allvariants->size();
	char* str = new char[n + 1];
	str[n] = '\0';

	for (int i = 0; i < (int)allvariants->size(); i++)
	{
		if (this->contains(allvariants->at(i)))
		{
			str[i] = '1';
		}
		else
		{
			str[i] = '0';
		}
	}

	return str;
}