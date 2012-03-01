#include "model_cmp.h"

size_t ModelCmp::operator()(const Model* a) const
{
	Model* ma = (Model*)a;
	return ma->hashcode;
}
bool ModelCmp::operator()(const Model* a, const Model* b) const
{
	Model* ma = (Model*)a;
	Model* mb = (Model*)b;

	return ma->compare(mb) < 0;
}