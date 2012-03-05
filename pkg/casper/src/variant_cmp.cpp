#include "variant_cmp.h"

size_t VariantCmp::operator()(const Variant* a) const
{
	Variant* va = (Variant*)a;
	return va->hashcode;
}
bool VariantCmp::operator()(const Variant* a, const Variant* b) const
{
	Variant* va = (Variant*)a;
	Variant* vb = (Variant*)b;

	return va->compare(vb) < 0;
}
