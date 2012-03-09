#include "variant.h"

Variant::Variant(Gene* gene, vector<Exon*>* exons)
{
	this->id = -1;
	this->gene = gene;
	//this->strand = strand;

	this->exonCount = exons->size();
	this->exons = new Exon*[exonCount];
	this->positions = new int[exonCount + 1];
	this->positions[0] = 1;
		
    this->codelen = (exonCount + 31) / 32;
    this->codes = new int[codelen];
	for (int i = 0; i < codelen; i++)
	{
		this->codes[i] = 0;
	}

	int i = 0;
	vector<Exon*>::const_iterator ei;
	for (ei = exons->begin(); ei != exons->end(); ei++)
	{
		Exon* exon = *ei;
		
		this->idmap[exon->id] = i;
		this->exons[i] = exon;
		this->positions[i + 1] = this->positions[i] + exon->length;
		
        int j = gene->indexOf(exon);
        this->codes[j / 32] |= 1 << (j % 32);
		i++;
	}

	this->length = this->positions[i] - 1;

	this->hashcode = gethash();
}

int Variant::indexOf(int exonid)
{
	return this->idmap[exonid];
}

bool Variant::contains(Exon* e)
{
	return this->idmap.count(e->id) > 0;
}

bool Variant::contains(Fragment* frag)
{
	for (int l = 0; l < frag->leftc; l++)
	{
		if (idmap.count(frag->left[l]) == 0)
		{
			return false;
		}
	}
	for (int r = 0; r < frag->rightc; r++)
	{
		if (idmap.count(frag->right[r]) == 0)
		{
			return false;
		}
	}
	return true;
}

int Variant::compare(const Variant* other)
{
	/*if (this->strand && !other->strand)
	{
		return -1;
	}
	if (!this->strand && other->strand)
	{
		return +1;
	}*/

	if (this->exonCount < other->exonCount) {
	  return -1;
	} else if (this->exonCount > other->exonCount) {
	  return +1;
	}

	for (int c = 0; c < this->codelen; c++) {
	  if (this->codes[c] < other->codes[c]) {
            return -1;
	  } else if (this->codes[c] > other->codes[c]) {
	    return +1;
	  }
	}

	return 0;
}
int Variant::gethash()
{
	int h = 0;

	for (int c = 0; c < codelen; c++)
	{
		h ^= codes[c];
	}
	h += gene->id;

	return h;
}
