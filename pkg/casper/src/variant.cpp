#include "stdafx.h"
#include "fragment.h"
#include "variant.h"

Variant::Variant(int id, int exonCount)
{
	this->id = id;
	this->exonCount = exonCount;
	this->exons = new int[exonCount];
	this->positions = new int[exonCount + 1];
	this->positions[0] = 1;
	this->curexon = 0;
}

void Variant::addExon(Exon* e)
{
	Variant::addExon(e->id, e->length);
}
void Variant::addExon(int id, int length)
{
	this->idmap[id] = curexon;
	this->exons[curexon] = id;
	this->positions[curexon + 1] = this->positions[curexon] + length;

	this->length = this->positions[curexon + 1] - 1;
	this->curexon++;
}

int Variant::indexOf(int exonid)
{
	return this->idmap[exonid];
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
