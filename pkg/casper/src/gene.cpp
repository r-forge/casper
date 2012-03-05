#include "gene.h"

Gene::Gene(int id)
{
    this->id = id;
	this->length = 0;
}

void Gene::addExon(Exon* e)
{
	int i = this->exons.size();
	this->idmap[e] = i;
	this->exons.push_back(e);
	this->length += e->length;
}
bool Gene::contains(Exon* e)
{
	return (this->idmap.count(e) > 0);
}
int Gene::indexOf(Exon* e)
{
	if (this->idmap.count(e) == 0)
	{
		return -1;
	}
	return this->idmap[e];
}
