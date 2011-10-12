#include "stdafx.h"
#include "exon.h"

Exon::Exon(int id, int start, int end)
{
	this->id = id;
	this->start = start;
	this->end = end;
	this->length = end - start;
}
