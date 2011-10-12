#include <map>
#include "exon.h"

using namespace std;

class Variant
{
	friend class Casper;

public:
	int id;
	int exonCount;
	int* exons;
	int* positions;
	int length;

	Variant(int id, int exonCount);

	void addExon(Exon* e);
	void addExon(int exonid, int exonlength);
	int indexOf(int exonid);
	bool contains(Fragment* frag);
	
private:
	map<int, int> idmap;
	int num;
	int curexon;
};
