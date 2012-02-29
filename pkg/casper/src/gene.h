#include "exon.h"
#include <vector>
#include <map>

using namespace std;

class Gene
{
public:
	int id;
	// exons of this gene in their correct order
	vector<Exon*> exons;
	int length;

	map<Exon*, int> idmap;

	Gene(int id);
	void addExon(Exon* e);
	bool contains(Exon* e);
	int indexOf(Exon* e);
};
