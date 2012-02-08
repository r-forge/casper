#include "exon.h"
#include <vector>

#ifdef __GNUC__
#include <ext/hash_map>
#else
#include <hash_map>
#endif

using namespace std;

class Gene
{
public:
	int id;
	// exons of this gene in their correct order
	vector<Exon*> exons;
	int length;

	hash_map<Exon*, int> idmap;

	Gene(int id);
	void addExon(Exon* e);
	bool contains(Exon* e);
	int indexOf(Exon* e);
};