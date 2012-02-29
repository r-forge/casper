#include <R.h>
#include <Rinternals.h>
#include "stdafx.h"
#include "casper.h"

extern "C"
{
	SEXP fun_fragsta;

	double cumu_fragsta(double x)
	{
		SEXP sval;
		PROTECT(sval = allocVector(REALSXP, 1));
		//double* val = REAL(sval);

		SEXP R_fcall;
		PROTECT(R_fcall = lang2(fun_fragsta, R_NilValue));
		SETCADR(R_fcall, sval);

		SEXP funval = eval(R_fcall, R_NilValue);
		double res = REAL(funval)[0];
		UNPROTECT(1);

		return res;
	}
	double cumu_fragstaLIN(double x)
	{
		return x;
	}

	SEXP calc(SEXP exons, SEXP transcripts, SEXP pathCounts, SEXP fragsta, SEXP fraglen) 
	{
		fun_fragsta = fragsta;
		
		// DiscreteDF
		SEXP lims;
		PROTECT(lims = getAttrib(fraglen, R_DimSymbol));
		int ln = INTEGER(lims)[0];
		double* ld = REAL(fraglen);
		UNPROTECT(1);

		DiscreteDF* lenfun = new DiscreteDF(ld, ln);

		// Main
		DataFrame* df = new DataFrame(lenfun, cumu_fragsta);
		
		// Exons
		SEXP dims;
		PROTECT(dims = getAttrib(exons, R_DimSymbol));
		int ne = INTEGER(dims)[0];
		UNPROTECT(1);

		int* me = INTEGER(exons);

		REprintf("Reading %i exons..\n", ne);
		for (int i = 0; i < ne; i++)
		{
			int id = me[0*ne+i];
			int length = me[1*ne+i];

			df->addExon(new Exon(id, length));
		}

		// PathCounts
		int np = LENGTH(pathCounts);

		SEXP pnames = getAttrib(pathCounts, R_NamesSymbol);
		int* pvalus = INTEGER(pathCounts);

		REprintf("Reading %i pathCounts...\n", np);
		for (int i = 0; i < np; i++)
		{
			const char* pname = CHAR(STRING_ELT(pnames, i));
			int count = pvalus[i];

			char* left = new char[strlen(pname) + 1];
			strcpy(left, pname);
			char* mid = strchr(left, '-');
			if (mid == NULL)
			{
				continue;
			}
			mid[0] = '\0';
			char* right = mid+1;

			int leftc = 0;
			for (int l = strlen(left)-1; l >= 0; l--)
			{
				if (left[l] == '.')
				{
					leftc++;
				}
			}
			int rightc = 0;
			for (int r = strlen(right)-1; r >= 0; r--)
			{
				if (right[r] == '.')
				{
					rightc++;
				}
			}

			left = left+1;
			right[strlen(right)-1] = '\0';

			Fragment* f = new Fragment(leftc, rightc, count);

			char* item;
			item = strtok(left, ".");
			for (int i = 0; item != NULL; i++)
			{
				int eid = atoi(item);
				f->left[i] = eid;

				item = strtok(NULL, ".");
			}

			item = strtok(right, ".");
			for (int i = 0; item != NULL; i++)
			{
				int eid = atoi(item);
				f->right[i] = eid;

				item = strtok(NULL, ".");
			}

			df->addFragment(f);
		}

		// Genes
		// TODO: add genes here like this of positive strand
		Gene* g = new Gene(1);
		g->addExon(df->exons[1]);
		g->addExon(df->exons[2]);
		g->addExon(df->exons[3]);
		df->addGene(g);

		// Variants
		int nt = LENGTH(transcripts);
		vector<Variant*>* vm = new vector<Variant*>();
		SEXP tnames = getAttrib(transcripts, R_NamesSymbol);

		REprintf("Reading %i transcripts...\n", nt);
		for (int i = 0; i < nt; i++)
		{
			int tid = atoi(CHAR(STRING_ELT(tnames, i)));

			SEXP trow = VECTOR_ELT(transcripts, i);
			int ntsub = LENGTH(trow);
			int* tvals = INTEGER(trow);

			vector<Exon*>* el = new vector<Exon*>();

			for (int s = 0; s < ntsub; s++)
			{
				int eid = tvals[s];

				Exon* ex = c->getExon(eid);
				el->push_back(ex);
			}

			// TODO: read strand, positive or negative?
			bool strand = true;

			// TODO: read gene id of variant here
			int gid = -1;
			Gene* gene = df->genes[gid];

			Variant* v = new Variant(gene, strand, el);
			v->id = tid;

			vm->push_back(v);
		}
		Model* mdl = new Model(vm);

		
/*
		// Debug
		map<Variant*, map<Fragment*, double> >::const_iterator mi;
		for (mi = c->memvprobs.begin(); mi != c->memvprobs.end(); mi++)
		{
			map<Fragment*, double>::const_iterator ni;
			for (ni = mi->second.begin(); ni != mi->second.end(); ni++) 
			{
				if (ni->first->leftc == 1 && ni->first->rightc == 2)
				{
				//	REprintf("%i\t%i %i %i\t%f\n", mi->first->id, ni->first->left[0], ni->first->right[0], ni->first->right[1], ni->second);
				}
			}
		}
//		REprintf("\n");
		
		map<int, Exon*>::const_iterator ei;
		for (ei = c->exons.begin(); ei != c->exons.end(); ei++)
		{
			Exon* e = ei->second;
//			REprintf("%i\t%i\n", e->id, e->length);
		}
//		REprintf("\n");
		
		list<Fragment*>::const_iterator fi;
		for (fi = c->fragments.begin(); fi != c->fragments.end(); fi++)
		{
			Fragment* f = *fi;
		//	REprintf("%i\t%i\t%i\n", f->leftc, f->rightc, f->count);
			for (int l = 0; l < f->leftc; l++)
			{
		//		REprintf("%i\n", f->left[l]);
			}
			for (int r = 0; r < f->rightc; r++)
			{
		//		REprintf("%i\n", f->right[r]);
			}
		}
	//	REprintf("\n");

		vector<Variant*>::const_iterator vi;
		for (vi = c->variants.begin(); vi != c->variants.end(); vi++)
		{
			Variant* v = *vi;
	//		REprintf("%i\t%i\n", v->id, v->exonCount);
			for (int e = 0; e < v->exonCount; e++)
			{
	//			REprintf("%i\n", v->exons[e]);
			}
		}
	//	REprintf("\n");
*/

		// Output
		Casper* c = new Casper(mdl, df);
		int vc = c->model->count();
		double* em = c->calculateMode();

		SEXP Rc;
		PROTECT(Rc = allocVector(REALSXP, vc));
		double* res = REAL(Rc);
		for (int i = 0; i < vc; i++)
		{
			res[i] = em[i];
		}
		UNPROTECT(1);

		return(Rc);
	}
}
