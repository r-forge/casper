#include <R.h>
#include <Rinternals.h>
#include "stdafx.h"
#include "casper.h"
#include "rcasper.h"


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
  UNPROTECT(2);

  return res;
}

//double cumu_fragstaLIN(double x) { return x; }

extern "C"
{

  //Compute expression for known variants.
  SEXP calc(SEXP exons, SEXP transcripts, SEXP pathCounts, SEXP fragsta, SEXP fraglen, SEXP lenvals, SEXP readLength) 
	{
	        ::fun_fragsta = fragsta;
		
		// DiscreteDF
		SEXP lims;
		PROTECT(lims = getAttrib(fraglen, R_DimSymbol));
		int ln = INTEGER(lims)[0];
		double* ld = REAL(fraglen);
		UNPROTECT(1);
		int* lv= INTEGER(lenvals);

		DiscreteDF* lenfun = new DiscreteDF(ld, lv, ln);

		// Main
		DataFrame* df = new DataFrame(lenfun, cumu_fragsta);
                df->frag_readlen= INTEGER(readLength)[0];
		
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

			df->addData(f);
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

				Exon* ex = df->exons[eid];
				el->push_back(ex);
			}

			// TODO: read strand, positive or negative?
			bool strand = true;

			// TODO: read gene id of variant here
			int gid = -1;
			Gene* gene = df->genes[gid];

			//Variant* v = new Variant(gene, strand, el);
			Variant* v = new Variant(gene, el);
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


  SEXP calcDenovo(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorprobR) {
  //Perform de novo isoform discovery and estimate expression for a single gene
  //Input
  // - exons: vector with exon ids
  // - exonwidth: vector exon widths
  // - geneid: integer with gene id
  // - transcripts: list of transcripts. Names indicate transcript id. Each element contains list of exon ids. Used to initialize model
  // - pathCounts: vector with path counts. Names indicate series of visited exons e.g. e1.e2-e5
  // - fragsta: function that returns start distrib cdf
  // - fraglen: vector with fragment length distrib, i.e. P(length=lenvals[0]),P(length=lenvals[1]),... up to max length
  // - lenvals: vector with possibles values for length
  // - readLength: read length
  // - priorprob: unfinished. This should somehow allow to compute prior prob on model space

    int *exons= INTEGER(exonsR), *exonwidth= INTEGER(exonwidthR), geneid= INTEGER(geneidR)[0], *lenvals= INTEGER(lenvalsR);
    double *fraglen= REAL(fraglenR); 
    SEXP nexonsR = getAttrib(exonsR,R_DimSymbol), nfraglenR= getAttrib(fraglenR,R_DimSymbol);
    int nexons= INTEGER(nexonsR)[0], nfraglen= INTEGER(nfraglenR)[0], readLength= INTEGER(readLengthR)[0];

    Casper* casp = initCasper(exons, exonwidth, transcriptsR, geneid, nexons, pathCountsR, fraglen, lenvals, nfraglen, readLength, fragstaR);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0]= 1.0;
    UNPROTECT(1);
    return(ans);
  }

}


Casper* initCasper(int *exons, int *exonwidth, SEXP transcriptsR, int geneid, int nexons, SEXP pathCountsR, double *fraglen, int *lenvals, int nfraglen, int readLength, SEXP fragstaR) {

  //Define fragment length/start distributions and create DataFrame
  DiscreteDF* fraglen_dist = new DiscreteDF(fraglen, lenvals, nfraglen);

  ::fun_fragsta = fragstaR;  //Manuel: cumu_fragsta does not seem to use its argument x. Is this correct?
  DataFrame* df = new DataFrame(fraglen_dist, cumu_fragsta);
  df->frag_readlen= readLength;

  //Add exons & gene to DataFrame
  Gene* g= new Gene(geneid);
  for (int i=0; i<nexons; i++) {
    df->addExon(new Exon(exons[i], exonwidth[i]));
    g->addExon(new Exon(exons[i], exonwidth[i]));  //Manuel: do we really need to store exon info in the gene? It seems that df has this info twice
  }
  df->addGene(g);

  //Add exon path counts
  int np = LENGTH(pathCountsR);
  SEXP pnames = getAttrib(pathCountsR, R_NamesSymbol);
  int* pathCounts = INTEGER(pathCountsR);  

  for (int i=0; i<np; i++) {
    //Set counts
    int count = pathCounts[i];

    //Set nb of left & right visited exons
    const char* pname = CHAR(STRING_ELT(pnames, i));
    char* left = new char[strlen(pname) + 1];
    strcpy(left, pname);
    char* mid = strchr(left, '-');
    if (mid == NULL) continue;
    mid[0] = '\0';
    char* right = mid+1;
    int leftc = 0, rightc= 0;
    for (int l = strlen(left)-1; l >= 0; l--) { if (left[l] == '.') leftc++; }  
    for (int r = strlen(right)-1; r >= 0; r--) { if (right[r] == '.') rightc++; } 

    Fragment* f = new Fragment(leftc, rightc, count);

    //Set sequence of visited exons
    left = left+1;
    right[strlen(right)-1] = '\0';
    char* item;
    item = strtok(left, ".");  //split string
    for (int i = 0; item != NULL; i++) {
      int eid = atoi(item);
      f->left[i] = eid;
      item = strtok(NULL, ".");
    }
    item = strtok(right, ".");
    for (int i = 0; item != NULL; i++) {
      int eid = atoi(item);
      f->right[i] = eid;
      item = strtok(NULL, ".");
    }

    df->addData(f);
  }

  //Set initial variants
  int nt = LENGTH(transcriptsR);
  vector<Variant*>* initvars = new vector<Variant*>();
  SEXP tnames = getAttrib(transcriptsR, R_NamesSymbol);
  Gene* gene;

  for (int i = 0; i < nt; i++) {
    int tid = atoi(CHAR(STRING_ELT(tnames, i)));
    SEXP trow = VECTOR_ELT(transcriptsR, i);
    int ntsub = LENGTH(trow);
    int* tvals = INTEGER(trow);

    vector<Exon*>* el = new vector<Exon*>();
    for (int s = 0; s < ntsub; s++) {
      int eid = tvals[s];
      Exon* ex = df->exons[eid];
      el->push_back(ex);
    }

    gene = df->genes[geneid];
    Variant* v = new Variant(gene, el);
    v->id = tid;

    initvars->push_back(v);
  }

  //Add variants to initial set so that all observed exon paths have positive probability
  int deletedPath= 0;
  list<Fragment*>::iterator it;
  vector<Variant*>::iterator it2;
  for (it=df->data.begin(); it != df->data.end(); it++) {
    int posprob=0, newct=0;
    if (deletedPath==1) { it--; deletedPath=0; }
    for (it2= initvars->begin(); (it2 != initvars->end()) && posprob==0; it2++) {
      if (df->probability(*it2, *it) > 0) posprob=1;  //probability(Variant *v, Fragment *it)
    }
    //Add new variant containing observed exon path (plus all exons before & after the path)
    if (posprob==0) {  
      int eid; Exon* ex;
      vector<Exon*>* el = new vector<Exon*>();
      for (int i=0; i< (*it)->left[0]; i++) {
        ex= df->exons[i];
        el->push_back(ex);
      }
      for (int i=0; i< (*it)->leftc; i++) {
        eid = (*it)->left[i];
        ex = df->exons[eid];
        el->push_back(ex);
      }
      if (eid != (*it)->right[0]) {
        eid = (*it)->right[0];
        ex = df->exons[eid];
        el->push_back(ex);
      }
      for (int i=1; i< (*it)->rightc; i++) {
        eid = (*it)->right[i];
        ex = df->exons[eid];
        el->push_back(ex);
      }
      int nexon= gene->exons.size();
      for (int i=eid; i<nexon; i++) {
        ex= df->exons[i];
        el->push_back(ex);
      }
      Variant* v = new Variant(gene, el);
      v->id = nt+newct;
      newct++;
      initvars->push_back(v);
      //If path count still has 0 prob, remove path count
      if (df->probability(v, *it)==0) { 
        it= df->data.erase(it);
        deletedPath= 1;
        Rprintf("Observed path count could not be explained, so it was removed\n"); 
      }  
    }
  }

  Model* model = new Model(initvars);
  Casper* casp = new Casper(model, df);

  return casp;
}

