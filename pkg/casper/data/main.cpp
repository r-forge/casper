#include <R.h>
#include <Rinternals.h>

extern "C"
{

SEXP calc(SEXP exons, SEXP transcripts, SEXP pathCounts) 
{
	// Exons
	SEXP dims;
	PROTECT(dims = getAttrib(exons, R_DimSymbol));
	int ne = INTEGER(dims)[0];
	UNPROTECT(1);
	
	int* me = INTEGER(exons);

	REprintf("Reading %i exons..\n", ne);
	for (int i = 0; i < 10; i++)
	{
		int id = me[0*ne+i];
		int start = me[1*ne+i];
		int end = me[2*ne+i];
	}

	// Variants
	int nt = LENGTH(transcripts);

	SEXP tnames = getAttrib(transcripts, R_NamesSymbol);

	REprintf("Reading %i transcripts...\n", nt);
	for (int i = 0; i < 10; i++)
	{
		int tid = atoi(CHAR(STRING_ELT(tnames, i)));

		SEXP trow = VECTOR_ELT(transcripts, i);
		int ntsub = LENGTH(trow);
		int* tvals = INTEGER(trow);

		REprintf("AAAAAA %i:\n", tid);
		for (int c = 0; c < ntsub; c++)
		{
			int eid = tvals[c];
			REprintf("%i\n", eid);
		}
	}

	// PathCounts
	int np = LENGTH(pathCounts);

	SEXP pnames = getAttrib(pathCounts, R_NamesSymbol);

	REprintf("Reading %i pathCounts...\n", np);
	for (int i = 0; i < 10; i++)
	{
		const char* pname = CHAR(STRING_ELT(pnames, i));
		int count = INTEGER(VECTOR_ELT(transcripts, i))[0];
		
		char* left = new char[strlen(pname)];
		strcpy(left, pname);
		char* mid = strchr(left, '-');
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

		REprintf("%i %s xxx %i %s xxx %i\n", leftc, left, rightc, right, count);

		char* item;
		item = strtok(left, ".");
		while (item != NULL)
		{
			int eid = atoi(item);
			REprintf("L %i\n", eid);
			item = strtok(NULL, ".");
		}

		item = strtok(right, ".");
		while (item != NULL)
		{
			int eid = atoi(item);
			REprintf("R %i\n", eid);
			item = strtok(NULL, ".");
		}
	}
	
	// Output
	int vc = 3;

	SEXP Rc;
	PROTECT(Rc = allocVector(REALSXP, vc));
	double* res = REAL(Rc);
	res[0] = 0.7;
	res[1] = 0.2;
	res[2] = 0.1;
	UNPROTECT(1);

	return(Rc);
}

}
