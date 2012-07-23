#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "header.h"
#include "hash.h"
#include "functions.h"


SEXP procBam(SEXP qname, SEXP flags, SEXP chr, SEXP start, SEXP cigar, SEXP strand, SEXP totFrags, SEXP totReads, SEXP len, SEXP strs, SEXP flag, SEXP key, SEXP chrom, SEXP rid, SEXP rstrand){
	read_t *frags;
	int totF, j, l, hashSize, i, *frags_size, *reads_size;
	hash_t *fragsHashPtr, fragsHash;
	hash_node_t *bucket;
	verbose=0;
    
	
	// Alloc memory for fragments
	PROTECT(totFrags = coerceVector(totFrags, INTSXP));
	frags_size = INTEGER(totFrags);
	frags = malloc((frags_size[0] + 1) * sizeof(read_t));
	
    	fragsHashPtr = &fragsHash;		
	hashSize=frags_size[0];
	hash_init(fragsHashPtr, hashSize);
    
	PROTECT(qname = coerceVector(qname, STRSXP));
	PROTECT(flags = coerceVector(flags, INTSXP));
	int *p_flags=INTEGER(flags);
	PROTECT(chr = coerceVector(chr, STRSXP));
	PROTECT(start = coerceVector(start, INTSXP));
	int *p_start=INTEGER(start);
	PROTECT(cigar = coerceVector(cigar, STRSXP));
	PROTECT(totReads = coerceVector(totReads, INTSXP));
	PROTECT(strand = coerceVector(strand, INTSXP));
	int *p_strand=INTEGER(strand);
	reads_size = INTEGER(totReads);
	
	totF=0;
		// Find paired reads and fill fragment's object 
	for (i=0; i<frags_size[0]; i++) {
   	  l=hash_lookup(fragsHashPtr, CHAR(STRING_ELT(qname, i)));
	  if(l!=HASH_FAIL) {
            addRead2Frag(CHAR(STRING_ELT(qname, i)), p_flags[i], CHAR(STRING_ELT(chr, i)), p_start[i], p_strand[i], CHAR(STRING_ELT(cigar, i)), l, frags, 2);
	    //if(verbose) printf("qname: %s orig_cig %s st: %d fl: %d cig: %s l: %d\n", qname[i], cigar[i], frags[l].st_2, frags[l].flag_2, frags[l].cigar_2, l);
	  }
	  else {
            hash_insert(fragsHashPtr, CHAR(STRING_ELT(qname, i)), totF);
            addRead2Frag(CHAR(STRING_ELT(qname, i)), p_flags[i], CHAR(STRING_ELT(chr, i)), p_start[i], p_strand[i], CHAR(STRING_ELT(cigar, i)), totF, frags, 1);
//   if(verbose) printf("qname: %s orig_cig: %s st: %d fl: %d cig: %s\n", qname[i], cigar[i], frags[totF].st_1, frags[totF].flag_1, frags[totF].cigar_1);
            totF++;
	  }
	}
	if(verbose) printf("%d %s %d %d %d\n", totF, frags[0].cigar_1, frags[0].st_1, frags[0].flag_1, frags[0].nreads);
	
	//SEXP len, strs, flag, key, chrom, rid;	
	int tmp, *cigs, counter=0;
	
	PROTECT(len = coerceVector(len, INTSXP));
	int *p_len=INTEGER(len);
	PROTECT(flag = coerceVector(flag, INTSXP));
	int *p_flag=INTEGER(flag);
	PROTECT(chrom = coerceVector(chrom, STRSXP));
	PROTECT(strs = coerceVector(strs, INTSXP));
	int *p_strs=INTEGER(strs);
	PROTECT(key = coerceVector(key, STRSXP));
	PROTECT(rid = coerceVector(rid, INTSXP));
	int *p_rid=INTEGER(rid);
        PROTECT(rstrand = coerceVector(rstrand, INTSXP));
	int *p_rstrand = INTEGER(rstrand);
	for(i=0; i<fragsHash.size; i++) {
	  if(verbose) printf("%d %d\n", i, fragsHash.size);
	  if(fragsHash.bucket[i]!=NULL)  {
            bucket=fragsHash.bucket[i];
            while(bucket) {
	      tmp=bucket->data;
	      cigs=malloc(50 * sizeof(int));
	      if(verbose) printf("nreads %d\n", frags[tmp].nreads);
	      if(frags[tmp].nreads==2) {
		cigs=procCigar(m_strdup(frags[tmp].cigar_1), cigs);
		frags[tmp].len_1=frags[tmp].st_1;
		if(verbose) printf("proc cigs %s %d", frags[tmp].cigar_1, cigs[0]);
		for(j=1; j<cigs[0]+1;j=j+2) {
		  SET_STRING_ELT(key, counter, mkChar(bucket->key));
		  SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_1));
		  p_strs[counter] = frags[tmp].len_1;
		  p_len[counter] = frags[tmp].len_1+(cigs[j]-1);
		  p_flag[counter] = frags[tmp].flag_1;
		  p_rid[counter] = 1;
		  p_rstrand[counter] = frags[tmp].strand_1;
		  frags[tmp].len_1+=cigs[j];
		  if((cigs[0]>1)&&(j<cigs[0]-1)) frags[tmp].len_1+=cigs[j+1];
		  if(verbose) printf("1 %d %d %d %d %d\n", p_strs[counter], p_flag[counter], p_len[counter], counter, p_strand[counter]);
		  counter++;
		}
		cigs=procCigar(m_strdup(frags[tmp].cigar_2), cigs);
		frags[tmp].len_2=frags[tmp].st_2;
		for(j=1; j<cigs[0]+1;j=j+2) {        
		  SET_STRING_ELT(key, counter, mkChar(bucket->key));
		  SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_2));
		  p_strs[counter] = frags[tmp].len_2;
		  p_len[counter] = frags[tmp].len_2+(cigs[j]-1);
		  p_flag[counter] = frags[tmp].flag_2;
		  p_rid[counter] = 2;
		  p_rstrand[counter] = frags[tmp].strand_2;
		  frags[tmp].len_2+=cigs[j];
		  if((cigs[0]>1)&&(j<cigs[0]-1)) frags[tmp].len_2+=cigs[j+1];
		  if(verbose) printf("2 %d %d %d %d\n", p_strs[counter], p_flag[counter], p_len[counter], counter);
		  counter++;
                  
		}
		free(frags[tmp].chr_2);
		free(frags[tmp].cigar_2);      
	      }
	      if(verbose) printf("%s %d %d %s\n", bucket->key, frags[tmp].st_1, frags[tmp].flag_1, frags[tmp].cigar_1);
	      free(frags[tmp].qname);
	      free(frags[tmp].chr_1);
	      free(frags[tmp].cigar_1);
	      free(cigs);
	      bucket=bucket->next;
            }
	  }
	}
    
	SEXP reads;
	PROTECT(reads = allocVector(VECSXP, 7));
	SET_VECTOR_ELT(reads, 0, len);
	SET_VECTOR_ELT(reads, 1, strs);
	SET_VECTOR_ELT(reads, 2, flag);
	SET_VECTOR_ELT(reads, 3, key);
	SET_VECTOR_ELT(reads, 4, chrom);
	SET_VECTOR_ELT(reads, 5, rid);
	SET_VECTOR_ELT(reads, 6, rstrand);

  	free(frags);
	hash_destroy(fragsHashPtr);
	UNPROTECT(16);
	return(reads);
}
