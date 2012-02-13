#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "header.h"
#include "hash.h"
#include "functions.h"
#include "fragFunc.h"

SEXP procBam(SEXP qname, SEXP flags, SEXP chr, SEXP start, SEXP cigar, SEXP totFrags, SEXP totReads){
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
	reads_size = INTEGER(totReads);
	
    totF=0;
    	// Find paired reads and fill fragment's object 
	for (i=0; i<frags_size[0]; i++) {
   	  l=hash_lookup(fragsHashPtr, CHAR(STRING_ELT(qname, i)));
	  if(l!=HASH_FAIL) {
            addRead2Frag(CHAR(STRING_ELT(qname, i)), p_flags[i], CHAR(STRING_ELT(chr, i)), p_start[i], CHAR(STRING_ELT(cigar, i)), l, frags, 2);
// if(verbose) printf("qname: %s orig_cig %s st: %d fl: %d cig: %s l: %d\n", qname[i], cigar[i], frags[l].st_2, frags[l].flag_2, frags[l].cigar_2, l);
		}
        else {
            hash_insert(fragsHashPtr, CHAR(STRING_ELT(qname, i)), totF);
            addRead2Frag(CHAR(STRING_ELT(qname, i)), p_flags[i], CHAR(STRING_ELT(chr, i)), p_start[i], CHAR(STRING_ELT(cigar, i)), totF, frags, 1);
//   if(verbose) printf("qname: %s orig_cig: %s st: %d fl: %d cig: %s\n", qname[i], cigar[i], frags[totF].st_1, frags[totF].flag_1, frags[totF].cigar_1);
            totF++;
        }
    }
	if(verbose) printf("%d %s %d %d %d\n", totF, frags[0].cigar_1, frags[0].st_1, frags[0].flag_1, frags[0].nreads);
    
	SEXP len, strs, flag, key, chrom, rid;	
	int tmp, *cigs, counter=0, *p_len, *p_strs, *p_flag, *p_rid;
	
	PROTECT(len = allocVector(INTSXP, reads_size[0]));
	PROTECT(strs = allocVector(INTSXP, reads_size[0]));
	PROTECT(flag = allocVector(INTSXP, reads_size[0]));
	PROTECT(rid = allocVector(INTSXP, reads_size[0]));
	PROTECT(key = allocVector(STRSXP, reads_size[0]));
	PROTECT(chrom = allocVector(STRSXP, reads_size[0]));
    
	p_len=INTEGER(len);
	p_strs=INTEGER(strs);
	p_flag=INTEGER(flag);
	p_rid=INTEGER(rid);
    
	SEXP len1, len2, st1, st2, flag1, flag2, chrom1, chrom2;
	int *p_len1, *p_len2, *p_st1, *p_st2, *p_flag1, *p_flag2;
	int count2=0;
    
	PROTECT(len1 = allocVector(INTSXP, frags_size[0]+1));
	PROTECT(len2 = allocVector(INTSXP, frags_size[0]+1));
	PROTECT(st1 = allocVector(INTSXP, frags_size[0]+1));
	PROTECT(st2 = allocVector(INTSXP, frags_size[0]+1));
	PROTECT(flag1 = allocVector(INTSXP, frags_size[0]+1));
	PROTECT(flag2 = allocVector(INTSXP, frags_size[0]+1));
	PROTECT(chrom1 = allocVector(STRSXP, frags_size[0]+1));
	PROTECT(chrom2 = allocVector(STRSXP, frags_size[0]+1));
    
	p_len1=INTEGER(len1);
	p_len2=INTEGER(len2);
	p_st1=INTEGER(st1);
	p_st2=INTEGER(st2);
	p_flag1=INTEGER(flag1);
	p_flag2=INTEGER(flag2);
    
	
    for(i=0; i<fragsHash.size; i++) {
	  if(verbose) printf("%d %d\n", i, fragsHash.size);
	  if(fragsHash.bucket[i]!=NULL)  {
            bucket=fragsHash.bucket[i];
            while(bucket) {
	      tmp=bucket->data;
	      if(frags[tmp].nreads==2) {
		cigs=procCigar(m_strdup(frags[tmp].cigar_1));
		frags[tmp].len_1=frags[tmp].st_1;
		for(j=1; j<cigs[0]+1;j=j+2) {
		  SET_STRING_ELT(key, counter, mkChar(bucket->key));
		  SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_1));
		  p_strs[counter] = frags[tmp].len_1;
		  p_len[counter] = frags[tmp].len_1+(cigs[j]-1);
		  p_flag[counter] = frags[tmp].flag_1;
		  p_rid[counter] = 1;
		  frags[tmp].len_1+=cigs[j];
		  if((cigs[0]>1)&&(j<cigs[0]-1)) frags[tmp].len_1+=cigs[j+1];
		  if(verbose) printf("1 %d %d %d %d\n", p_strs[counter], p_flag[counter], p_len[counter], counter);
		  counter++;
		}
		cigs=procCigar(m_strdup(frags[tmp].cigar_2));
		frags[tmp].len_2=frags[tmp].st_2;
		for(j=1; j<cigs[0]+1;j=j+2) {        
		  SET_STRING_ELT(key, counter, mkChar(bucket->key));
		  SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_2));
		  p_strs[counter] = frags[tmp].len_2;
		  p_len[counter] = frags[tmp].len_2+(cigs[j]-1);
		  p_flag[counter] = frags[tmp].flag_2;
		  p_rid[counter] = 2;
		  frags[tmp].len_2+=cigs[j];
		  if((cigs[0]>1)&&(j<cigs[0]-1)) frags[tmp].len_2+=cigs[j+1];
	if(verbose) printf("2 %d %d %d %d\n", p_strs[counter], p_flag[counter], p_len[counter], counter);
            counter++;
		}
		if(verbose) printf("%s %d %d %s count %d \n", bucket->key, frags[tmp].st_1, frags[tmp].flag_1, frags[tmp].cigar_1, counter);
		if(verbose) printf("%s %d %d %s %d %d count %d \n", bucket->key, frags[tmp].st_2, frags[tmp].flag_2, frags[tmp].cigar_2, cigs[0], cigs[1], counter);
		if(verbose) if((frags[tmp].len_1<0)||(frags[tmp].len_2<0)) printf("caca %s %d %d %d %d %d count %d \n", bucket->key, frags[tmp].st_1, frags[tmp].len_1, frags[tmp].st_2, frags[tmp].len_2, j, counter);
		SET_STRING_ELT(chrom1, count2, mkChar(frags[tmp].chr_1));
		SET_STRING_ELT(chrom2, count2, mkChar(frags[tmp].chr_2));
		p_st1[count2] = frags[tmp].st_1;
		p_st2[count2] = frags[tmp].st_2;
		p_len1[count2] = frags[tmp].len_1;
		p_len2[count2] = frags[tmp].len_2;
		p_flag1[count2] = frags[tmp].flag_1;
		p_flag2[count2] = frags[tmp].flag_2;
		free(frags[tmp].qname);
		count2++;
              //if(count2>frags_size[0]){
                //  printf("Too many fragments\n");
                 // break;
             // }
	      }
	      if(verbose) printf("%s %d %d %s %d\n", bucket->key, frags[tmp].st_1, frags[tmp].flag_1, frags[tmp].cigar_1, count2);
	      bucket=bucket->next;
            }
	  }
	}
    
	SEXP ans, totFr, totRe, lengths, reads;	
	int *p_totFr, *p_totRe;
    
	PROTECT(ans = allocVector(VECSXP, 2));
	PROTECT(lengths = allocVector(VECSXP, 9));
	PROTECT(reads = allocVector(VECSXP, 7));
	PROTECT(totFr = allocVector(INTSXP, 1));
	PROTECT(totRe = allocVector(INTSXP, 1));
	p_totFr = INTEGER(totFr);
	p_totFr[0] = counter;
	p_totRe = INTEGER(totRe);
	p_totRe[0] = count2;
    
	SET_VECTOR_ELT(ans, 0, reads);
	SET_VECTOR_ELT(ans, 1, lengths);
	
	SET_VECTOR_ELT(reads, 0, len);
	SET_VECTOR_ELT(reads, 1, strs);
	SET_VECTOR_ELT(reads, 2, flag);
	SET_VECTOR_ELT(reads, 3, key);
	SET_VECTOR_ELT(reads, 4, chrom);
	SET_VECTOR_ELT(reads, 5, rid);
	SET_VECTOR_ELT(reads, 6, totFr);
    
	SET_VECTOR_ELT(lengths, 0, st1);
	SET_VECTOR_ELT(lengths, 1, len1);
	SET_VECTOR_ELT(lengths, 2, st2);
	SET_VECTOR_ELT(lengths, 3, len2);
	SET_VECTOR_ELT(lengths, 4, flag1);
	SET_VECTOR_ELT(lengths, 5, flag2);
	SET_VECTOR_ELT(lengths, 6, chrom1);
	SET_VECTOR_ELT(lengths, 7, chrom2);
	SET_VECTOR_ELT(lengths, 8, totRe);
    
	free(frags);
	UNPROTECT(26);
	return(ans);
}
