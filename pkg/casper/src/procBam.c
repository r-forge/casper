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

SEXP procBam(SEXP bamFileName, SEXP samtools, SEXP chromo){
    //void procBam(){
	FILE *bamFile;
	char *pch, *holder, tab[]="\t", comm[1000];
	read_t *frags;
	int totF, j, l, hashSize, i, left;
	hash_t *fragsHashPtr, fragsHash;
	hash_node_t *bucket;
	verbose=0;
    
	fragsHashPtr = &fragsHash;		
	hashSize=pow(2,10);
	hash_init(fragsHashPtr, hashSize);
	holder=malloc(2000*sizeof(char));
    
	
	// Alloc memory for fragments
	frags = malloc(10000000 * sizeof(read_t));
	
	char tmpqname[50];
    
	// Open bamFile
	PROTECT(bamFileName = coerceVector(bamFileName, STRSXP));	
	PROTECT(samtools = coerceVector(samtools, STRSXP));	
	PROTECT(chromo = coerceVector(chromo, STRSXP));
	if(strlen(CHAR(STRING_ELT(chromo,0)))>0) sprintf(comm, "%s/samtools view -f 2 %s %s ", CHAR(STRING_ELT(samtools,0)), CHAR(STRING_ELT(bamFileName,0)), CHAR(STRING_ELT(chromo,0)));
	else sprintf(comm, "%s/samtools view -f 2 %s ", CHAR(STRING_ELT(samtools,0)), CHAR(STRING_ELT(bamFileName,0)));
	printf("Reading file with: %s\n", comm);
	bamFile=popen(comm, "r");
    
	if (bamFile == NULL) perror("Error opening bamFile file");  

	// Read BAM file and find paired reads  
	totF=0;
	while(fgets (holder, 1000, bamFile)!=0){
	  if(verbose) printf("Read %d\n", totF);
	  j=(int)strlen(holder);
	  holder[j]='\0';
	  holder[j-1]='\0';
	  pch=strtok(m_strdup(holder), tab);
	  strcpy(tmpqname, pch);
	  l=hash_lookup(fragsHashPtr, tmpqname);
	  if(l!=HASH_FAIL) {
            addRead2Frag(holder, l, frags, 2);
            if(verbose) printf("hold: %s st: %d fl: %d cig: %s\n", holder, frags[l].st_2, frags[l].flag_2, frags[l].cigar_2);
        }
        else {
            hash_insert(fragsHashPtr, tmpqname, totF);
            addRead2Frag(holder, totF, frags, 1);
            if(verbose) printf("hold: %s st: %d fl: %d cig: %s\n", holder, frags[totF].st_1, frags[totF].flag_1, frags[totF].cigar_1);
        }
        totF++;
    }
    
	SEXP len, strs, flag, key, chrom, rid;	
	int tmp, *cigs, counter=0, *p_len, *p_strs, *p_flag, *p_rid;
	
	PROTECT(len = allocVector(INTSXP, totF*4));
	PROTECT(strs = allocVector(INTSXP, totF*4));
	PROTECT(flag = allocVector(INTSXP, totF*4));
	PROTECT(rid = allocVector(INTSXP, totF*4));
	PROTECT(key = allocVector(STRSXP, totF*4));
	PROTECT(chrom = allocVector(STRSXP, totF*4));
    
	p_len=INTEGER(len);
	p_strs=INTEGER(strs);
	p_flag=INTEGER(flag);
	p_rid=INTEGER(rid);
    
	SEXP len1, len2, st1, st2, flag1, flag2, chrom1, chrom2;
	int *p_len1, *p_len2, *p_st1, *p_st2, *p_flag1, *p_flag2;
	int count2=0;
    
	PROTECT(len1 = allocVector(INTSXP, totF));
	PROTECT(len2 = allocVector(INTSXP, totF));
	PROTECT(st1 = allocVector(INTSXP, totF));
	PROTECT(st2 = allocVector(INTSXP, totF));
	PROTECT(flag1 = allocVector(INTSXP, totF));
	PROTECT(flag2 = allocVector(INTSXP, totF));
	PROTECT(chrom1 = allocVector(STRSXP, totF));
	PROTECT(chrom2 = allocVector(STRSXP, totF));
    
	p_len1=INTEGER(len1);
	p_len2=INTEGER(len2);
	p_st1=INTEGER(st1);
	p_st2=INTEGER(st2);
	p_flag1=INTEGER(flag1);
	p_flag2=INTEGER(flag2);
    
    
	for(i=0; i<totF; i++) {
	  if(verbose) printf("%d %d\n", i, fragsHash.size);
	  if(fragsHash.bucket[i]!=NULL)  {
            bucket=fragsHash.bucket[i];
            while(bucket) {
	      tmp=bucket->data;
	      if(frags[tmp].nreads==2) {
		left=2;
		if(frags[tmp].st_1<frags[tmp].st_2) left=1;
		cigs=procCigar(m_strdup(frags[tmp].cigar_1));
		frags[tmp].len_1=frags[tmp].st_1;
		for(j=1; j<cigs[0]+1;j=j+2) {
		  SET_STRING_ELT(key, counter, mkChar(bucket->key));
		  SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_1));
		  p_len[counter] = frags[tmp].len_1;
		  p_strs[counter] = frags[tmp].len_1+cigs[j];
		  p_flag[counter] = frags[tmp].flag_1;
		  p_rid[counter] = 1;
		  frags[tmp].len_1+=cigs[j];
		  if(cigs[0]>1) frags[tmp].len_1+=cigs[j+1];
		  counter++;
		}
		cigs=procCigar(m_strdup(frags[tmp].cigar_2));
		frags[tmp].len_2=frags[tmp].st_2;
		for(j=1; j<cigs[0]+1;j=j+2) {        
		  SET_STRING_ELT(key, counter, mkChar(bucket->key));
		  SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_2));
		  p_len[counter] = frags[tmp].len_2;
		  p_strs[counter] = frags[tmp].len_2+cigs[j];
		  p_flag[counter] = frags[tmp].flag_2;
		  p_rid[counter] = 2;
		  frags[tmp].len_2+=cigs[j];
		  if(cigs[0]>1) frags[tmp].len_2+=cigs[j+1];
		  counter++;
		}
		if(verbose) printf("%s %d %d %s\n", bucket->key, frags[tmp].st_1, frags[tmp].flag_1, frags[tmp].cigar_1);
		if(verbose) printf("%s %d %d %s %d %d\n", bucket->key, frags[tmp].st_2, frags[tmp].flag_2, frags[tmp].cigar_2, cigs[0], cigs[1]);
		SET_STRING_ELT(chrom1, count2, mkChar(frags[tmp].chr_1));
		SET_STRING_ELT(chrom2, count2, mkChar(frags[tmp].chr_2));
		p_st1[count2] = frags[tmp].st_1;
		p_st2[count2] = frags[tmp].st_2;
		p_len1[count2] = frags[tmp].len_1;
		p_len2[count2] = frags[tmp].len_2;
		p_flag1[count2] = frags[tmp].flag_1;
		p_flag2[count2] = frags[tmp].flag_2;
		// if((int)p_len2[count2]*(int)p_len1[count2]<0) printf("%s %d %d %d %d %d %d %s %d\n", bucket->key, p_st1[count2], p_st2[count2], p_len1[count2],p_len2[count2], p_flag1[count2], p_flag2[count2], frags[tmp].cigar_1, (int)p_len2[count2]*(int)p_len1[count2]);
		count2++;
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
	free(holder);
	UNPROTECT(22);
	return(ans);
}
