#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"
#include "header.h"

//********** Functions ***********

int *procCigar(char *cigar){
	char *pch, *tab="DIMN";
	int *cigs;
	cigs=malloc(50 * sizeof(int));
	pch=strtok(cigar, tab);
	cigs[0]=0;
	while(pch!=NULL){
		sscanf(pch, "%d", &cigs[cigs[0]+1]);
		pch=strtok(NULL, tab);
		cigs[0]++;
	}
	return(cigs);
}


void addRead2Frag(const char *qname, int flag, const char *chr, int start, const char *cigar, int totF, read_t *frags, int read){
	if(read==1){
        frags[totF].qname = malloc(strlen(qname) * sizeof(char));
		strcpy(frags[totF].qname, qname);
        frags[totF].flag_1=flag;
        strcpy(frags[totF].chr_1, chr);
     	frags[totF].st_1=start;
     	strcpy(frags[totF].cigar_1, cigar);
  		frags[totF].nreads=1;
	} else {    
        frags[totF].flag_2=flag;
        strcpy(frags[totF].chr_2, chr);
		frags[totF].st_2=start;
     	strcpy(frags[totF].cigar_2, cigar);
		frags[totF].nreads=2;
	}   
}
