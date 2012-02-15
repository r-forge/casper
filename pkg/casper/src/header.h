typedef struct read_t {
	int len;
	int st_1;
	int end_1;
	int flag_1;
	int len_1;
	char *cigar_1;
	char *chr_1;
	int st_2;
	int end_2;
	int flag_2;
	int len_2;
	char *cigar_2;
	char *chr_2;
	char *qname;
	int nreads;
} read_t;

int verbose;
