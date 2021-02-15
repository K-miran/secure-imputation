#ifndef __NUCLEOTIDE__
#define __NUCLEOTIDE__

// Basic functionality about the nucleotide naming, pairing, numbering, etc.
bool check_rna_pairing(char nuc1, char nuc2);
char num_2_nuc(int num);
int nuc_2_num(char nuc);
bool is_valid_nuc_char(char nuc);
char get_transcribed_rna_nuc_per_dna_nuc(char dna_nuc);
char get_transcribing_dna_nuc_per_rna_nuc(char rna_nuc);

char get_complementary_dna_nuc_per_dna_nuc(char dna_nuc);
char get_complementary_rna_nuc_per_rna_nuc(char rna_nuc);

int* get_per_char_as_nuc_2_num_coding_array();

void reverse_complement_seq(char* seq);

// Convert nucleotide symbols into indices: XACGUI -> 012345 
// refer to IUPAC nucleotide symbols for more information:
// http://www.mun.ca/biochem/courses/3107/symbols.html
void map_nuc_IUPAC_code(char raw_nuc, 
								char& trans_nuc, 
								int& num, 
								bool& force_unpaired);

//TTT	Phe		TCT	Ser		TAT	Tyr		TGT	Cys
//TTC	Phe		TCC	Ser		TAC	Tyr		TGC	Cys
//TTA	Leu		TCA	Ser		TAA	STOP	TGA	STOP
//TTG	Leu		TCG	Ser		TAG	STOP	TGG	Trp
//CTT	Leu		CCT	Pro		CAT	His		CGT	Arg
//CTC	Leu		CCC	Pro		CAC	His		CGC	Arg
//CTA	Leu		CCA	Pro		CAA	Gln		CGA	Arg
//CTG	Leu		CCG	Pro		CAG	Gln		CGG	Arg
//ATT	Ile		ACT	Thr		AAT	Asn		AGT	Ser
//ATC	Ile		ACC	Thr		AAC	Asn		AGC	Ser
//ATA	Ile		ACA	Thr		AAA	Lys		AGA	Arg
//ATG	Met		ACG	Thr		AAG	Lys		AGG	Arg
//GTT	Val		GCT	Ala		GAT	Asp		GGT	Gly
//GTC	Val		GCC	Ala		GAC	Asp		GGC	Gly
//GTA	Val		GCA	Ala		GAA	Glu		GGA	Gly
//GTG	Val		GCG	Ala		GAG	Glu		GGG	Gly 
static char codons[64][4] = {"AAA",
"AAC",
"AAG",
"AAT",
"ACA",
"ACC",
"ACG",
"ACT",
"AGA",
"AGC",
"AGG",
"AGT",
"ATA",
"ATC",
"ATG",
"ATT",
"CAA",
"CAC",
"CAG",
"CAT",
"CCA",
"CCC",
"CCG",
"CCT",
"CGA",
"CGC",
"CGG",
"CGT",
"CTA",
"CTC",
"CTG",
"CTT",
"GAA",
"GAC",
"GAG",
"GAT",
"GCA",
"GCC",
"GCG",
"GCT",
"GGA",
"GGC",
"GGG",
"GGT",
"GTA",
"GTC",
"GTG",
"GTT",
"TAA",
"TAC",
"TAG",
"TAT",
"TCA",
"TCC",
"TCG",
"TCT",
"TGA",
"TGC",
"TGG",
"TGT",
"TTA",
"TTC",
"TTG",
"TTT"};

static char AAs[64][5] = {"Lys",
"Asn",
"Lys",
"Asn",
"Thr",
"Thr",
"Thr",
"Thr",
"Arg",
"Ser",
"Arg",
"Ser",
"Ile",
"Ile",
"Met",
"Ile",
"Gln",
"His",
"Gln",
"His",
"Pro",
"Pro",
"Pro",
"Pro",
"Arg",
"Arg",
"Arg",
"Arg",
"Leu",
"Leu",
"Leu",
"Leu",
"Glu",
"Asp",
"Glu",
"Asp",
"Ala",
"Ala",
"Ala",
"Ala",
"Gly",
"Gly",
"Gly",
"Gly",
"Val",
"Val",
"Val",
"Val",
"STOP",
"Tyr",
"STOP",
"Tyr",
"Ser",
"Ser",
"Ser",
"Ser",
"STOP",
"Cys",
"Trp",
"Cys",
"Leu",
"Phe",
"Leu",
"Phe"};

char* get_AA_code_per_codon(char* codon);
int codon_2_num(char* codon);

char get_dna_pair(char nuc);

#endif // __NUCLEOTIDE__