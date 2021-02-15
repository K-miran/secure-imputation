#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmpt_variation_tools.h"
#include "hmpt_genomics_coords.h"
#include "hmpt_annot_region_tools.h"
#include "hmpt_utils.h"
#include "hmpt_genomics_coords.h"
#include "hmpt_ansi_string.h"
#include "hmpt_nucleotide.h"
#include <string.h>
#include <ctype.h>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmpt_annot_region_tools.h"
#include "hmpt_nomenclature.h"
#include "hmpt_ansi_string.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

bool __DUMP_VARIATION_TOOLS_MSGS__ = false;

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

vector<t_annot_region*>* load_variant_signal_regions_wrapper(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp)
{
	vector<t_annot_region*>* var_sig_regs = NULL;

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	if (sample_ids == NULL)
	{
		fprintf(stderr, "Could not load the sample ids from %s", sample_ids_list_fp);
		exit(0);
	}

	if (t_string::compare_strings(geno_sig_regs_BED_fp, "stdin") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed") || 
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt.gz") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed.gz"))
	{
		var_sig_regs = load_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else if (t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed.gz"))
	{
		var_sig_regs = load_binarized_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else
	{
		fprintf(stderr, "%s(%d): Could not identify variant signal regions file type: %s\n", __FILE__, __LINE__, geno_sig_regs_BED_fp);
		exit(0);
	}

	t_string::clean_string_list(sample_ids);

	return(var_sig_regs);
}

void extract_genotype_signals_per_region_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* regions_BED_fp, char* op_fp)
{
	vector<t_annot_region*>* roi_regs = load_BED(regions_BED_fp);
	for (int i_reg = 0; i_reg < roi_regs->size(); i_reg++)
	{
		roi_regs->at(i_reg)->data = NULL;
	} // i_reg loop.
	fprintf(stderr, "Extracting genotype signals for %d ROIs.\n", roi_regs->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", geno_sig_regs->size(), sample_ids->size());

	fprintf(stderr, "Intersecting ROI's with genotype signal regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(roi_regs, geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_roi_reg = int_info->src_reg;
		t_annot_region* cur_geno_sig_reg = int_info->dest_reg;

		if (t_string::compare_strings(cur_roi_reg->chrom, cur_geno_sig_reg->chrom) &&
			cur_roi_reg->start == cur_geno_sig_reg->start &&
			cur_roi_reg->end == cur_geno_sig_reg->end)
		{
			if (cur_roi_reg->data != NULL)
			{
				fprintf(stderr, "Sanity check failed: Already assigned: %s:%d-%d\n", cur_geno_sig_reg->chrom, cur_geno_sig_reg->start, cur_geno_sig_reg->end);
				exit(0);
			}

			cur_roi_reg->data = cur_geno_sig_reg->data;
		}
	} // i_int loop.

	// Count the number of ROI with signal regions.
	vector<t_annot_region*>* roi_regs_w_signals = new vector<t_annot_region*>();
	int n_roi_w_signal = 0;
	for (int i_reg = 0; i_reg < roi_regs->size(); i_reg++)
	{
		if (roi_regs->at(i_reg)->data != NULL)
		{
			roi_regs_w_signals->push_back(roi_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Matched genotype signals to %d regions.\n", roi_regs_w_signals->size());

	// Save the regions with signals on them.
	binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
}

// Following function extracts the genotype signals for samples from an input VCF file. It focuses on a chromosome and a subset of regions to decrease memory usage.
void extract_genotype_signals_per_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
										char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
										char* var_regions_BED_fp,			// Regions to focus on while extracting.
										char* chr_id_2_process,				// Chromosome to process.
										char* bin_seq_dir,					// Sequences directory.
										bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
										bool match_region_names_flag,		// This is a flag to enforce matching of region names.
										bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
										char* op_fp)						// Output file path.
{
	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(0);
	}

	vector<t_annot_region*>* all_var_regions = load_BED(var_regions_BED_fp);

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", vcf_sample_ids->size());

	vector<char*>* chr_ids = get_chr_ids(all_var_regions);

	vector<t_annot_region*>* var_regions = NULL;
	int chrom_i_2_process = t_string::get_i_str(chr_ids, chr_id_2_process);
	if (chrom_i_2_process == chr_ids->size())
	{
		var_regions = all_var_regions;
	}
	else
	{
		var_regions = get_regions_per_chromosome(all_var_regions, chr_ids->at(chrom_i_2_process));
	}

	chr_ids = get_chr_ids(var_regions);

	fprintf(stderr, "Extracting signals on %d regions.\n", var_regions->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}
	else
	{
	}

	char** per_chr_seq = new char*[chr_ids->size() + 2];

	if (match_ref_alleles_flag)
	{
		//for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		//{
		//	int l_chrom = 0;
		//	fprintf(stderr, "Loading %s.\n", chr_ids->at(i_chr));
		//	char cur_chr_seq_fp[1000];
		//	sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));

		//	if (check_file(cur_chr_seq_fp))
		//	{
		//		//per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
		//		per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
		//	}
		//	else
		//	{
		//		sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));

		//		if (check_file(cur_chr_seq_fp))
		//		{
		//			per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
		//		}
		//		else
		//		{
		//			fprintf(stderr, "Could not load the sequence for %s\n", chr_ids->at(i_chr));
		//			exit(0);
		//		}
		//	}
		//} // i_chr loop.
	}
	else
	{
		fprintf(stderr, "Skipping ref genome loading.\n");
	}
	  // Set the genotype signal array on all regions.
	for (int i_v = 0; i_v < var_regions->size(); i_v++)
	{
		char* pooled_var_alleles = var_regions->at(i_v)->name;

		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = NULL;
		cur_reg_info[1] = NULL;
		var_regions->at(i_v)->data = cur_reg_info;
	} // i_v loop.

	FILE* f_vcf = open_f(vcf_fp, "r");

	char* buff = new char[100 * 1000];
	vector<t_annot_region*>* vcf_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if (vcf_regs->size() % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region.           \r", vcf_regs->size());
		}

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT
		char chrom[100];
		char posn_str[100];
		char id[100];
		char ref[100];
		char alt[100];
		char qual[100];
		char filter[100];
		char info[10000];
		char format[100];

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(format, buff);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int GT_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];
		while (t_string::get_next_token(format, cur_tok, l_format, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, "GT"))
			{
				GT_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "Found GT entry @ %d\n", GT_entry_i);
		}

		if (GT_entry_i == -1)
		{
			fprintf(stderr, "Could not find GT entry in %s, skipping\n", cur_line);
			delete[] cur_line;
			continue;
		}

		if (GT_entry_i != 0)
		{
			fprintf(stderr, "Format string is not as expected: %s\n", format);
		}
		
		// Make sure we have the chromosome id match and we are using alleles.
		char* cur_var_geno_sig = new char[vcf_sample_ids->size() + 2];

		// Initialize everything to not set: This contains the variant signals for each individual.
		for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
		{
			cur_var_geno_sig[i_s] = -1;
		} // i_s loop.

		bool correctly_parsed_all_genotypes = true;
		for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);

			//if (buff[3] != 0)
			//{
			//	fprintf(stderr, "Failed to read the %d. genotype correctly for line:%s\n", i_s, cur_line);
			//	exit(0);
			//}

			// Check if all the variants are bi-allelic, only.
			if ((buff[0] != '.' && buff[0] != '.' && buff[0] != '0' && buff[0] != '1') ||
				(buff[1] != '|' && buff[1] != '/') ||
				(buff[2] != '0' && buff[2] != '1' && buff[2] != '.' && buff[2] != '.'))
			{
				correctly_parsed_all_genotypes = false;
				fprintf(stderr, "Failed to read the %d. genotype correctly for entry %s in line: %s\n", i_s, buff, cur_line);
				//exit(0);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[1] == '|')
				{
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
			}
			else
			{
				cur_var_geno_sig[i_s] = 0;

				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[0] == '1' && buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 2;
				}
				else if (buff[0] == '1' || buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 1;
				}
			}

			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(0);
		}

		if (correctly_parsed_all_genotypes)
		{
			// Intersect with the variant regions.
			t_annot_region* cur_vcf_reg = get_empty_region();
			cur_vcf_reg->chrom = t_string::copy_me_str(chrom);
			normalize_chr_id(cur_vcf_reg->chrom);
			cur_vcf_reg->start = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);

			int l_ref_allele = t_string::string_length(ref);
			cur_vcf_reg->end = cur_vcf_reg->start + l_ref_allele - 1;
			cur_vcf_reg->name = t_string::copy_me_str(id);
			cur_vcf_reg->strand = '+';

			// Set the signal to data.
			void** vcf_reg_info = new void*[2];

			// 1st data: signals.
			vcf_reg_info[0] = cur_var_geno_sig;

			// 2nd data: alleles.
			char** ref_alt_alleles = new char*[2];
			ref_alt_alleles[0] = t_string::copy_me_str(ref);
			ref_alt_alleles[1] = t_string::copy_me_str(alt);
			vcf_reg_info[1] = ref_alt_alleles;

			// If matching ref alleles is requested, check the alleles.
			if (l_ref_allele < 100 &&
				match_ref_alleles_flag)
			{
				int i_chr_cur_line = t_string::get_i_str(chr_ids, cur_vcf_reg->chrom);
				for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
				{
					if (toupper(per_chr_seq[i_chr_cur_line][i]) != toupper(ref[i - cur_vcf_reg->start]))
					{
						fprintf(stderr, "Could not match the reference allele to chromosome sequence @ %s:%d: %c, %c\n",
							cur_vcf_reg->chrom, cur_vcf_reg->start,
							toupper(per_chr_seq[i_chr_cur_line][i]), toupper(ref[i - cur_vcf_reg->start]));

						exit(0);
					}
				} // i loop.

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Match @ %s:%d-%d (%s): \n%s\n", cur_vcf_reg->chrom, cur_vcf_reg->start, cur_vcf_reg->end, cur_vcf_reg->name, ref);
					for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
					{
						fprintf(stderr, "%c", toupper(per_chr_seq[i_chr_cur_line][i]));
					} // i loop.
					fprintf(stderr, "\n");
				}
			} // allele check is done for variants smaller than 100 bp.

			cur_vcf_reg->data = vcf_reg_info;

			vcf_regs->push_back(cur_vcf_reg);
		} // genotype parsing check.

		delete[] cur_line;
	} // vcf file reading loop.

	vector<t_annot_region*>* intersects = intersect_annot_regions(var_regions, vcf_regs, true);
	fprintf(stderr, "Processing %d intersections using coordinate matching.               \n", intersects->size());

	if (match_region_names_flag)
	{
		fprintf(stderr, "Enforcing substring based name matching.\n");
	}

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* pooled_var_reg = cur_int_info->src_reg;
		t_annot_region* vcf_reg = cur_int_info->dest_reg;

		bool coord_check_pass = true;
		if (pooled_var_reg->start != vcf_reg->start ||
			pooled_var_reg->end != vcf_reg->end)
		{
			coord_check_pass = false;
		}

		// Do we need to have name matching?
		bool name_check_pass = true;
		if (match_region_names_flag &&
			pooled_var_reg->name != NULL &&
			vcf_reg->name != NULL)
		{
			int i_char = 0;

			// If there is a dot for the name, skip.
			if (!t_string::compare_strings(pooled_var_reg->name, ".") &&
				!t_string::compare_strings(vcf_reg->name, ".") &&
				!t_string::compare_substrings_ci(pooled_var_reg->name, vcf_reg->name, i_char) &&
				!t_string::compare_substrings_ci(vcf_reg->name, pooled_var_reg->name, i_char))
			{
				name_check_pass = false;
			}
		}

		if (name_check_pass &&
			coord_check_pass)
		{
			int i_chr = t_string::get_i_str(chr_ids, pooled_var_reg->chrom);
			char* chr_seq = per_chr_seq[i_chr];

			// Get the genotype signal for the vcf region.
			void** vcf_reg_info = (void**)(vcf_reg->data);
			char* cur_vcf_reg_geno_sig = (char*)(vcf_reg_info[0]);

			// Set the genotype signal.
			void** cur_var_reg_info = (void**)(pooled_var_reg->data);
			cur_var_reg_info[0] = cur_vcf_reg_geno_sig;
		}

		delete cur_int_info;
	} // i_int loop.
	delete_annot_regions(vcf_regs);

	// Copy the current haplotype alleles.
	for (int i_reg = 0; i_reg < var_regions->size(); i_reg++)
	{
		void** cur_var_reg_info = (void**)(var_regions->at(i_reg)->data);

		if (cur_var_reg_info[0] == NULL)
		{
			char* cur_var_geno_sig = new char[vcf_sample_ids->size()];
			for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
			{
				cur_var_geno_sig[i_s] = -1;
			} // i_s loop.

			cur_var_reg_info[0] = cur_var_geno_sig;
		}
	} // i_reg loop.

	// Save the binarized signals.
	binarize_variant_genotype_signal_regions(var_regions, NULL, vcf_sample_ids, op_fp);

	//// Dump the plain matrix, too.
	//dump_geno_sig_regs_plain(var_regions, vcf_sample_ids, "op.bed");
	
	//// Now dump the variant signal profiles.
	//FILE* f_geno_sig = open_f("op.bed", "w");
	//for (int i_v = 0; i_v < var_regions->size(); i_v++)
	//{
	//	fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
	//		var_regions->at(i_v)->chrom, 
	//		translate_coord(var_regions->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//		translate_coord(var_regions->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
	//		var_regions->at(i_v)->name);

	//	void** cur_var_reg_info = (void**)(var_regions->at(i_v)->data);
	//	int* cur_var_signal = (int*)(cur_var_reg_info[0]);
	//	for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
	//	{
	//		fprintf(f_geno_sig, "\t%d", cur_var_signal[i_s]);
	//	} // i_s loop.

	//	fprintf(f_geno_sig, "\n");
	//} // i_v loop.
	//fclose(f_geno_sig);
}

int get_genotype_per_haplocoded_genotype(char haplocoded_geno)
{
	int all1 = get_allele_per_haplotype(haplocoded_geno, 0);
	int all2 = get_allele_per_haplotype(haplocoded_geno, 1);
	return(all1 + all2);
}

int get_allele_per_haplotype(char geno, int hap_i)
{
	int allele = ((geno & (1 << hap_i)) >> hap_i);
	return(allele);
}

void dump_geno_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, bool dump_geno_sig_regs_plain, const char* op_fp)
{
	if (dump_geno_sig_regs_plain)
	{
		fprintf(stderr, "Saving regions only.\n");
	}

	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);

		if (!dump_geno_sig_regs_plain)
		{
			void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
			char* cur_var_signal = (char*)(cur_var_reg_info[0]);
			for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
			{
				fprintf(f_geno_sig, "\t%d", (int)cur_var_signal[i_s]);
			} // i_s loop.
		}

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	fclose(f_geno_sig);
}

void extract_genotype_signals_per_subsample_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* subsample_ids_list_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", geno_sig_regs->size(), sample_ids->size());

	vector<char*>* geno_subsample_ids = buffer_file(subsample_ids_list_fp);
	if (geno_subsample_ids == NULL)
	{
		fprintf(stderr, "Could not load subsample id's list from %s\n", subsample_ids_list_fp);
		exit(0);
	}
	else
	{
		fprintf(stderr, "Extracting genotype signals for %d individuals.\n", geno_subsample_ids->size());
	}

	vector<int>* per_subsample_geno_sample_i = new vector<int>();
	for (int i_ss = 0; i_ss < geno_subsample_ids->size(); i_ss++)
	{
		int cur_ss_sample_i = t_string::get_i_str(sample_ids, geno_subsample_ids->at(i_ss));
		per_subsample_geno_sample_i->push_back(cur_ss_sample_i);
	} // i_sub_s loop.

	vector<t_annot_region*>* geno_sig_regs_w_subsample_signals = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < geno_sig_regs->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Extracting %d. region            \r", i_reg);
		}

		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		char* cur_reg_geno_sigs = (char*)(cur_reg_info[0]);

		// Copy the signals.
		char* cur_copy_reg_geno_sigs = new char[geno_subsample_ids->size() + 2];
		for (int i_ss = 0; i_ss < geno_subsample_ids->size(); i_ss++)
		{			
			if (per_subsample_geno_sample_i->at(i_ss) < sample_ids->size())
			{
				cur_copy_reg_geno_sigs[i_ss] = cur_reg_geno_sigs[per_subsample_geno_sample_i->at(i_ss)];
			}
			else
			{
				cur_copy_reg_geno_sigs[i_ss] = -1;
			}
		} // i_ss loop.

		void** cur_copy_reg_info = new void*[2];		
		cur_copy_reg_info[0] = cur_copy_reg_geno_sigs;

		t_annot_region* copy_reg = duplicate_region(geno_sig_regs->at(i_reg));
		copy_reg->data = cur_copy_reg_info;
		geno_sig_regs_w_subsample_signals->push_back(copy_reg);
	} // i_reg loop.

	// Save the regions.
	binarize_variant_genotype_signal_regions(geno_sig_regs_w_subsample_signals, NULL, geno_subsample_ids, op_fp);
}

// Following can dump the supplied regions or load them then dump.
void binarize_variant_genotype_signal_regions(vector<t_annot_region*>* genotype_signal_regions, char* variant_geno_sig_regs_BED_fp, vector<char*>* geno_sample_ids, const char* op_fp)
{
	if (genotype_signal_regions == NULL)
	{
		if (t_string::compare_strings(variant_geno_sig_regs_BED_fp, "stdin") || check_file(variant_geno_sig_regs_BED_fp))
		{
			genotype_signal_regions = load_variant_genotype_signal_regions(variant_geno_sig_regs_BED_fp, geno_sample_ids);
		}		
		else
		{
			fprintf(stderr, "Signal regions are not supplied and could not load them using %s.\n", variant_geno_sig_regs_BED_fp);
			exit(0);
		}
	}

	vector<char*>* chr_ids = get_chr_ids(genotype_signal_regions);

	FILE* f_op = open_f(op_fp, "wb");
	int n_chrs = chr_ids->size();
	fwrite(&n_chrs, sizeof(int), 1, f_op);
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chr[1000];
		strcpy(cur_chr, chr_ids->at(i_chr));
		fwrite(cur_chr, sizeof(char), 1000, f_op);
	} // i_chr loop.

	int sample_size = geno_sample_ids->size();
	int n_regs = genotype_signal_regions->size();

	fwrite(&sample_size, sizeof(int), 1, f_op);
	fwrite(&n_regs, sizeof(int), 1, f_op);
	for (int i_reg = 0; i_reg < genotype_signal_regions->size(); i_reg++)
	{
		// Write the chromosome index.
		int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		fwrite(&i_chr, sizeof(int), 1, f_op);
		fwrite(&(reg_BED_start), sizeof(int), 1, f_op);
		fwrite(&(reg_BED_end), sizeof(int), 1, f_op);

		// Write the regions's name.
		if (genotype_signal_regions->at(i_reg)->name == NULL)
		{
			genotype_signal_regions->at(i_reg)->name = t_string::copy_me_str(".");
		}

		int l_reg_name_str = t_string::string_length(genotype_signal_regions->at(i_reg)->name);
		fwrite(&l_reg_name_str, sizeof(int), 1, f_op);
		fwrite(genotype_signal_regions->at(i_reg)->name, sizeof(char), l_reg_name_str, f_op);

		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		fwrite(cur_reg_geno_sig, sizeof(char), geno_sample_ids->size(), f_op);
	} // i_reg loop.

	fprintf(stderr, "\nClosing the file.\n");
	close_f(f_op, op_fp);
	fprintf(stderr, "Finished writing. Doing a sanity check.\n");

	bool perform_check = true;

	if (!perform_check)
	{
		fprintf(stderr, "Skipping sanity check.\n");
		return;
	}

	clock_t cur_time = clock();
	vector<t_annot_region*>* loaded_sig_regs = load_binarized_variant_genotype_signal_regions(op_fp, geno_sample_ids);
	fprintf(stderr, "Loaded in %.4f seconds.\n", ((double)(clock() - cur_time)) / CLOCKS_PER_SEC);

	if (loaded_sig_regs->size() != genotype_signal_regions->size())
	{
		fprintf(stderr, "Dumping/Loading failed.\n");
		exit(0);
	}

	for (int i_reg = 0; i_reg < loaded_sig_regs->size(); i_reg++)
	{
		if (loaded_sig_regs->at(i_reg)->start != genotype_signal_regions->at(i_reg)->start ||
			loaded_sig_regs->at(i_reg)->end != genotype_signal_regions->at(i_reg)->end)
		{
			fprintf(stderr, "Non-match: %s:%d-%d, %s:%d-%d\n", 
					loaded_sig_regs->at(i_reg)->chrom, loaded_sig_regs->at(i_reg)->start, loaded_sig_regs->at(i_reg)->end, 
					genotype_signal_regions->at(i_reg)->chrom, genotype_signal_regions->at(i_reg)->start, genotype_signal_regions->at(i_reg)->end);
			exit(0);
		}

		// Go over the genotype signals.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

		void** cur_loaded_reg_info = (void**)(loaded_sig_regs->at(i_reg)->data);
		char* cur_loaded_reg_geno_sig = (char*)(cur_loaded_reg_info[0]);

		// Write the signal values.
		for (int i_s = 0; i_s < geno_sample_ids->size(); i_s++)
		{
			if (cur_reg_geno_sig[i_s] != cur_loaded_reg_geno_sig[i_s])
			{
				fprintf(stderr, "Non-match: %d. sample genotype: %d, %d\n",
						cur_reg_geno_sig[i_s], cur_loaded_reg_geno_sig[i_s]);
				exit(0);
			}
		} // i_s loop.
	} // i_reg loop.

	fprintf(stderr, "Check success!\n");
}

vector<t_annot_region*>* load_binarized_variant_genotype_signal_regions(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids)
{
	FILE* f_bin_geno_sig_regs = open_f(bin_geno_sig_bed_fp, "rb");

	// Load the chromosomes.
	int n_chrs = 0;
	fread(&n_chrs, sizeof(int), 1, f_bin_geno_sig_regs);
	vector<char*>* chr_ids = new vector<char*>();
	for (int i_chr = 0; i_chr < n_chrs; i_chr++)
	{
		char cur_chr[1000];
		fread(cur_chr, sizeof(char), 1000, f_bin_geno_sig_regs);
		chr_ids->push_back(t_string::copy_me_str(cur_chr));
	} // i_chr loop.

	fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());
	int sample_size = 0;
	fread(&sample_size, sizeof(int), 1, f_bin_geno_sig_regs);
	fprintf(stderr, "Reading sample size of %d.\n", sample_size);

	if (sample_size != geno_sample_ids->size())
	{
		fprintf(stderr, "Sanity check failed: Sample sizes do not match: %d, %d\n", sample_size, geno_sample_ids->size());
		exit(0);
	}

	int n_regs = 0;
	fread(&n_regs, sizeof(int), 1, f_bin_geno_sig_regs);
	fprintf(stderr, "Reading %d regions.\n", n_regs);

	vector<t_annot_region*>* geno_sig_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < n_regs; i_reg++)
	{
		//int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		//int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		//int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		int i_chr = 0;
		int reg_BED_start = 0;
		int reg_BED_end = 0;
		fread(&i_chr, sizeof(int), 1, f_bin_geno_sig_regs);
		fread(&(reg_BED_start), sizeof(int), 1, f_bin_geno_sig_regs);
		fread(&(reg_BED_end), sizeof(int), 1, f_bin_geno_sig_regs);

		// Read the region's name.
		int l_reg_name_str = 0;
		fread(&l_reg_name_str, sizeof(int), 1, f_bin_geno_sig_regs);
		char* cur_reg_name = new char[l_reg_name_str + 2];
		memset(cur_reg_name, 0, sizeof(char) * (l_reg_name_str + 2));
		fread(cur_reg_name, sizeof(char), l_reg_name_str, f_bin_geno_sig_regs);

		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
		reg->start = translate_coord(reg_BED_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		reg->end = translate_coord(reg_BED_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		reg->strand = '+';
		reg->name = cur_reg_name;

		void** cur_reg_info = new void*[2];
		char* cur_reg_geno_sig = new char[sample_size + 2];
		fread(cur_reg_geno_sig, sizeof(char), sample_size, f_bin_geno_sig_regs);

		cur_reg_info[0] = cur_reg_geno_sig;
		cur_reg_info[1] = NULL;

		reg->data = cur_reg_info;

		geno_sig_regs->push_back(reg);
	} // i_reg loop.

	// Close the file.
	close_f(f_bin_geno_sig_regs, bin_geno_sig_bed_fp);

	return(geno_sig_regs);
}

vector<t_annot_region*>* load_variant_genotype_signal_regions(char* link_variant_genotype_signal_fp, vector<char*>* geno_sample_ids)
{
	vector<t_annot_region*>* link_var_regs = new vector<t_annot_region*>();
	char tok_buff[1000];
	FILE* f_link_variant_genotype_signal = open_f(link_variant_genotype_signal_fp, "r");

	if (f_link_variant_genotype_signal == NULL)
	{
		fprintf(stderr, "Could not open %s\n", link_variant_genotype_signal_fp);
		exit(0);
	}

	int i_reg = 0;
	while (1)
	{
		char* cur_reg_line = getline(f_link_variant_genotype_signal);
		if (cur_reg_line == NULL)
		{
			break;
		}
		else
		{
			i_reg++;
		}

		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Parsing genotype signal for %d. region.                  \r", i_reg);
		}

		//char* cur_reg_line = (char*)(link_var_regs->at(i_reg)->data);
		//t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
		char* cur_var_geno_signal = new char[geno_sample_ids->size()];

		t_annot_region* cur_reg = get_empty_region();

		// Copy the name.
		int i_cur_char = 0;
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->chrom = t_string::copy_me_str(tok_buff);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->start = translate_coord(atoi(tok_buff), BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->end = translate_coord(atoi(tok_buff), BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_reg->strand = '+';

		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->name = t_string::copy_me_str(tok_buff);

		//fprintf(stderr, "%s\t%d\t%d: %s            \r", link_var_regs->at(i_reg)->chrom, link_var_regs->at(i_reg)->start, link_var_regs->at(i_reg)->end, link_var_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < geno_sample_ids->size(); i_s++)
		{
			if (t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char) == false)
			{
				fprintf(stderr, "Could not parse the %d. sample's genotype signal entry.\n", i_s);
				exit(0);
			}

			cur_var_geno_signal[i_s] = (char)(atoi(tok_buff));
		} // i_tok loop.

		void** link_var_reg_info = new void*[3];
		link_var_reg_info[0] = cur_var_geno_signal;
		link_var_reg_info[1] = NULL;

		cur_reg->data = link_var_reg_info;

		link_var_regs->push_back(cur_reg);
	} // i_reg loop.

	close_f(f_link_variant_genotype_signal, link_variant_genotype_signal_fp);

	return(link_var_regs);
}

