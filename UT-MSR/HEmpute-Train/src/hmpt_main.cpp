#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "hmpt_annot_region_tools.h"
#include "hmpt_variation_tools.h"
#include "hmpt_imputation_utils.h"
#include "hmpt_nomenclature.h"
#include "hmpt_genomics_coords.h"
#include "hmpt_nucleotide.h"
#include "hmpt_ansi_string.h"
#include "hmpt_utils.h"
#include <ctype.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <algorithm>

bool __DUMP_VARIATION_TOOLS_MAIN_MSGS__ = false;

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Imputation Model Training:\n\
	-train_vicinity_based_LMSE_imputation_model_multithreaded\n\
Genotype Extraction:\n\
	-extract_genotype_signals_per_VCF [1kG VCF file path]: This extracts genotype signal into one file.\n\
	-dump_plain_geno_signal_regions\n\
	-extract_genotype_signals_per_subsample_list [Genotype signal regions BED file path] [Sample ids list fp] [Subsample ids list fp]\n\
	-extract_genotype_signals_per_region_list\n\
	-binarize_genotype_signals_BED\n\
	-convert_haplocoded_2_genocoded\n", argv[0]);
		exit(0);
	}

	clock_t start_c = clock();
	time_t start_t = time(NULL);

	if (strcmp(argv[1], "-train_vicinity_based_LMSE_imputation_model_multithreaded") == 0)
	{
		if (argc != 15)
		{
			fprintf(stderr, "USAGE: %s -train_vicinity_based_LMSE_imputation_model_multithreaded [Start posn.] [End posn.] \
[# tag vars per target variant's left/right] \
[Tag genotype matrix file path (text)] [Tag variant matrix sample id's list file path] \
[Target genotype matrix file path (text)] [Target variant matrix sample id's list file path] \
[Testing Tag genotype matrix file path(text)][Testing Tag variant matrix sample id's list file path] \
[Testing Target genotype matrix file path(text)][Testing Target variant matrix sample id's list file path] [Model output directory] [# threads]\n", argv[0]);
			exit(0);
		}

		int start_pos = atoi(argv[2]);
		int end_pos = atoi(argv[3]);
		int n_tags_vars_per_side = atoi(argv[4]);
		char* tag_genotype_matrix_fp = argv[5];
		char* tag_sample_ids_list_fp = argv[6];
		char* target_genotype_matrix_fp = argv[7];
		char* target_sample_ids_list_fp = argv[8];

		char* testing_tag_genotype_matrix_fp = argv[9];
		char* testing_tag_sample_ids_list_fp = argv[10];
		char* testing_target_genotype_matrix_fp = argv[11];
		char* testing_target_sample_ids_list_fp = argv[12];

		char* op_dir = argv[13];

		int n_threads = atoi(argv[14]);

		train_vicinity_based_LMSE_imputation_model_multithreaded(start_pos, end_pos, n_threads,
			n_tags_vars_per_side,
			tag_genotype_matrix_fp,
			tag_sample_ids_list_fp,
			target_genotype_matrix_fp,
			target_sample_ids_list_fp,
			testing_tag_genotype_matrix_fp,
			testing_tag_sample_ids_list_fp,
			testing_target_genotype_matrix_fp,
			testing_target_sample_ids_list_fp,
			op_dir);
	} // -train_vicinity_based_LMSE_imputation_model_multithreaded option.
	else if (strcmp(argv[1], "-dump_plain_geno_signal_regions") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -dump_plain_geno_signal_regions [Genotype signals matBED file path] [sample ids list file path] [Save regions only? (0/1)] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		bool regs_only_flag = (argv[4][0] == '1');
		char* op_fp = argv[5];
		
		vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		dump_geno_sig_regs_plain(geno_sig_regs, sample_ids, regs_only_flag, op_fp);
	} // -dump_plain_geno_signal_regions option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_region_list") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_region_list [Genotype signals BED file path] [sample ids list file path] [BED file with regions] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* regions_BED_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_region_list(geno_sig_regs_BED_fp, sample_ids_list_fp, regions_BED_fp, op_fp);
	} // -extract_genotype_signals_per_region_list option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_subsample_list") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_subsample_list [Genotype signals BED file path] [sample ids list file path] [Subsample ids list file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* subsample_ids_list_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_subsample_list(geno_sig_regs_BED_fp, sample_ids_list_fp, subsample_ids_list_fp, op_fp);
	}
	else if (strcmp(argv[1], "-extract_genotype_signals_per_VCF") == 0)
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_VCF [VCF file path] [VCF sample ids list file path] \
[Variant regions BED file path] [chromosome id to process] \
[Match region names? (0/1)] \
[Haplotype encoded output? (0/1)] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* vcf_fp = argv[2];
		char* vcf_sample_ids_list_fp = argv[3];
		char* var_regions_BED_fp = argv[4];
		char* chr_id_2_process = argv[5];
		bool match_ref_alleles_flag = false;
		char* bin_seq_dir = NULL;
		bool match_region_names_flag = (argv[6][0] == '1');
		bool haplotype_specific_encoding = (argv[7][0] == '1');
		char* op_fp = argv[8];

		extract_genotype_signals_per_VCF(vcf_fp,
											vcf_sample_ids_list_fp,
											var_regions_BED_fp,
											chr_id_2_process,
											bin_seq_dir, 
											match_ref_alleles_flag, 
											match_region_names_flag,
											haplotype_specific_encoding,
											op_fp);
	} // -extract_genotype_signals_per_1kg_VCF option.
	else if (strcmp(argv[1], "-binarize_genotype_signals_BED") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -binarize_genotype_signals_BED [Genotype signals BED file path] [VCF sample ids list file path (Use EpiLeak to extract)] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* variant_geno_sig_regs_BED_fp = argv[2];
		char* vcf_sample_ids_list_fp = argv[3];
		char* op_fp = argv[4];

		vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
		fprintf(stderr, "Loaded %d sample ids.\n", vcf_sample_ids->size());

		binarize_variant_genotype_signal_regions(NULL, variant_geno_sig_regs_BED_fp, vcf_sample_ids, op_fp);
	} // -binarize_genotype_signals_BED option.
	else if (strcmp(argv[1], "-convert_haplocoded_2_genocoded") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -convert_haplocoded_2_genocoded [Haplocoded genotype matbed file path] \
[Haplocoded genotype matrix sample ids list path] \
[Output genotype signal regions bed file path]\n", argv[0]);
			exit(0);
		}

		char* haplocoded_geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_matbed_fp = argv[4];

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
		vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_geno_sig_regs_fp, sample_ids_list_fp);

		fprintf(stderr, "Loaded %d haplocoded genotype regions for %d samples.\n", haplocoded_geno_regs->size(), sample_ids->size());

		for (int i_reg = 0; i_reg < haplocoded_geno_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
			char* haplocoded_geno_sigs = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < sample_ids->size(); i_s++)
			{
				int genocoded_geno = get_genotype_per_haplocoded_genotype(haplocoded_geno_sigs[i_s]);
				haplocoded_geno_sigs[i_s] = genocoded_geno;
			} // i_s loop.
		} // i_reg loop.

		  // Save.
		fprintf(stderr, "Saving to %s.\n", op_matbed_fp);
		binarize_variant_genotype_signal_regions(haplocoded_geno_regs, NULL, sample_ids, op_matbed_fp);
	} // -convert_haplocoded_2_genocoded option.
	
	else 
	{
	}

	FILE* f_beacon = open_f("beacon.log", "a");
	clock_t end_c = clock();
	time_t end_t = time(NULL);
	fprintf(f_beacon, "%s finished option \"%s\" in %d (%d) seconds.\n", argv[0], argv[1], (int)(end_t - start_t), (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fprintf(stderr, "%s finished option \"%s\" in %d (%d) seconds.\n", argv[0], argv[1], (int)(end_t - start_t), (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fclose(f_beacon);

}