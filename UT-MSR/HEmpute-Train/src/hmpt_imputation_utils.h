#ifndef __IMPUTATION_UTILS__
#define __IMPUTATION_UTILS__

struct t_IMPUTE2_op_col_info
{
	int name_col_i;
	int posn_col_i;
	int geno_prob_starting_0_based_col_i;
};

enum { MATCH_BY_NAME, MATCH_BY_START_POSN, N_MATCH_BY_IDS };
char** get_match_by_identifiers();

void compute_imputation_stats_per_5col_genotype_probs_w_PR_curve_info(char* imputed_chr_id,
	char* imputed_5col_genotype_probs_fp,
	char* imputed_5col_genotype_probs_sample_ids_list_fp,
	char* known_genotype_regs_fp,
	char* known_genotype_regs_sample_ids_list_fp,
	char* match_by_str,
	double min_score_threshold,
	double max_score_threshold,
	double score_delta,
	char* op_stats_fp);

void compute_imputation_stats_per_5col_genotype_probs(char* imputed_chr_id,
	char* imputed_5col_genotype_probs_fp,
	char* imputed_5col_genotype_probs_sample_ids_list_fp,
	char* known_genotype_regs_fp,
	char* known_genotype_regs_sample_ids_list_fp,
	char* match_by_str,
	char* op_stats_fp);

void train_vicinity_based_LMSE_imputation_model_multithreaded(int start_pos, int end_pos,
	int n_threads,
	int n_tags_vars_per_side,

	// The training data.
	char* tag_genotype_matrix_fp,
	char* tag_sample_ids_list_fp,
	char* target_genotype_matrix_fp,
	char* target_sample_ids_list_fp,

	// The testing data.
	char* testing_tag_genotype_matrix_fp,
	char* testing_tag_sample_ids_list_fp,
	char* testing_target_genotype_matrix_fp,
	char* testing_target_sample_ids_list_fp,

	char* op_dir);

void train_vicinity_based_LMSE_imputation_model(int start_pos, int end_pos,
	int n_tags_vars_per_side,

	// The training data.
	char* tag_genotype_matrix_fp,
	char* tag_sample_ids_list_fp,
	char* target_genotype_matrix_fp,
	char* target_sample_ids_list_fp,

	// The testing data.
	char* testing_tag_genotype_matrix_fp,
	char* testing_tag_sample_ids_list_fp,
	char* testing_target_genotype_matrix_fp,
	char* testing_target_sample_ids_list_fp,

	char* op_dir);

void extract_genotype_probability_info_per_GT_entry_id_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	char* target_geno_sig_regs_fp,
	char* target_sample_ids_list_fp,
	char* GT_entry_id,
	char* op_prefix);

void extract_genotype_probability_info_per_IMPUTE2_probabilities_output(char* IMPUTE2_output_fp,
	t_IMPUTE2_op_col_info* IMPUTE2_op_col_info,
	char* target_geno_sig_regs_fp,
	char* target_sample_ids_list_fp,
	char* chrom,
	char* op_prefix);

void extract_genotyped_values_per_baseline_regression_model_continuous_genotypes_output(char* baseline_model_output_file_path,
	char* input_genotype_matrix_sample_ids_list_fp,
	char* output_genotype_matrix_bed_fp);

void extract_IMPUTE2_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp,
	char* g_option_fp,
	char* strand_g_option_fp);

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp);

void extract_MACH_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* d_option_fp,
	char* p_option_fp,
	char* h_option_fp,
	char* s_option_fp);

void compute_imputation_stats_per_IMPUTED_genotypes(char* imputed_genotyped_regs_fp, char* imputed_genotyped_regs_sample_ids_list_fp,
	char* known_genotype_regs_fp, char* known_genotype_regs_sample_ids_list_fp);

void extract_genotyped_values_per_IMPUTE2_probabilities_output(char* IMPUTE2_output_fp, char* input_genotype_matrix_sample_ids_list_fp, char* output_genotype_matrix_bed_fp,
	double min_prob_2_genotype,
	t_IMPUTE2_op_col_info* IMPUTE2_op_col_info,
	char* chrom);

#endif // __IMPUTATION_UTILS__
