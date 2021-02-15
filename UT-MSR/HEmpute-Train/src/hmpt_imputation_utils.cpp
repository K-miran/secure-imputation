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
#include "hmpt_nucleotide.h"
#include "hmpt_nomenclature.h"
#include "hmpt_ansi_string.h"
#include "hmpt_ansi_thread.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "hmpt_imputation_utils.h"

#ifdef __unix__
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#endif 

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

bool __DUMP_INPUTATION_UTILS_MSGS__ = false;

struct t_LMSE_imputation_thread_info
{
	int which;
	int outof;

	double target_normalizer;
	double tag_normalizer;

	//vector<t_annot_region*>* testing_tag_genotype_signal_regs;
	vector<t_annot_region*>* tag_genotype_signal_regs;
	//vector<t_annot_region*>* testing_target_genotype_signal_regs;
	vector<t_annot_region*>* target_genotype_signal_regs;

	t_restr_annot_region_list* restr_target_genotype_signal_regs;
	t_restr_annot_region_list* restr_tag_genotype_signal_regs;

	int testing_flag;
	int start_pos;
	int end_pos;
	char* op_dir;

	vector<char*>* tag_sample_ids;
	vector<char*>* target_sample_ids;

	vector<char*>* testing_tag_sample_ids;
	vector<char*>* testing_target_sample_ids;

	int n_tags_vars_per_side;
};

void get_stats(double* vals, int n_pts, double& mean, double& var)
{
	// Get the mean.
	double total = 0.0;
	for (int i = 0; i < n_pts; i++)
	{
		total += vals[i];
	} // i loop.

	mean = total / n_pts;

	// Get the variance.
	var = 0.0;
	for (int i = 0; i < n_pts; i++)
	{
		var += (vals[i] - mean) * (vals[i] - mean);
	} // i loop.

	var = var / (n_pts - 1);
}

void get_stats(vector<double>* energies, double& mean, double& std_dev)
{
	mean = 0.0;
	for (int i = 0; i < energies->size(); i++)
	{
		mean += energies->at(i);
	} // i loop.

	mean /= energies->size();

	std_dev = 0.0;
	for (int i = 0; i < energies->size(); i++)
	{
		std_dev += (energies->at(i) - mean) * (energies->at(i) - mean);
	} // i loop.

	  // Do unbiased computation for sample variance estimate.
	std_dev /= (energies->size() - 1);

	std_dev = pow(std_dev, 0.5);
}

static void* LMSE_imputation_thread_callback(void* thread_info_ptr)
{
	t_LMSE_imputation_thread_info* impute_thread_info = (t_LMSE_imputation_thread_info*)(thread_info_ptr);

	// Following are the thread info.
	double target_normalizer = impute_thread_info->target_normalizer;
	double tag_normalizer = impute_thread_info->tag_normalizer;

	//vector<t_annot_region*>* testing_tag_genotype_signal_regs = impute_thread_info->testing_tag_genotype_signal_regs;
	vector<t_annot_region*>* tag_genotype_signal_regs = impute_thread_info->tag_genotype_signal_regs;
	//vector<t_annot_region*>* testing_target_genotype_signal_regs = impute_thread_info->testing_target_genotype_signal_regs;
	vector<t_annot_region*>* target_genotype_signal_regs = impute_thread_info->target_genotype_signal_regs;

	t_restr_annot_region_list* restr_target_genotype_signal_regs = impute_thread_info->restr_target_genotype_signal_regs;
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = impute_thread_info->restr_tag_genotype_signal_regs;

	int testing_flag = impute_thread_info->testing_flag;
	int start_pos = impute_thread_info->start_pos;
	int end_pos = impute_thread_info->end_pos;
	char* op_dir = impute_thread_info->op_dir;

	vector<char*> * tag_sample_ids = impute_thread_info->tag_sample_ids;
	vector<char*>* target_sample_ids = impute_thread_info->target_sample_ids;

	vector<char*>* testing_tag_sample_ids = impute_thread_info->testing_tag_sample_ids;
	vector<char*>* testing_target_sample_ids = impute_thread_info->testing_target_sample_ids;

	int n_tags_vars_per_side = impute_thread_info->n_tags_vars_per_side;
	int which = impute_thread_info->which;
	int outof = impute_thread_info->outof;

	fprintf(stderr, "Started Imputation thread %d/%d\n", which, outof);

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// This is the number of data points in training.
	int training_sample_size = target_sample_ids->size();

	// Allocate the fitting data.
	double xi, yi, ei, chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	gsl_vector *testing_x;

	// SVD associated matrices/vectors.
	gsl_matrix* fullX = gsl_matrix_alloc(training_sample_size, n_params_w_intercept);
	gsl_matrix * svdV = gsl_matrix_alloc(n_params_w_intercept, n_params_w_intercept);
	gsl_vector * svdS = gsl_vector_alloc(n_params_w_intercept);
	gsl_vector * svdWork = gsl_vector_alloc(n_params_w_intercept);

	// Tag usage flag at the variant coordinates; this is necessary to save file and write per variant messages.
	int* tag_usage_flag = new int[n_params_w_intercept + 1];
	memset(tag_usage_flag, 0, sizeof(int) * (n_params_w_intercept + 1));

	char per_var_stats_fp[1000];
	sprintf(per_var_stats_fp, "per_var_stats_thread_%d.txt", which);
	FILE* f_per_var_stats = open_f(per_var_stats_fp, "w");

	char pooled_params_fp[1000];
	sprintf(pooled_params_fp, "pooled_params_%d.txt", which);
	FILE* f_pooled_params = open_f(pooled_params_fp, "w");

	// Process all the target variants.
	for (int target_chr_i = 0;
		target_chr_i < restr_target_genotype_signal_regs->chr_ids->size();
		target_chr_i++)
	{
		fprintf(stderr, "Building the models on chromosome %s (%d/%d)\n", restr_target_genotype_signal_regs->chr_ids->at(target_chr_i), which, outof);

		// Get the chromosome index for testring variant regions.
		int tag_chr_i = t_string::get_i_str(restr_tag_genotype_signal_regs->chr_ids, restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		vector<t_annot_region*>* target_var_regs = restr_target_genotype_signal_regs->regions_per_chrom[target_chr_i];
		vector<t_annot_region*>* tag_var_regs = restr_tag_genotype_signal_regs->regions_per_chrom[tag_chr_i];

		// Following sets up n data points. We set this up over the samples.
		for (int target_var_i = 0; target_var_i < target_var_regs->size(); target_var_i++)
		{
			if (start_pos < 0 ||
				end_pos < 0 ||
				target_var_regs->at(target_var_i)->start >= start_pos && target_var_regs->at(target_var_i)->start <= end_pos)
			{

			}
			else
			{
				continue;
			}

			if (target_var_i % outof != which)
			{
				continue;
			}

			int closest_tag_var_left_i = -1;
			for (int tag_var_i = 0; (tag_var_i + 1)< tag_var_regs->size(); tag_var_i++)
			{
				if (tag_var_regs->at(tag_var_i)->start < target_var_regs->at(target_var_i)->start &&
					tag_var_regs->at(tag_var_i + 1)->start > target_var_regs->at(target_var_i)->start)
				{
					closest_tag_var_left_i = tag_var_i;
					break;
				}
			} // tag_var_i loop.

			if (closest_tag_var_left_i == -1)
			{
				continue;
			}

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Closest tag to target %d: (%d, %d)\n", target_var_regs->at(target_var_i)->start,
					tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);
			}

			// Add the tag variant to the right.
			vector<t_annot_region*>* tag_var_regs_per_cur_target_var = new vector<t_annot_region*>();
			for (int tag_var_i = MAX(0, closest_tag_var_left_i - n_tags_vars_per_side + 1);
				tag_var_i <= closest_tag_var_left_i;
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			for (int tag_var_i = closest_tag_var_left_i + 1;
				tag_var_i <= (closest_tag_var_left_i + n_tags_vars_per_side) && tag_var_i < tag_var_regs->size();
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			  // Check the number of tag SNVs around the target SNVs, make sure they satisfy the requested number.
			if (tag_var_regs_per_cur_target_var->size() != 2 * n_tags_vars_per_side)
			{
				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					fprintf(stderr, "Could not find the correct number of tag variants for %s:%d (%d), will skip.\n",
						target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start,
						tag_var_regs_per_cur_target_var->size());
				}
			}
			else
			{
				// Set the tag usage flag.
				memset(tag_usage_flag, 0, sizeof(int) * n_params_w_intercept);
				int n_tags_2_use = 0;
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					tag_usage_flag[par_i] = 1;
				} // par_i loop.

				  // Fill fullX.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(fullX, sample_i, 0, 1.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						// Set the current entry for the full X matrix.
						double tagg_ij = cur_tag_genotypes[sample_i];

						// Normalize: Use tag normalization factor.
						tagg_ij /= tag_normalizer;

						gsl_matrix_set(fullX, sample_i, tag_dim_i, tagg_ij);
					} // tag_dim_i loop.
				} // sample_i loop.			

				  // https://stackoverflow.com/questions/12304963/using-eigenvalues-to-test-for-singularity-identifying-collinear-columns
				  // Generate the SVD of X.
				gsl_linalg_SV_decomp(fullX, svdV, svdS, svdWork);

				// Look for the low singular values that can create problems.
				double top_singular_value = gsl_vector_get(svdS, 0);
				for (int dim_i = 0; dim_i < n_params_w_intercept; dim_i++)
				{
					// Evaluate the current singular value.
					double cur_singular_val = gsl_vector_get(svdS, dim_i);
					//fprintf(stderr, "%d: %d: %lf\n", target_var_regs->at(target_var_i)->start, i, cur_singular_val);

					double sv_tolerance_fraction = 0.001;
					double min_component_strength = 0.01;

					// If the singular value is 0, identify which variants have components on the eigenvectors with 0 singular:
					if (cur_singular_val < top_singular_value * sv_tolerance_fraction &&
						cur_singular_val < 0.001)
					{
						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "Negligible singular value: %lf (dim: %d)\n", cur_singular_val, dim_i);
						}

						// Check the entries that are on this eigenvector.
						vector<int>* colinear_var_par_i = new vector<int>();
						for (int par_i = 1; par_i < n_params_w_intercept; par_i++)
						{
							double cur_comp = gsl_matrix_get(svdV, par_i, dim_i);

							if (fabs(cur_comp) > min_component_strength)
							{
								if (__DUMP_INPUTATION_UTILS_MSGS__)
								{
									fprintf(stderr, "%d has component on 0 eigenvalue: %lf\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, cur_comp);
								}

								colinear_var_par_i->push_back(par_i - 1);

								tag_usage_flag[par_i - 1] = 0;
							}
						} // par_i loop.

						  //// Following performs variant selection to optimize the selections.
						  //int n_colinear_vars = colinear_var_par_i->size();
						  //int tag_i_2_set_unused = n_colinear_vars - 1;
						  //while (tag_i_2_set_unused >= 0)
						  //{
						  //	if (tag_usage_flag[colinear_var_par_i->at(tag_i_2_set_unused)] == 0)
						  //	{
						  //		tag_i_2_set_unused--;
						  //	}
						  //	else
						  //	{
						  //		tag_usage_flag[colinear_var_par_i->at(tag_i_2_set_unused)] = 0;
						  //		fprintf(stderr, "%d colinear tags; setting %d to be unused: %lf\n",
						  //				n_colinear_vars,
						  //				tag_var_regs_per_cur_target_var->at(colinear_var_par_i->at(tag_i_2_set_unused))->start);
						  //		break;
						  //	}
						  //}
						delete colinear_var_par_i;
					}
				} // i loop.

				n_tags_2_use = 0;
				for (int par_i = 0; par_i < n_params_w_intercept - 1; par_i++)
				{
					if (tag_usage_flag[par_i] == 1)
					{
						n_tags_2_use++;
					}
				} // par_i loop.

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					fprintf(stderr, "%s:%d (%d tags)\n", target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start, n_tags_2_use);
				}

				int n_tags_2_use_w_intercept = n_tags_2_use + 1;

				// Allocate the structures.
				X = gsl_matrix_alloc(training_sample_size, n_tags_2_use_w_intercept);
				testing_x = gsl_vector_alloc(n_tags_2_use_w_intercept);

				gsl_matrix * svdX = gsl_matrix_alloc(training_sample_size, n_tags_2_use_w_intercept);
				gsl_matrix * V = gsl_matrix_alloc(n_tags_2_use_w_intercept, n_tags_2_use_w_intercept);
				gsl_vector * S = gsl_vector_alloc(n_tags_2_use_w_intercept);
				gsl_vector * svd_work = gsl_vector_alloc(n_tags_2_use_w_intercept);

				// Target SNP genotypes vector.
				y = gsl_vector_alloc(training_sample_size);

				// Weights of each data point.
				w = gsl_vector_alloc(training_sample_size);

				// This is the coefficients matrix.
				c = gsl_vector_alloc(n_tags_2_use_w_intercept);
				cov = gsl_matrix_alloc(n_tags_2_use_w_intercept, n_tags_2_use_w_intercept);

				// Set the matrices and vectors for regression.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(X, sample_i, 0, 1.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					int matrix_dim_i = 1;
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						if (tag_usage_flag[tag_dim_i - 1] == 1)
						{
							void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
							char* cur_tag_genotypes = (char*)(reg_data[0]);

							// Set the current X matrix value.
							double tagg_ij = cur_tag_genotypes[sample_i];

							// Normalize with tag normalizer.
							tagg_ij /= tag_normalizer;

							//gsl_matrix_set(X, sample_i, tag_dim_i, tagg_ij);
							gsl_matrix_set(X, sample_i, matrix_dim_i, tagg_ij);
							matrix_dim_i++;
						}
					} // tag_dim_i loop.

					if (matrix_dim_i != n_tags_2_use_w_intercept)
					{
						fprintf(stderr, "Could not match the matrix dimension to the tags 2 use: %d, %d\n", matrix_dim_i, n_tags_2_use_w_intercept);
						exit(0);
					}

					// Set the target genotype vector entry.
					void** reg_data = (void**)(target_var_regs->at(target_var_i)->data);
					char* target_var_genotypes = (char*)(reg_data[0]);
					double targetg_ij = target_var_genotypes[sample_i];

					// Normalize with target normalizer.
					targetg_ij /= target_normalizer;

					gsl_vector_set(y, sample_i, targetg_ij);
				} // sample_i loop.				

				  // At this point, the target genotype array and the tag genotype matrix are setup.
				  // Now do the fitting.
				  // Do fitting.
				{
					// Do the fit.
					gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_tags_2_use_w_intercept);
					//gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_params_w_intercept);
					//gsl_multifit_linear_tsvd(X, y, c, cov, &chisq, &sv_rank, work); // Does not link!
					gsl_multifit_linear(X, y, c, cov, &chisq, work);
					gsl_multifit_linear_free(work);
				} // fitting block.

				  // Save the current parameters.
				  //#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

			
				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					fprintf(stderr, "Saving Parameters for %s:%d (%d-%d): %d tags used.          \r",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start, n_tags_2_use_w_intercept);

					int used_par_i = 0;
					for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
					{
						if (par_i == 0)
						{
							fprintf(stderr, "Intercept\t%.17f\n", gsl_vector_get(c, used_par_i));
							used_par_i++;
						}
						else if (tag_usage_flag[par_i - 1] == 1)
						{
							fprintf(stderr, "%d\t%.17f\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, gsl_vector_get(c, used_par_i));
							used_par_i++;
						}
						else
						{
							fprintf(stderr, "%d\tNA\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start);
						}
					} // par_i loop.
				}

				// Save the parameters file.
				char op_fp[1000];
				sprintf(op_fp, "%s/%d.params", op_dir, target_var_regs->at(target_var_i)->start);
				FILE* f_op = open_f(op_fp, "w");
				int used_par_i = 0;
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					if (par_i == 0)
					{
						fprintf(f_op, "%.17f\n", gsl_vector_get(c, used_par_i));
						fprintf(f_pooled_params, "%.17f\t", gsl_vector_get(c, used_par_i));
						used_par_i++;
					}
					else if (tag_usage_flag[par_i - 1] == 1)
					{
						fprintf(f_op, "%.17f\n", gsl_vector_get(c, used_par_i));
						fprintf(f_pooled_params, "%.17f\t", gsl_vector_get(c, used_par_i));
						used_par_i++;
					}
					else
					{
						fprintf(f_pooled_params, "NA\t");
						fprintf(f_op, "NA\n");
					}
				} // par_i loop.

				// Save the coordinates: For some weird formatting problem, we have to use the start-1 as the starting position.
				for (int tag_i = 0; tag_i < tag_var_regs_per_cur_target_var->size(); tag_i++)
				{
					fprintf(f_op, "%d\n", tag_var_regs_per_cur_target_var->at(tag_i)->start - 1);
					fprintf(f_pooled_params, "%d\t", tag_var_regs_per_cur_target_var->at(tag_i)->start - 1);
				} // par_i loop.

				fprintf(f_op, "%d\n", target_var_regs->at(target_var_i)->start - 1);
				fprintf(f_pooled_params, "%d\n", target_var_regs->at(target_var_i)->start - 1);
				fclose(f_op);

				// Test if available.
				if (testing_flag)
				{
					vector<double>* per_testing_sample_errors = new vector<double>();
					vector<double>* per_testing_sample_non_ref_errors = new vector<double>();
					for (int testing_sample_i = 0; testing_sample_i < testing_tag_sample_ids->size(); testing_sample_i++)
					{
						// Reset the testing vector x: First element of x is 2.
						gsl_vector_set(testing_x, 0, 1.0);

						int matrix_dim_i = 1;
						for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
						{
							if (tag_usage_flag[tag_dim_i - 1] == 1)
							{
								void** cur_tag_reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);

								char* cur_testing_tag_genotypes = (char*)(cur_tag_reg_data[1]);
								if (cur_testing_tag_genotypes == NULL)
								{
									fprintf(stderr, "Testing Target region is null, exiting.\n");
									exit(0);
								}

								double cur_testing_sample_geno = (double)(cur_testing_tag_genotypes[testing_sample_i]);

								// Normalize with tag normalizer.
								cur_testing_sample_geno /= tag_normalizer;

								gsl_vector_set(testing_x, matrix_dim_i, cur_testing_sample_geno);
								matrix_dim_i++;
							} // tag usage flag check.
						} // par_i loop.

						if (matrix_dim_i != n_tags_2_use_w_intercept)
						{
							fprintf(stderr, "Could not match the matrix dimensions: %d, %d\n",
								matrix_dim_i, n_tags_2_use_w_intercept);
							exit(0);
						}

						double geno_est = 0;
						double geno_err_est = 0;
						gsl_multifit_linear_est(testing_x, c, cov, &geno_est, &geno_err_est);

						// Estimate the real error.
						void** target_reg_data = (void**)(target_var_regs->at(target_var_i)->data);
						char* testing_target_var_genotypes = (char*)(target_reg_data[1]);
						if (testing_target_var_genotypes == NULL)
						{
							fprintf(stderr, "Testing Target region is null, exiting.\n");
							exit(0);
						}

						// Set the predicted genotype.
						char* pred_genoypes = (char*)(target_reg_data[2]);
						pred_genoypes[testing_sample_i] = (char)(100 * geno_est);

						double cur_testing_sample_target_geno = (double)(testing_target_var_genotypes[testing_sample_i]);

						// Normalize.
						cur_testing_sample_target_geno /= target_normalizer;

						double real_err = fabs(cur_testing_sample_target_geno - geno_est);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
								target_var_regs->at(target_var_i)->chrom,
								target_var_regs->at(target_var_i)->start,
								target_var_regs->at(target_var_i)->end,
								target_var_regs->at(target_var_i)->name,
								real_err,
								geno_err_est);
						}

						per_testing_sample_errors->push_back(real_err);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d: Sample %d: %.4f\n", target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start, testing_sample_i, real_err);
						}

						if (cur_testing_sample_target_geno != 0)
						{
							per_testing_sample_non_ref_errors->push_back(real_err);
						}
					} // testing_sample_i loop.

					  // Get the error statistics.
					double mean_err = 0;
					double std_dev_err = 0;
					get_stats(per_testing_sample_errors, mean_err, std_dev_err);

					double mean_non_ref_err = 0;
					double stddev_non_ref_err = 0;
					get_stats(per_testing_sample_non_ref_errors, mean_non_ref_err, stddev_non_ref_err);

					if (__DUMP_INPUTATION_UTILS_MSGS__)
					{
						fprintf(stderr, "%s:%d-%d (%s): Avg Error: %lf (%d); Avg Non-ref Error: %lf (%d)          \n",
							target_var_regs->at(target_var_i)->chrom,
							target_var_regs->at(target_var_i)->start,
							target_var_regs->at(target_var_i)->end,
							target_var_regs->at(target_var_i)->name,
							mean_err,
							per_testing_sample_errors->size(),
							mean_non_ref_err,
							per_testing_sample_non_ref_errors->size());
					}

					fprintf(f_per_var_stats, "%s\t%d\t%d\t%s\t%lf\t%d\t%lf\t%d\t%d\n",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						target_var_regs->at(target_var_i)->end,
						target_var_regs->at(target_var_i)->name,
						mean_err,
						per_testing_sample_errors->size(),
						mean_non_ref_err,
						per_testing_sample_non_ref_errors->size(),
						n_tags_2_use_w_intercept);

					delete per_testing_sample_non_ref_errors;
					delete per_testing_sample_errors;
				} // testing check.

				  // Free memory.
				gsl_matrix_free(X);
				gsl_vector_free(testing_x);
				gsl_matrix_free(svdX);
				gsl_matrix_free(V);
				gsl_vector_free(S);
				gsl_vector_free(svd_work);
				gsl_vector_free(y);
				gsl_vector_free(w);
				gsl_vector_free(c);
				gsl_matrix_free(cov);

				//printf("# covariance matrix:\n");
				//printf("[ %+.5e, %+.5e, %+.5e  \n",
				//	COV(0, 0), COV(0, 1), COV(0, 2));
				//printf("  %+.5e, %+.5e, %+.5e  \n",
				//	COV(1, 0), COV(1, 1), COV(1, 2));
				//printf("  %+.5e, %+.5e, %+.5e ]\n",
				//	COV(2, 0), COV(2, 1), COV(2, 2));
				//printf("# chisq = %g\n", chisq);
			} // tag snp count check.

			delete tag_var_regs_per_cur_target_var;
		} // target_var_i loop.
	} // target_chr_i loop.

	fclose(f_per_var_stats);
	fclose(f_pooled_params);
}

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

	char* op_dir)
{
	fprintf(stderr, "Building LMSE models using matrices %s, %s with sample id's in %s, %s for SNVs in range [%d-%d] using %d threads.\n",
		tag_genotype_matrix_fp, target_genotype_matrix_fp,
		tag_sample_ids_list_fp, target_sample_ids_list_fp,
		start_pos, end_pos, n_threads);

	vector<char*>* tag_sample_ids = buffer_file(tag_sample_ids_list_fp);
	fprintf(stderr, "Loading tag genotype matrix with %d training samples.\n", tag_sample_ids->size());
	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_genotype_signal_regions(tag_genotype_matrix_fp, tag_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d training samples.\n", tag_genotype_signal_regs->size(), tag_sample_ids->size());
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = restructure_annot_regions(tag_genotype_signal_regs);

	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loading target genotype matrix with %d training samples.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_genotype_signal_regs = load_variant_genotype_signal_regions(target_genotype_matrix_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d training samples.\n", target_genotype_signal_regs->size(), target_sample_ids->size());
	t_restr_annot_region_list* restr_target_genotype_signal_regs = restructure_annot_regions(target_genotype_signal_regs);

	if (target_sample_ids->size() != tag_sample_ids->size())
	{
		fprintf(stderr, "Target and tag sample id's are not the same size: %d, %d\n", target_sample_ids->size(), tag_sample_ids->size());
		exit(0);
	}

	// Make sure the sample id's are matching.
	for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
	{
		if (!t_string::compare_strings(target_sample_ids->at(i_s), tag_sample_ids->at(i_s)))
		{
			fprintf(stderr, "Target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", target_sample_ids->at(i_s), tag_sample_ids->at(i_s));
			exit(0);
		}
	} // i_s loop.

	  // Decide on haplocoded signal; tag and target normalizers.
	double tag_normalizer = -1;
	double target_normalizer = -1;
	for (int var_i = 0; var_i < tag_genotype_signal_regs->size(); var_i++)
	{
		// Get the maximum.
		void** cur_var_sig_info = (void**)(tag_genotype_signal_regs->at(var_i)->data);
		char* cur_var_sig = (char*)(cur_var_sig_info[0]);

		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			tag_normalizer = MAX(cur_var_sig[i_s], tag_normalizer);
		} // i_s loop.
	} // var_i loop.

	for (int var_i = 0; var_i < target_genotype_signal_regs->size(); var_i++)
	{
		// Get the maximum.
		void** cur_var_sig_info = (void**)(target_genotype_signal_regs->at(var_i)->data);
		char* cur_var_sig = (char*)(cur_var_sig_info[0]);

		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			target_normalizer = MAX(cur_var_sig[i_s], target_normalizer);
		} // i_s loop.
	} // var_i loop.

	bool testing_flag = false;
	vector<char*>* testing_tag_sample_ids = NULL;
	vector<char*>* testing_target_sample_ids = NULL;
	if (check_file(testing_tag_genotype_matrix_fp) &&
		check_file(testing_target_genotype_matrix_fp))
	{
		fprintf(stderr, "Testing data exists, loading and setting it.\n");
		testing_flag = true;

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_tag_sample_ids = buffer_file(testing_tag_sample_ids_list_fp);
		fprintf(stderr, "Loading testing tag genotype matrix with %d testing samples.\n", testing_tag_sample_ids->size());
		vector<t_annot_region*>* testing_tag_genotype_signal_regs = load_variant_genotype_signal_regions(testing_tag_genotype_matrix_fp, testing_tag_sample_ids);
		fprintf(stderr, "Loaded %d testing tag genotype signal regions with %d testing samples.\n", testing_tag_genotype_signal_regs->size(), testing_tag_sample_ids->size());

		vector<t_annot_region*>* tag_intersects = intersect_annot_regions(tag_genotype_signal_regs, testing_tag_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d tag intersects.\n", tag_intersects->size());
		for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
			t_annot_region* tag_reg = int_info->src_reg;
			t_annot_region* testing_tag_reg = int_info->dest_reg;

			if (t_string::compare_strings(tag_reg->name, testing_tag_reg->name))
			{
				void** cur_tag_reg_info = (void**)(tag_reg->data);
				void** cur_testing_tag_reg_info = (void**)(testing_tag_reg->data);
				cur_tag_reg_info[1] = cur_testing_tag_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(tag_intersects);

		fprintf(stderr, "Testing assignment of tag regions to all training tag regions.\n");
		for (int i_reg = 0; i_reg < tag_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_tag_reg_info = (void**)(tag_genotype_signal_regs->at(i_reg)->data);

			if (cur_tag_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing tag genotype data for %s:%d-%d (%s)\n",
					tag_genotype_signal_regs->at(i_reg)->chrom,
					tag_genotype_signal_regs->at(i_reg)->start,
					tag_genotype_signal_regs->at(i_reg)->end,
					tag_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_target_sample_ids = buffer_file(testing_target_sample_ids_list_fp);
		fprintf(stderr, "Loading testing target genotype matrix with %d testing samples.\n", testing_target_sample_ids->size());
		vector<t_annot_region*>* testing_target_genotype_signal_regs = load_variant_genotype_signal_regions(testing_target_genotype_matrix_fp, testing_target_sample_ids);
		fprintf(stderr, "Loaded %d testing target genotype signal regions with %d testing samples.\n", testing_target_genotype_signal_regs->size(), testing_target_sample_ids->size());

		vector<t_annot_region*>* target_intersects = intersect_annot_regions(target_genotype_signal_regs, testing_target_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d target intersects.\n", target_intersects->size());
		for (int i_int = 0; i_int < target_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
			t_annot_region* target_reg = int_info->src_reg;
			t_annot_region* testing_target_reg = int_info->dest_reg;

			if (t_string::compare_strings(target_reg->name, testing_target_reg->name))
			{
				void** cur_target_reg_info = (void**)(target_reg->data);
				void** cur_testing_target_reg_info = (void**)(testing_target_reg->data);
				cur_target_reg_info[1] = cur_testing_target_reg_info[0];

				void** updated_info = new void*[5];
				updated_info[0] = cur_target_reg_info[0];
				updated_info[1] = cur_target_reg_info[1];
				updated_info[2] = new char[testing_target_sample_ids->size() + 2];
				memset(updated_info[2], 0, sizeof(char)  * (testing_target_sample_ids->size() + 2));

				target_reg->data = updated_info;
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(target_intersects);

		fprintf(stderr, "Testing assignment of target regions to all training target regions.\n");
		for (int i_reg = 0; i_reg < target_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_target_reg_info = (void**)(target_genotype_signal_regs->at(i_reg)->data);

			if (cur_target_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing target genotype data for %s:%d-%d (%s)\n",
					target_genotype_signal_regs->at(i_reg)->chrom,
					target_genotype_signal_regs->at(i_reg)->start,
					target_genotype_signal_regs->at(i_reg)->end,
					target_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		if (testing_target_sample_ids->size() != testing_tag_sample_ids->size())
		{
			fprintf(stderr, "Testing target and tag sample id's are not the same size: %d, %d\n", testing_target_sample_ids->size(), testing_tag_sample_ids->size());
			exit(0);
		}

		// Make sure the sample id's are matching.
		for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Testing target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s));
				exit(0);
			}
		} // i_s loop.

		  // Go over the testing data too for updating normalizers.
		for (int var_i = 0; var_i < testing_tag_genotype_signal_regs->size(); var_i++)
		{
			// Get the maximum.
			void** cur_var_sig_info = (void**)(testing_tag_genotype_signal_regs->at(var_i)->data);
			char* cur_var_sig = (char*)(cur_var_sig_info[0]);

			for (int i_s = 0; i_s < testing_tag_sample_ids->size(); i_s++)
			{
				tag_normalizer = MAX(cur_var_sig[i_s], tag_normalizer);
			} // i_s loop.
		} // var_i loop.

		for (int var_i = 0; var_i < testing_target_genotype_signal_regs->size(); var_i++)
		{
			// Get the maximum.
			void** cur_var_sig_info = (void**)(testing_target_genotype_signal_regs->at(var_i)->data);
			char* cur_var_sig = (char*)(cur_var_sig_info[0]);

			for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
			{
				target_normalizer = MAX(cur_var_sig[i_s], target_normalizer);
			} // i_s loop.
		} // var_i loop.
	} // Testing data existence check.
	else
	{
		fprintf(stderr, "Not going to perform testing statistics.\n");
	}

	if (tag_normalizer < 2 || tag_normalizer > 3 ||
		target_normalizer < 2 || target_normalizer > 3)
	{
		fprintf(stderr, "tag or target normalizer is illegal or very improbably unlucky: %.1f, %.1f\n", tag_normalizer, target_normalizer);
		exit(0);
	}

	fprintf(stderr, "----------------------------------\n");
	if (tag_normalizer == 2)
	{
		fprintf(stderr, "Tag genotypes are unphased.\n");
	}
	else if (tag_normalizer == 3)
	{
		fprintf(stderr, "Tag genotypes are phased.\n");
	}

	if (target_normalizer == 2)
	{
		fprintf(stderr, "Target genotypes are unphased.\n");
	}
	else if (target_normalizer == 3)
	{
		fprintf(stderr, "Target genotypes are phased.\n");
	}
	fprintf(stderr, "----------------------------------\n");

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// This is the number of data points in training.
	int training_sample_size = target_sample_ids->size();

	vector<t_ansi_thread*>* impute_threads = new vector<t_ansi_thread*>();
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		t_LMSE_imputation_thread_info* cur_thread_info = new t_LMSE_imputation_thread_info();
		cur_thread_info->which = i_thread;
		cur_thread_info->outof = n_threads;

		cur_thread_info->target_normalizer = target_normalizer;
		cur_thread_info->tag_normalizer = tag_normalizer;

		//cur_thread_info->testing_tag_genotype_signal_regs = testing_tag_genotype_signal_regs;
		cur_thread_info->tag_genotype_signal_regs = tag_genotype_signal_regs;
		//cur_thread_info->testing_target_genotype_signal_regs = testing_target_genotype_signal_regs;
		cur_thread_info->target_genotype_signal_regs = target_genotype_signal_regs;

		cur_thread_info->restr_target_genotype_signal_regs = restr_target_genotype_signal_regs;
		cur_thread_info->restr_tag_genotype_signal_regs = restr_tag_genotype_signal_regs;

		cur_thread_info->testing_flag = testing_flag;
		cur_thread_info->start_pos = start_pos;
		cur_thread_info->end_pos = end_pos;
		cur_thread_info->op_dir = op_dir;

		cur_thread_info->tag_sample_ids = tag_sample_ids;
		cur_thread_info->target_sample_ids = target_sample_ids;

		cur_thread_info->testing_tag_sample_ids = testing_tag_sample_ids;
		cur_thread_info->testing_target_sample_ids = testing_target_sample_ids;

		cur_thread_info->n_tags_vars_per_side = n_tags_vars_per_side;

		t_ansi_thread* new_thread = new t_ansi_thread(LMSE_imputation_thread_callback, cur_thread_info);

		new_thread->run_thread();
		impute_threads->push_back(new_thread);
	} // thread_i loop.

	// Wait for the threads.
	for (int i_thread = 0; i_thread < impute_threads->size(); i_thread++)
	{
		impute_threads->at(i_thread)->wait_thread();
	} // i_thread loop.

	// Save results.
	char geno_op_fp[1000];
	sprintf(geno_op_fp, "%s/LMSE_imputed_genotypes.sig", op_dir);

	// Following sets up n data points. We set this up over the samples.
	vector<t_annot_region*>* target_regs_2_dump = new vector<t_annot_region*>();
	for (int target_var_i = 0; target_var_i < target_genotype_signal_regs->size(); target_var_i++)
	{
		void** target_reg_info = (void**)(target_genotype_signal_regs->at(target_var_i)->data);

		if (target_reg_info[1] == NULL || target_reg_info[2] == NULL)
		{

		}
		else
		{
			t_annot_region* cur_reg = duplicate_region(target_genotype_signal_regs->at(target_var_i));
			void** cur_reg_info = new void*[5];
			cur_reg_info[0] = (char*)(target_reg_info[2]);
			cur_reg->data = cur_reg_info;

			target_regs_2_dump->push_back(cur_reg);
		}
	} // target_var_i loop.

	fprintf(stderr, "Saving %d regions to %s\n", target_regs_2_dump->size(), geno_op_fp);
	dump_geno_sig_regs_plain(target_regs_2_dump, testing_tag_sample_ids, false, geno_op_fp);
}

// http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_35.html
void train_vicinity_based_LMSE_imputation_model_ALL_FEATURES(int n_tags_vars_per_side,
												char* tag_genotype_matrix_fp,
												char* tag_sample_ids_list_fp,
												char* target_genotype_matrix_fp,
												char* target_sample_ids_list_fp,

												// The testing data.
												char* testing_tag_genotype_matrix_fp,
												char* testing_tag_sample_ids_list_fp,
												char* testing_target_genotype_matrix_fp,
												char* testing_target_sample_ids_list_fp,

												char* op_dir)
{
	fprintf(stderr, "Building LMSE models using matrices %s, %s with sample id's in %s, %s\n", 
			tag_genotype_matrix_fp, target_genotype_matrix_fp, 
			tag_sample_ids_list_fp, target_sample_ids_list_fp);

	vector<char*>* tag_sample_ids = buffer_file(tag_sample_ids_list_fp);
	fprintf(stderr, "Loading tag genotype matrix with %d testing samples.\n", tag_sample_ids->size());
	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_genotype_signal_regions(tag_genotype_matrix_fp, tag_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d testing samples.\n", tag_genotype_signal_regs->size(), tag_sample_ids->size());
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = restructure_annot_regions(tag_genotype_signal_regs);

	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loading target genotype matrix with %d testing samples.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_genotype_signal_regs = load_variant_genotype_signal_regions(target_genotype_matrix_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d testing samples.\n", target_genotype_signal_regs->size(), target_sample_ids->size());
	t_restr_annot_region_list* restr_target_genotype_signal_regs = restructure_annot_regions(target_genotype_signal_regs);

	if (target_sample_ids->size() != tag_sample_ids->size())
	{
		fprintf(stderr, "Target and tag sample id's are not the same size: %d, %d\n", target_sample_ids->size(), tag_sample_ids->size());
		exit(0);
	}

	// Make sure the sample id's are matching.
	for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
	{
		if (!t_string::compare_strings(target_sample_ids->at(i_s), tag_sample_ids->at(i_s)))
		{
			fprintf(stderr, "Target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", target_sample_ids->at(i_s), tag_sample_ids->at(i_s));
			exit(0);
		}
	} // i_s loop.

	bool testing_flag = false;
	vector<char*>* testing_tag_sample_ids = NULL;
	vector<char*>* testing_target_sample_ids = NULL;
	if (check_file(testing_tag_genotype_matrix_fp) &&
		check_file(testing_target_genotype_matrix_fp))
	{
		fprintf(stderr, "Testing data exists, loading and setting it.\n");
		testing_flag = true;

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_tag_sample_ids = buffer_file(testing_tag_sample_ids_list_fp);
		fprintf(stderr, "Loading testing tag genotype matrix with %d testing samples.\n", testing_tag_sample_ids->size());
		vector<t_annot_region*>* testing_tag_genotype_signal_regs = load_variant_genotype_signal_regions(testing_tag_genotype_matrix_fp, testing_tag_sample_ids);
		fprintf(stderr, "Loaded %d testing tag genotype signal regions with %d testing samples.\n", testing_tag_genotype_signal_regs->size(), testing_tag_sample_ids->size());

		vector<t_annot_region*>* tag_intersects = intersect_annot_regions(tag_genotype_signal_regs, testing_tag_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d tag intersects.\n", tag_intersects->size());
		for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
			t_annot_region* tag_reg = int_info->src_reg;
			t_annot_region* testing_tag_reg = int_info->dest_reg;

			if (t_string::compare_strings(tag_reg->name, testing_tag_reg->name))
			{
				void** cur_tag_reg_info = (void**)(tag_reg->data);
				void** cur_testing_tag_reg_info = (void**)(testing_tag_reg->data);
				cur_tag_reg_info[1] = cur_testing_tag_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(tag_intersects);

		fprintf(stderr, "Testing assignment of tag regions to all training tag regions.\n");
		for (int i_reg = 0; i_reg < tag_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_tag_reg_info = (void**)(tag_genotype_signal_regs->at(i_reg)->data);

			if (cur_tag_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing tag genotype data for %s:%d-%d (%s)\n", 
					tag_genotype_signal_regs->at(i_reg)->chrom, 
					tag_genotype_signal_regs->at(i_reg)->start, 
					tag_genotype_signal_regs->at(i_reg)->end,
					tag_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_target_sample_ids = buffer_file(testing_target_sample_ids_list_fp);
		fprintf(stderr, "Loading testing target genotype matrix with %d testing samples.\n", testing_target_sample_ids->size());
		vector<t_annot_region*>* testing_target_genotype_signal_regs = load_variant_genotype_signal_regions(testing_target_genotype_matrix_fp, testing_target_sample_ids);
		fprintf(stderr, "Loaded %d testing target genotype signal regions with %d testing samples.\n", testing_target_genotype_signal_regs->size(), testing_target_sample_ids->size());

		vector<t_annot_region*>* target_intersects = intersect_annot_regions(target_genotype_signal_regs, testing_target_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d target intersects.\n", target_intersects->size());
		for (int i_int = 0; i_int < target_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
			t_annot_region* target_reg = int_info->src_reg;
			t_annot_region* testing_target_reg = int_info->dest_reg;

			if (t_string::compare_strings(target_reg->name, testing_target_reg->name))
			{
				void** cur_target_reg_info = (void**)(target_reg->data);
				void** cur_testing_target_reg_info = (void**)(testing_target_reg->data);
				cur_target_reg_info[1] = cur_testing_target_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(target_intersects);

		fprintf(stderr, "Testing assignment of target regions to all training target regions.\n");
		for (int i_reg = 0; i_reg < target_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_target_reg_info = (void**)(target_genotype_signal_regs->at(i_reg)->data);

			if (cur_target_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing target genotype data for %s:%d-%d (%s)\n",
					target_genotype_signal_regs->at(i_reg)->chrom,
					target_genotype_signal_regs->at(i_reg)->start,
					target_genotype_signal_regs->at(i_reg)->end,
					target_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		if (testing_target_sample_ids->size() != testing_tag_sample_ids->size())
		{
			fprintf(stderr, "Testing target and tag sample id's are not the same size: %d, %d\n", testing_target_sample_ids->size(), testing_tag_sample_ids->size());
			exit(0);
		}

		// Make sure the sample id's are matching.
		for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Testing target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s));
				exit(0);
			}
		} // i_s loop.

	} // Testing data existence check.

	// Allocate the fitting data.
	//int i, n;
	double xi, yi, ei, chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	gsl_vector *testing_x;

	// This is the number of data points in training.
	//n = atoi(argv[1]);
	//int n = training_sample_ids->size();
	int training_sample_size = target_sample_ids->size();

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// Tag SNP genotypes
	//X = gsl_matrix_alloc(n, 3);
	X = gsl_matrix_alloc(training_sample_size, n_params_w_intercept);
	testing_x = gsl_vector_alloc(n_params_w_intercept);

	// Target SNP genotypes vector.
	y = gsl_vector_alloc(training_sample_size);
	
	// Weights of each data point.
	w = gsl_vector_alloc(training_sample_size);

	// This is the coefficients matrix.
	c = gsl_vector_alloc(n_params_w_intercept);
	cov = gsl_matrix_alloc(n_params_w_intercept, n_params_w_intercept);

	// Process all the target variants.
	for (int target_chr_i = 0;
		target_chr_i < restr_target_genotype_signal_regs->chr_ids->size();
		target_chr_i++)
	{
		fprintf(stderr, "Building the models on chromosome %s\n", restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		// Get the chromosome index for testring variant regions.
		int tag_chr_i = t_string::get_i_str(restr_tag_genotype_signal_regs->chr_ids, restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		vector<t_annot_region*>* target_var_regs = restr_target_genotype_signal_regs->regions_per_chrom[target_chr_i];
		vector<t_annot_region*>* tag_var_regs = restr_tag_genotype_signal_regs->regions_per_chrom[tag_chr_i];

		// Following sets up n data points. We set this up over the samples.
		for (int target_var_i = 0; target_var_i < target_var_regs->size(); target_var_i++)
		{
			int closest_tag_var_left_i = -1;
			for (int tag_var_i = 0; (tag_var_i+1)< tag_var_regs->size(); tag_var_i++)
			{
				if (tag_var_regs->at(tag_var_i)->start < target_var_regs->at(target_var_i)->start &&
					tag_var_regs->at(tag_var_i+1)->start > target_var_regs->at(target_var_i)->start)
				{
					closest_tag_var_left_i = tag_var_i;
					break;
				}
			} // tag_var_i loop.

			if (closest_tag_var_left_i == -1)
			{
				continue;
			}

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Closest tag to target %d: (%d, %d)\n", target_var_regs->at(target_var_i)->start,
						tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);
			}

			// Add the tag variant to the right.
			vector<t_annot_region*>* tag_var_regs_per_cur_target_var = new vector<t_annot_region*>();
			for (int tag_var_i = MAX(0, closest_tag_var_left_i - n_tags_vars_per_side + 1);
				tag_var_i <= closest_tag_var_left_i;
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.
					
			for (int tag_var_i = closest_tag_var_left_i + 1;
				tag_var_i <= (closest_tag_var_left_i + n_tags_vars_per_side) && tag_var_i < tag_var_regs->size();
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			// Check the number of tag SNVs around the target SNVs, make sure they satisfy the requested number.
			if (tag_var_regs_per_cur_target_var->size() != 2 * n_tags_vars_per_side)
			{
				fprintf(stderr, "Could not find the correct number of tag variants for %s:%d (%d), will skip.\n",
						target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start,
						tag_var_regs_per_cur_target_var->size());
			}
			else
			{
				// Check the tag variant genotype correlation among samples.
				double cur_n_per_geno[3];	
				double* dim_i_genotypes = new double[testing_tag_sample_ids->size() + 2];
				double* dim_j_genotypes = new double[testing_tag_sample_ids->size() + 2];
				int* tag_usage_flag = new int[2 * n_tags_vars_per_side + 3];
				memset(tag_usage_flag, 0, sizeof(int) * (2 * n_tags_vars_per_side + 3));
				for (int tag_dim_i = 0; tag_dim_i < 2 * n_tags_vars_per_side; tag_dim_i++)
				{
					memset(cur_n_per_geno, 0, sizeof(double) * 3);

					for (int sample_i = 0;
						sample_i < tag_sample_ids->size();
						sample_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						dim_i_genotypes[sample_i] = cur_tag_genotypes[sample_i];

						// Set the 
						cur_n_per_geno[cur_tag_genotypes[sample_i]]++;
					} // sample_i loop.

					  // Analyze the distribution.
					if (cur_n_per_geno[0] == tag_sample_ids->size() ||
						cur_n_per_geno[1] == tag_sample_ids->size() ||
						cur_n_per_geno[2] == tag_sample_ids->size())
					{
						fprintf(stderr, "Dimension %d provides no info.\n",
							tag_var_regs_per_cur_target_var->at(tag_dim_i)->chrom, tag_var_regs_per_cur_target_var->at(tag_dim_i)->start);
					}

					double min_distance_dim_i = 1000000;
					for (int tag_dim_j = 0; tag_dim_j < tag_dim_i; tag_dim_j++)
					{
						for (int sample_i = 0;
							sample_i < tag_sample_ids->size();
							sample_i++)
						{
							void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_j)->data);
							char* cur_tag_genotypes = (char*)(reg_data[0]);

							dim_j_genotypes[sample_i] = cur_tag_genotypes[sample_i];
						} // sample_i loop.

						// Compute the distance.
						double cur_dist = 0;
						for (int sample_i = 0; sample_i < tag_sample_ids->size(); sample_i++)
						{
							cur_dist += (dim_i_genotypes[sample_i] != dim_j_genotypes[sample_i]);
						} // samepl_i loop.

						// Update the min distance.
						min_distance_dim_i = MIN(min_distance_dim_i, cur_dist);
					} // tag_dim_j;

					// Do greedy selection.
					if (min_distance_dim_i > 0)
					{
						tag_usage_flag[tag_dim_i] = 1;
					}
				} // tag_dim_i;

				delete[] dim_i_genotypes;
				delete[] dim_j_genotypes;

				// Set the matrices and vectors for regression.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(X, sample_i, 0, 2.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						// Set the 
						double tagg_ij = cur_tag_genotypes[sample_i];
						gsl_matrix_set(X, sample_i, tag_dim_i, tagg_ij);
					} // tag_dim_i loop.

					// Set the target genotype vector entry.
					void** reg_data = (void**)(target_var_regs->at(target_var_i)->data);
					char* target_var_genotypes = (char*)(reg_data[0]);
					double targetg_ij = target_var_genotypes[sample_i];
					gsl_vector_set(y, sample_i, targetg_ij);	
				} // sample_i loop.

				// At this point, the target genotype array and the tag genotype matrix are setup.
				// Now do the fitting.
				// Do fitting.
				{
					size_t sv_rank;

					gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_params_w_intercept);
					//gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
					//gsl_multifit_linear_tsvd(X, y, c, cov, &chisq, &sv_rank, work);
					gsl_multifit_linear(X, y, c, cov, &chisq, work);
					gsl_multifit_linear_free(work);
				} // fitting block.

				// Save the current parameters.
//#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
				
				fprintf(stderr, "Saving Parameters for %s:%d (%d-%d)          \r", 
						target_var_regs->at(target_var_i)->chrom, 
						target_var_regs->at(target_var_i)->start,
						tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
					{
						if (par_i == 0)
						{
							fprintf(stderr, "Intercept\t%.17f\n", gsl_vector_get(c, par_i));
						}
						else
						{
							fprintf(stderr, "%d\t%.17f\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, gsl_vector_get(c, par_i));
						}

						//printf("# best fit: Y = %g + %g X + %g X^2\n",
						//	C(0), C(1), C(2));
					} // par_i loop.
				}

				// Save the parameters file.
				char op_fp[1000];
				sprintf(op_fp, "%s/%d.params", op_dir, target_var_regs->at(target_var_i)->start);
				FILE* f_op = open_f(op_fp, "w");
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					fprintf(f_op, "%.17f\n", gsl_vector_get(c, par_i));
				} // par_i loop.

				for (int tag_i = 0; tag_i < tag_var_regs_per_cur_target_var->size(); tag_i++)
				{
					fprintf(f_op, "%d\n", tag_var_regs_per_cur_target_var->at(tag_i)->start);
				} // par_i loop.

				fprintf(f_op, "%d\n", target_var_regs->at(target_var_i)->start);
				fclose(f_op);

				// Test if available.
				if (testing_flag)
				{
					vector<double>* per_testing_sample_errors = new vector<double>();
					for (int testing_sample_i = 0; testing_sample_i < testing_tag_sample_ids->size(); testing_sample_i++)
					{
						// Reset the testing vector x: First element of x is 2.
						gsl_vector_set(testing_x, 0, 2);

						for (int reg_i = 0; reg_i < tag_var_regs_per_cur_target_var->size(); reg_i++)
						{
							void** cur_tag_reg_data = (void**)(tag_var_regs_per_cur_target_var->at(reg_i)->data);
							char* cur_testing_tag_genotypes = (char*)(cur_tag_reg_data[1]);
							if (cur_testing_tag_genotypes == NULL)
							{
								fprintf(stderr, "Testing Target region is null, exiting.\n");
								exit(0);
							}

							double cur_testing_sample_geno = (double)(cur_testing_tag_genotypes[testing_sample_i]);
							gsl_vector_set(testing_x, reg_i+1, cur_testing_sample_geno);	
						} // par_i loop.

						double geno_est = 0;
						double geno_err_est = 0;
						gsl_multifit_linear_est(testing_x, c, cov, &geno_est, &geno_err_est);

						// Estimate the real error.
						void** target_reg_data = (void**)(target_var_regs->at(target_var_i)->data);
						char* testing_target_var_genotypes = (char*)(target_reg_data[1]);
						if (testing_target_var_genotypes == NULL)
						{
							fprintf(stderr, "Testing Target region is null, exiting.\n");
							exit(0);
						}

						double cur_testing_sample_target_geno = (double)(testing_target_var_genotypes[testing_sample_i]);
						double real_err = fabs(cur_testing_sample_target_geno - geno_est);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
								target_var_regs->at(target_var_i)->chrom,
								target_var_regs->at(target_var_i)->start,
								target_var_regs->at(target_var_i)->end,
								target_var_regs->at(target_var_i)->name,
								real_err,
								geno_err_est);
						}

						per_testing_sample_errors->push_back(real_err);
					} // testing_sample_i loop.

					// Get the error statistics.
					double mean_err = 0;
					double std_dev_err = 0;
					get_stats(per_testing_sample_errors, mean_err, std_dev_err);

					fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						target_var_regs->at(target_var_i)->end,
						target_var_regs->at(target_var_i)->name,
						mean_err,
						std_dev_err);
				} // testing check.

				//printf("# covariance matrix:\n");
				//printf("[ %+.5e, %+.5e, %+.5e  \n",
				//	COV(0, 0), COV(0, 1), COV(0, 2));
				//printf("  %+.5e, %+.5e, %+.5e  \n",
				//	COV(1, 0), COV(1, 1), COV(1, 2));
				//printf("  %+.5e, %+.5e, %+.5e ]\n",
				//	COV(2, 0), COV(2, 1), COV(2, 2));
				//printf("# chisq = %g\n", chisq);
			} // tag snp count check.
		} // target_var_i loop.
	} // target_chr_i loop.

	// Free memory.
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	
//#endif 
}

