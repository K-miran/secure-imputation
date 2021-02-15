# Set up the path to HEmpute-Train below.
HEMPUTE_BIN=HEmpute_Train

# Download the Illumina Duo V3 manifest.
PARSE_ILLUMINADUOV3_MANIFEST=0
if [ ${PARSE_ILLUMINADUOV3_MANIFEST} == 1 ]
then
	wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/c1859532-43c5-4df2-a1b4-c9f0374476d2/Human1M-Duov3_H.csv

	head -n 8 Human1M-Duov3_H.csv | tail -n 1 | awk {'n_ents=split($0, arr, ",");for(i=1;i<=n_ents;i++){print i":"arr[i]}'} > col_ids.list

	echo "Parsing SNV loci from manifest"
	grep -v cnv Human1M-Duov3_H.csv | awk 'BEGIN{FS=","}{if(NR>18 && $1!="CHR" && $10!=0 && $11>1)print $10"\t"$11-1"\t"$11"\t"$1" "$2" "$17"\t.\t"$22}' > Human1M-Duov3_H.csv_SNV_loci.bed
fi

# Download and pare 1kG data, extract the tag and target SNVs.
PARSE_TAG_TARGET_SNPs=0
if [ ${PARSE_TAG_TARGET_SNPs} == 1 ]
then
	# Extract the common SNVs.
	min_MAF=0.005

	# Extract the tag and target SNP matrices.
	echo "Extracting the genotypes for the tag SNVs"
	$HEMPUTE_BIN -extract_genotype_signals_per_region_list ../1kG/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz ../1kG/1kg_sample_ids.list Human1M-Duov3_H.csv_SNV_loci.bed tag_snvs_haplocoded.matbed.gz	

	echo "Saving target variants"
	$HEMPUTE_BIN -dump_plain_geno_signal_regions ../1kG/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz ../1kG/1kg_sample_ids.list 0 temp.txt
	annot_region_tools -exclude temp.txt Human1M-Duov3_H.csv_SNV_loci.bed n
	rm -f temp.txt
	mv excluded.bed target_snvs.matbed.gz.txt
	awk 'BEGIN{FS="\t"}{n_toks=split($4, arr, "_");curAF=arr[n_toks];curMAF=curAF;if(curMAF>.5){curMAF=1-curMAF;}if(curMAF>0.05)print $0}' target_snvs.matbed.gz.txt > temp.txt
	mv temp.txt target_snvs.matbed.gz.txt

	echo "Binarizing target variants."
	$HEMPUTE_BIN -binarize_genotype_signals_BED target_snvs.matbed.gz.txt ../1kG/1kg_sample_ids.list target_snvs_haplocoded.matbed.gz
fi

BUILD_TAG_TARGET_TRAINING_TESTING_DATA=0
if [ ${BUILD_TAG_TARGET_TRAINING_TESTING_DATA} == 1 ]
then
	echo "Separating data"

	GET_NEW_TRAINING_TESTING=1
	if [ ${GET_NEW_TRAINING_TESTING} == 1 ]
	then
		n_training_samples=1500
		shuf ../1kG/1kg_sample_ids.list | head -n ${n_training_samples} > training_sample_ids.list
		grep -v -f training_sample_ids.list ../1kG/1kg_sample_ids.list > testing_sample_ids.list
	fi

	$HEMPUTE_BIN -extract_genotype_signals_per_subsample_list tag_snvs_haplocoded.matbed.gz ../1kG/1kg_sample_ids.list training_sample_ids.list tag_snvs_haplocoded.matbed.gz_training.matbed.gz
	$HEMPUTE_BIN -convert_haplocoded_2_genocoded tag_snvs_haplocoded.matbed.gz_training.matbed.gz training_sample_ids.list tag_snvs_genocoded.matbed.gz_training.matbed.gz	
	$HEMPUTE_BIN -dump_plain_geno_signal_regions tag_snvs_genocoded.matbed.gz_training.matbed.gz training_sample_ids.list 0 tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt
	sort -n -k2,2 tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt > sorted_tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt

        $HEMPUTE_BIN -extract_genotype_signals_per_subsample_list tag_snvs_haplocoded.matbed.gz ../1kG/1kg_sample_ids.list testing_sample_ids.list tag_snvs_haplocoded.matbed.gz_testing.matbed.gz
        $HEMPUTE_BIN -convert_haplocoded_2_genocoded tag_snvs_haplocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list tag_snvs_genocoded.matbed.gz_testing.matbed.gz
	$HEMPUTE_BIN -dump_plain_geno_signal_regions tag_snvs_genocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list 0 tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt
	sort -n -k2,2 tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt > sorted_tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt

        $HEMPUTE_BIN -extract_genotype_signals_per_subsample_list target_snvs_haplocoded.matbed.gz ../1kG/1kg_sample_ids.list training_sample_ids.list target_snvs_haplocoded.matbed.gz_training.matbed.gz
        $HEMPUTE_BIN -convert_haplocoded_2_genocoded target_snvs_haplocoded.matbed.gz_training.matbed.gz training_sample_ids.list target_snvs_genocoded.matbed.gz_training.matbed.gz
	$HEMPUTE_BIN -dump_plain_geno_signal_regions target_snvs_genocoded.matbed.gz_training.matbed.gz training_sample_ids.list 0 target_snvs_genocoded.matbed.gz_training.matbed.gz.txt
	sort -n -k2,2 target_snvs_genocoded.matbed.gz_training.matbed.gz.txt > sorted_target_snvs_genocoded.matbed.gz_training.matbed.gz.txt

        $HEMPUTE_BIN -extract_genotype_signals_per_subsample_list target_snvs_haplocoded.matbed.gz ../1kG/1kg_sample_ids.list testing_sample_ids.list target_snvs_haplocoded.matbed.gz_testing.matbed.gz
        $HEMPUTE_BIN -convert_haplocoded_2_genocoded target_snvs_haplocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list target_snvs_genocoded.matbed.gz_testing.matbed.gz
	$HEMPUTE_BIN -dump_plain_geno_signal_regions target_snvs_genocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list 0 target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt
	sort -n -k2,2 target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt > sorted_target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt
fi

EXTRACT_SUPERPOP_DATA=0
if [ ${EXTRACT_SUPERPOP_DATA} == 1 ]
then
	POP_DATA_DIR=POP_DATA
	mkdir ${POP_DATA_DIR}
	super_pops=('EUR' 'AFR' 'AMR')
	for sup_pop in ${super_pops[@]}
	do
		echo "Extracting data for "$sup_pop
		awk -v sup_pop=$sup_pop 'BEGIN{FS="\t"}{if($3==sup_pop)print $1}' ../1kG/super_pop_info.txt > temp_cur_super_pop_pops.list
		grep -f temp_cur_super_pop_pops.list ../1kG/population_info.list | awk 'BEGIN{FS="\t"}{print $1}' > temp_cur_super_pop_sample_ids.list

		# Clean and create the directory.
		rm -f -r ${POP_DATA_DIR}/$sup_pop
		mkdir ${POP_DATA_DIR}/$sup_pop

		grep -f temp_cur_super_pop_sample_ids.list training_sample_ids.list > ${POP_DATA_DIR}/$sup_pop/${sup_pop}_training_sample_ids.list
		grep -f temp_cur_super_pop_sample_ids.list testing_sample_ids.list > ${POP_DATA_DIR}/$sup_pop/${sup_pop}_testing_sample_ids.list

		cp ${POP_DATA_DIR}/$sup_pop/${sup_pop}_training_sample_ids.list .
		cp ${POP_DATA_DIR}/$sup_pop/${sup_pop}_testing_sample_ids.list .

		n_train=`awk 'END{print NR}' ${sup_pop}_training_sample_ids.list`
		n_test=`awk 'END{print NR}' ${sup_pop}_testing_sample_ids.list`
		echo "${n_train} training ,${n_test} testing samples."

		$HEMPUTE_BIN -extract_genotype_signals_per_subsample_list tag_snvs_genocoded.matbed.gz_training.matbed.gz training_sample_ids.list ${sup_pop}_training_sample_ids.list ${POP_DATA_DIR}/${sup_pop}/training_tag.matbed.gz
		$HEMPUTE_BIN -dump_plain_geno_signal_regions ${POP_DATA_DIR}/${sup_pop}/training_tag.matbed.gz ${sup_pop}_training_sample_ids.list 0 ${POP_DATA_DIR}/${sup_pop}/tag_snvs_genocoded.matbed.gz_training.matbed.gz_${sup_pop}.txt
		sort -n -k2,2 ${POP_DATA_DIR}/${sup_pop}/tag_snvs_genocoded.matbed.gz_training.matbed.gz_${sup_pop}.txt > ${POP_DATA_DIR}/${sup_pop}/sorted_tag_snvs_genocoded.matbed.gz_training.matbed.gz_${sup_pop}.txt

		$HEMPUTE_BIN -extract_genotype_signals_per_subsample_list tag_snvs_genocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list ${sup_pop}_testing_sample_ids.list ${POP_DATA_DIR}/${sup_pop}/testing_tag.matbed.gz
                $HEMPUTE_BIN -dump_plain_geno_signal_regions ${POP_DATA_DIR}/${sup_pop}/testing_tag.matbed.gz ${sup_pop}_testing_sample_ids.list 0 ${POP_DATA_DIR}/${sup_pop}/tag_snvs_genocoded.matbed.gz_testing.matbed.gz_${sup_pop}.txt
                sort -n -k2,2 ${POP_DATA_DIR}/${sup_pop}/tag_snvs_genocoded.matbed.gz_testing.matbed.gz_${sup_pop}.txt > ${POP_DATA_DIR}/${sup_pop}/sorted_tag_snvs_genocoded.matbed.gz_testing.matbed.gz_${sup_pop}.txt

		$HEMPUTE_BIN -extract_genotype_signals_per_subsample_list target_snvs_genocoded.matbed.gz_training.matbed.gz training_sample_ids.list ${sup_pop}_training_sample_ids.list ${POP_DATA_DIR}/${sup_pop}/training_target.matbed.gz
		$HEMPUTE_BIN -dump_plain_geno_signal_regions ${POP_DATA_DIR}/${sup_pop}/training_target.matbed.gz ${sup_pop}_training_sample_ids.list 0 ${POP_DATA_DIR}/${sup_pop}/target_snvs_genocoded.matbed.gz_training.matbed.gz_${sup_pop}.txt
		sort -n -k2,2 ${POP_DATA_DIR}/${sup_pop}/target_snvs_genocoded.matbed.gz_training.matbed.gz_${sup_pop}.txt > ${POP_DATA_DIR}/${sup_pop}/sorted_target_snvs_genocoded.matbed.gz_training.matbed.gz_${sup_pop}.txt

		$HEMPUTE_BIN -extract_genotype_signals_per_subsample_list target_snvs_genocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list ${sup_pop}_testing_sample_ids.list ${POP_DATA_DIR}/${sup_pop}/testing_target.matbed.gz
		$HEMPUTE_BIN -dump_plain_geno_signal_regions ${POP_DATA_DIR}/${sup_pop}/testing_target.matbed.gz ${sup_pop}_testing_sample_ids.list 0 ${POP_DATA_DIR}/${sup_pop}/target_snvs_genocoded.matbed.gz_testing.matbed.gz_${sup_pop}.txt
		sort -n -k2,2 ${POP_DATA_DIR}/${sup_pop}/target_snvs_genocoded.matbed.gz_testing.matbed.gz_${sup_pop}.txt > ${POP_DATA_DIR}/${sup_pop}/sorted_target_snvs_genocoded.matbed.gz_testing.matbed.gz_${sup_pop}.txt

	done
fi

LMSE_MODEL_RUN=1
if [ $LMSE_MODEL_RUN == 1 ]
then
        impute_reg_start=0
        impute_reg_end=52000000
        n_vic_snvs=32

	n_threads=16

	RUN_FULL_DATA=1
	if [ $RUN_FULL_DATA == 1 ]
	then
        	rm -f temp_train_cmds.sh

	        testing_target="$PWD/sorted_target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt $PWD/testing_sample_ids.list"
	        testing_tag="$PWD/sorted_tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt $PWD/testing_sample_ids.list"
	        #testing_target="NONE NONE"
	        #testing_tag="NONE NONE"

	        cur_n_vic=${n_vic_snvs}
	        mkdir LMSE_models_${cur_n_vic}
	        start_pos=-1
	        end_pos=-1
	        echo "/usr/bin/time -o timing_memory.log -f %e\""\\t"\"%M $HEMPUTE_BIN -train_vicinity_based_LMSE_imputation_model_multithreaded $start_pos $end_pos ${cur_n_vic} \
$PWD/sorted_tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt $PWD/training_sample_ids.list \
$PWD/sorted_target_snvs_genocoded.matbed.gz_training.matbed.gz.txt $PWD/training_sample_ids.list \
$testing_tag \
$testing_target \
$PWD/LMSE_models_${cur_n_vic} \
$n_threads ;\
tar -cvjf $PWD/LMSE_models_${cur_n_vic}.tar.bz2 $PWD/LMSE_models_${cur_n_vic}" >> temp_train_cmds.sh
        	chmod 755 temp_train_cmds.sh
	fi
fi

