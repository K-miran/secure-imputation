# Download and pare 1kG data, extract the tag and target SNVs.

HEMPUTE_BIN=/home/aharmanci1/dir/codebase/genomics-codebase/Genomics/CryptoSeq/HEmpute_Deployment/bin/HEmpute_Train

KG_PARSE=1
if [ ${KG_PARSE} == 1 ]
then
	echo "Downloading data, parsing sample ids."
	#wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	#gzip -cd ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk 'BEGIN{FS="\t"}{if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i};exit}}' > 1kg_sample_ids.list

	# Extract the common SNVs.
	echo "Parsing min MAF regions"
	min_MAF=0.005
	#gzip -cd ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk -v min_MAF=${min_MAF} {'curAF=0;split($8, arr, ";");for(i=1; i <= length(arr); i++){split(arr[i], arr2, "=");if(arr2[1]=="AF"){curAF=arr2[2];break;}};curMAF=curAF;if(curAF>0.5){curMAF=1-curAF;};if(curMAF>min_MAF && $1==22 && length($4)==1 && length($5)==1)print $1"\t"$2-1"\t"$2"\t"$3"_"$4"_"$5"_"curAF"\t.\t+"'} > snv_region_${min_MAF}.bed

	# Extract the matrix.
	echo "Extracting MAF > ${min_MAF} variants on snv_region_${min_MAF}.bed"
	gzip -cd ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk {'if($1==22 && length($4)==1 && length($5)==1)print $0'} | $HEMPUTE_BIN -extract_genotype_signals_per_VCF stdin 1kg_sample_ids.list snv_region_${min_MAF}.bed 0 0 1 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz
fi

