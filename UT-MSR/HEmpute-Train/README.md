
<html>
<font face="arial">
<h1>HEmpute-Train: Plantext linear model training for Secure Imputation</h1>

HEmpute-Train builds the linear models for secure imputation pipelines of UTHealth-CKKS and UTHealth-BFV. HEmpute-Train also provides options for processing VCF data to build the text files that the secure imputation methods take as input.<br>

The code takes as input:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
- Tag variant genotypes for the reference population panel, <br>
- Target variant genotypes for the reference population panel, <br>
- Optional: the known tag variant genotypes for the study population panel, <br>
- Optional: the known genotypes of the target variant genotypes for performing accuracy benchmarks. <br>
</div><br>
and outputs:<br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
- The paramter set for the mean squared-error optimized linear model weights for each target variant -- one file each variant. <br>
- Optional: The imputed variant genotypes if the study variant genotypes are supplied. <br>
</div>

The study tag and target variants can be used to test the models at the server side to test the accuracy of the plaintext imputation models.<br>

The output is formatted as an extended BED file with multiple columns. Please refer below for the specification of output file format.

<h2>Download and Installation</h2>
You can download HEmpute-Train code by clicking on the green "Clone or Download" button and downloading zip or by checking out via git. After that, navigate to the checkout (or download) directory. If you downloaded the zip file of the source, unzip it using "unzip master.zip".<br><br>

You need to have g++, gzip, and the GSL libraries installed for building HEmpute-Train. These are installed in most Unix distributions but if they are not installed, type:
```
sudo yum -y install gsl gsl-devel gcc-c++ gzip
```

Now HEmpute-Train can be built using:<br>
```
make clean
make
```
The executable is located under directory <font face="courier">bin/</font>. 

<h2>Example Run</h2>
We have included the imputation of common 1kG variants using the tag variants on the Illumina 1M Duo Version 3 array platform as an example. For this, use following:

```
cd Illumina_1M_Duo_1kG/1kG
chmod 755 setup.sh
./setup.sh
cd ../Illumina_1M_Duo
chmod 755 setup_run_models.sh
./setup_run_models.sh
```

These commands first download the 1000 Genomes Project genotype data, then downloads the array's variant loci, parses the 1000 Genomes Project variants. Then randomly splits the samples into 1500+1004 individuals and runs HEmpute-Train to build approximately 81,000 models.<br>

The parameter set for each variant are stored in files named "pooled_params_0.txt, pooled_params_1.txt, ..." where each line corresponds to a target SNP and the columns contain the linear model parameters. You can pool these files and get the whole parameters list in one file and sort the target variants at the same time:

```
cat pooled_params_*.txt | awk 'BEGIN{FS="\t"}{print $NF"\t"$0}' | sort -k1,1 -n | cut -f2- > pooled_params.txt
```

These commands can be used as a template to run other datasets.<br>

</html>
