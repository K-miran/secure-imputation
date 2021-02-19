# Ultra-Fast Homomorphic Encryption Models enable Secure Outsourcing of Genotype Imputation

Genotype imputation is a fundamental step in genomic data analysis such as GWAS, where missing variant genotypes are predicted using the existing genotypes of nearby ‘tag’ variants. Imputation greatly decreases the genotyping cost and provides high-quality estimates of common variant genotypes. As population panels increase, e.g., the TOPMED Project, genotype imputation is becoming more accurate, but it requires high computational power. Although researchers can outsource genotype imputation, privacy concerns may prohibit genetic data sharing with an untrusted imputation service. To address this problem, we developed the first fully secure genotype imputation by utilizing ultra-fast homomorphic encryption (HE) techniques that can evaluate millions of imputation models in seconds. In HE-based methods, the genotype data is end-to-end encrypted, i.e., encrypted in transit, at rest, and, most importantly, in analysis, and can be decrypted only by the data owner. We compared secure imputation with three other state-of-the-art non-secure methods under different settings. We found that HE-based methods provide full genetic data security with comparable or slightly lower accuracy. In addition, HE-based methods have time and memory requirements that are comparable and even lower than the non-secure methods. We provide five different implementations and workflows that make use of three cutting-edge HE schemes (BFV, CKKS, TFHE) developed by the top contestants of the iDASH19 Genome Privacy Challenge. Our results provide strong evidence that HE-based methods can practically perform resource-intensive computations for high throughput genetic data analysis. In addition, the publicly available codebases provide a reference for the development of secure genomic data analysis methods.

## Data Repository
You can see all the (compressed) data used for training and testing in the 'data' directory. 

## Secure genotype imputation protocol implementations 

You can download the main project with all contained submodules to your local computer by using the clone command: 
```sh
git clone --recurse-submodules https://github.com/K-miran/secure-imputation.git
```
You can specify a submodule to be initialized and cloned by runing the command:
```sh
git submodule update --init --recursive <pathspec>
```

Alternatively, you can also find an Implementation of secure genotype imputation protocol from each team as the following repository:

- **UTHealth & Microsoft Research**: [https://github.com/K-miran/HEmpute](https://github.com/K-miran/UTMSR_HEmpute)
- **Chimera-TFHE**: [https://github.com/ssmiler/idash2019_2](https://github.com/ssmiler/idash2019_2)
- **EPFL**:[https://github.com/ldsec/Lattigo_genotype_imputation](https://github.com/ldsec/Lattigo_genotype_imputation)
- **SNU**: [https://github.com/idashSNU/Imputation](https://github.com/idashSNU/Imputation)
