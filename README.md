# Secure and practical genotype imputation methods via homomorphic encryption

Genotype imputation is a fundamental step in genomic data analysis such as GWAS where missing variant genotypes are predicted using the existing genotypes of the nearby ‘tag’ variants. Imputation greatly decreases the genotyping cost while providing high quality estimates of the common variant genotypes. As the population panels are increasing, e.g. TOPMED Project, the genotype imputation is becoming more accurate but requires high computational power. While researchers can outsource the genotype imputation, privacy concerns may prohibit the genetic data sharing with an untrusted imputation service. To address this problem, we developed the first fully secure genotype imputation utilizing ultra-fast homomorphic encryption (HE) techniques. In HE-based methods, the genotype data is encrypted by the data owner and can never be decrypted in the course of imputation. We compare the secure imputation with 3 other state-of-the-art methods and found that HE-based methods provide full genetic data security at the expense of slightly lower accuracy. In addition, HE-based methods have comparable and even lower time and memory requirements than the state-of-the-art methods. We provide 5 different implementations and workflows that make use of 3 different HE schemes (BFV, CKKS, TFHE), developed by the top contestants of the iDASH19 Genome Privacy Challenge. Our results provide strong evidence that the HE-based methods can be practically used for large scale genetic data analysis. In addition, the publicly available codebases provide a reference for the development of secure genomic data analysis methods.

## Data Repository
---
This directory contains all the data used for training and testing. 

## Secure genotype imputation protocol implementations 
---
An Implementation of secure genotype imputation protocol from each team is available at the following repository:

- **UTHealth**: [https://github.com/K-miran/HEmpute](https://github.com/K-miran/HEmpute).
- **Chimera-TFHE**: [https://github.com/ssmiler/idash2019_2](https://github.com/ssmiler/idash2019_2).
- **EPFL**:
- **SNU**:
