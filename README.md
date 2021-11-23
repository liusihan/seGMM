# seGMM
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

## Background
Computational tools have been developed to infer sex for genotype array, WES and WGS data such as plink, seXY and XYalign. Plink calculated F coefficient with X chromosome heterozygosity to infer sex for genotype array data. seXY considered both X chromosome heterozygosity and Y chromosome missingness to infer sex in genotype array data by logistic regression. XYalign extract read count mapped to sex chromosomes and calculated the ratio of X and Y counts to infer sex in WES and WGS data. However, evaluation the accuracy of these methods in targeted gene panel data is not yet fully and improvements in sex inference from gene panel data are warranted.

## Description
`seGMM` is a tool that determines the gender of a sample from the called genotype data integrated BAM files and jointly considers information on the X and Y chromosomes in diverse genotype data including `panel data`. seGMM apply `Gaussian Mixture Model (GMM)` clustering to classify the samples into two clusters.<br>

Importantly, in clinical practice, individual patient or trio samples are usually sequenced to get a molecular diagnosis. Hence, seGMM permits users to provide an additional reference data, by combining the features from reference data, seGMM can ensure the accuracy for clinical application. Besides, seGMM can throw the exceptions with an uncertain classification, indicating potential events of sex chromosome abnormity.

![](https://github.com/liusihan/seGMM/blob/main/Workflow.GIF)  

## Installation
### Robust install
From PyPI:

```
pip install seGMM
```

Dependencies:
- Python >= 3
- Plink >= 1.9
- parallel
- samtools >= 1.9
- mosdepth
- R >= 3.5
- mclust
  
### Quick install
In order to install the software and dependencies, we recommend using a dedicated conda environment and `seGMM` is available on bioconda. (installation time: ~ 10min)
First you will need the `Conda` Python distribution and package manager. 

```
# Download conda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Set permissions to execute
chmod +x Miniconda3-latest-Linux-x86_64.sh 	

# Execute. Make sure to "yes" to add the conda to your PATH
sh ./Miniconda3-latest-Linux-x86_64.sh 		

# Add channels
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

After installing Miniconda, run the following commands to create an environment with seGMM's dependencies:

```
conda install -c bioconda seGMM
source activate seGMM
```

Once the above has completed, you can run:
```
seGMM -h
```
to print a list of all command-line options. 

## Usage
```
sh seGMM.sh -i vcf.gz -b bam.list -c x(y/b) -s n(y) -o output
```
- `-i`: input vcf file (vcf or vcf.gz).
- `-b`: file contain sample ID and directory of bam file(no header, space separated).
- `-c`: choose chromosome you want to use. Optional is x and y and b(both).
- `-s`: Including SRY gene or not.Optional is y and n.
- `-o`: Prefix of output file. (If not exist, seGMM will create it!)
- `-v`: Genome version(default: hg19. If your vcf data is mapping to hg38, please use this parameter with no value!). 

## License
MIT ©
