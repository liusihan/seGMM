# seGMM
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

## Background
Computational tools have been developed to infer sex for genotype array, WES and WGS data such as plink, seXY and SEXCMD. Plink calculated F coefficient with X chromosome heterozygosity to infer sex for genotype array data. seXY considered both X chromosome heterozygosity and Y chromosome missingness to infer sex in genotype array data by logistic regression. SEXCMD can extract read count which mapped to sex-specific marker sequences from syntenic regions and calculated ratio of X and Y counts to infer sex in WES and WGS data. However, evaluation the accuracy of these methods in panel data is not yet fully and improvements in sex inference from gene panel data are warranted.

## Description
`seGMM` is a tool that determines the gender of a sample from the called genotype data integrated BAM files and jointly considers information on the X and Y chromosomes in diverse genotype data including `panel data`. seGMM apply `Gaussian Mixture Model (GMM)` clustering to classify the samples into two clusters.<br>

![](https://github.com/liusihan/seGMM/blob/main/Workflow.GIF)  

### If you use seGMM, please cite our preprint (thanks!):
>Sihan Liu (2021) seGMM: a new tool to infer sex from massively parallel sequencing data. bioRxiv

## Installation
In order to download `seGMM`, you should clone this repository via the commands

```
git clone https://github.com/liusihan/seGMM
cd seGMM
```
In order to install the software and R packages dependencies, you will need the `Anaconda` Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with seGMM's dependencies:

```
conda env create -f environment.yaml
source activate seGMM
```

Once the above has completed, you can run:
```
sh seGMM.sh -h
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
