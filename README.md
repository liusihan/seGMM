# seGMM
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
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

## Usage
```
sh seGMM.sh vcf.gz bam.list output
```

## License
MIT © Richard McRichface
