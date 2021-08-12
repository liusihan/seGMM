# seGMM
## Description
seGMM is a tool that determines the gender of a sample from the called genotype data integrated BAM files and jointly considers information on the X and Y chromosomes in diverse genotype data including `panel data`. seGMM apply `Gaussian Mixture Model (GMM)` clustering to classify the samples into two clusters. The workflow of seGMM is shown in `Figule 1`.<br>
![](https://github.com/liusihan/seGMM/blob/main/Workflow.GIF)  
### If you use seGMM, please cite our preprint (thanks!):
>Sihan Liu (2021) seGMM: a new tool to infer sex from massively parallel sequencing data. bioRxiv

## Installation
### Quick installation
```
git clone https://github.com/liusihan/seGMM
```

### Robust installation (conda)
```
conda env create -f environment.yaml
```

## Usage
```
sh seGMM.sh vcf.gz bam.list output
```
