# seGMM
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

## Background
Computational tools have been developed to infer sex for genotype array, WES and WGS data such as plink, seXY and XYalign. Plink calculated F coefficient with X chromosome heterozygosity to infer sex for genotype array data. seXY considered both X chromosome heterozygosity and Y chromosome missingness to infer sex in genotype array data by logistic regression. XYalign extract read count mapped to sex chromosomes and calculated the ratio of X and Y counts to infer sex in WES and WGS data. However, evaluation the accuracy of these methods in targeted gene panel data is not yet fully and improvements in sex inference from gene panel data are warranted.

## Description
`seGMM` is a tool that determines the gender of a sample from the called genotype data integrated BAM files and jointly considers information on the X and Y chromosomes in diverse genotype data including `panel data`. seGMM apply `Gaussian Mixture Model (GMM)` clustering to classify the samples into two clusters.<br>

Importantly, in clinical practice, individual patient or trio samples are usually sequenced to get a molecular diagnosis. Hence, seGMM permits users to provide an additional reference data, by combining the features from reference data, seGMM can ensure the accuracy for clinical application. Besides, seGMM can throw the exceptions with an uncertain classification, indicating potential events of sex chromosome abnormity.

![](https://github.com/liusihan/seGMM/blob/main/Workflow.GIF)  

## Installation
### Quick install
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
  
### Robust install
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
to print a list of all command-line options. If these commands fail with an error, then something as gone wrong during the installation process.

Now we describe the different parameters needed in seGMM.

##Parameters
|Parameter|Type| Description|Required|
|---|---|---|---|
|--input/-i|character|Path of the input vcf files.|``true``|
|--bam/-b|character|Path of file contain the **sampleID and Full path of correspanding bam file. No header.**|``true``|
|--chromosome/-c|character|Sex chromosome to use collect features. **Optional is {xy,x,y}**|``false``|
|--type/-t|character|Study type. Note that if your **don't provide an additional reference data, you must use --type.** If the data type is WGS or WES, seGMM will automatic calculated all 5 features, otherwise if your **data type is TGS you have to choice which sex chromosome you want to use (--chromosome/-c) and tell seGMM the SRY gene is included or not (--SRY/-s)**|``false``|
|--output/-o                   character  Genome version. Default is hg19.                                                              ``false #{hg19,hg38}``
|--genome/-g                   character  Set to true if non standard chromosomes are to be kept for further analysis.                  ``true``                        
|--SRY/-s                      boolean    If false, seGMM will not calculate the mean depth of SRY gene.                                ``false``
|--reference/-r                character  The additional reference file contain features. We have provided two additinal files(1000G_WES and 1000G_WGS).      ``false``
|--uncertain_threshold/-u      numeric    The threshold for detecting outliers in GMM model. Default is 0.1. The range of threshold is 0-1.   ``false``
|--num_threshold/-n            numeric    Number of additional threads to use. Default is 1.                                           ``false``
|--Qulity/-q                   numeric    Mapping quality threshold of reads to count. Default is 30.                                  ``false``
|--XH/-x                       character  With a provided external reference data, using this parameter with no value, seGMM will calculated XH.                    ``false``
|--Xmap/-m                     character  With a provided external reference data, using this parameter with no value, seGMM will calculated Xmap.                                                              ``false``
|--Ymap/-y                     character  With a provided external reference data, using this parameter with no value, seGMM will calculated Ymap.                  ''false''
============================  =========  ============================================================================================  ======

## Citation
If you use the software or the LD Score regression intercept, please cite

## License
This project is licensed under GNU GPL v3.

## Authors
Sihan Liu (West china hospital)


