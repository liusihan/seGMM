# seGMM
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

## What's new
Version 1.3.0 fixes some critical bugs that affected the performance of the predicting Karyotype and changed the names of some parameters.

Also, the CRAM file as an input file is now possible.


## Background
Computational tools have been developed to infer sex for genotype array, WES, and WGS data such as plink, seXY, and XYalign. Plink calculated the F coefficient with X chromosome heterozygosity to infer sex for genotype array data. seXY considered both X chromosome heterozygosity and Y chromosome missingness to infer sex in genotype array data by logistic regression. XYalign extract read count mapped to sex chromosomes and calculated the ratio of X and Y counts to infer sex in WES and WGS data. However, evaluation of the accuracy of these methods in targeted gene panel data is not yet entirely and improvements in sex inference from gene panel data are warranted. In addition, PLINK, seXY, and XYalign could not report sex chromosome abnormality. 

## Description
`seGMM` applies unsupervised learning to determine the individual gender based on integrated information of the X and Y chromosomes from TGS, WES, or WGS data, providing the classification of six sex chromosome karyotypes (XX, XY, XYY, XXY, XXX, and X).<br>

In clinical practice, individual patient is usually sequenced to get a molecular diagnosis. Hence, seGMM permits users to provide additional reference data, by combining the features from reference data and testing data, seGMM can ensure the accuracy for clinical applications. 

## Installation
### Robust install
In order to install the software and dependencies, we recommend using a dedicated conda environment and `seGMM` is available on conda (installation time: ~ 15min).
First, you will need the `Conda` Python distribution and package manager. 

```shell
# Download conda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Set permissions to execute
chmod +x Miniconda3-latest-Linux-x86_64.sh 	

# Execute. Make sure to "yes" to add the conda to your PATH
sh ./Miniconda3-latest-Linux-x86_64.sh 		

# Add channels
conda config --add channels bioconda
conda config --add channels conda-forge
```

After installing Miniconda, run the following commands to install seGMM and seGMM's dependencies:

```
conda create -n seGMM
conda activate seGMM
conda install -c sihanliu segmm
```
### Quick install
If `Conda` has been installed before, and the dependencies of seGMM is already installed in your environment, then you can quickly install seGMM from PyPI:

```shell
pip install seGMM
```

Dependencies
- Programming languages:
  * [Python](https://www.python.org/) >3
  * [R](https://www.r-project.org/) >= 3.5
 
- Commandline tools and packages:
  * [Plink](https://www.cog-genomics.org/plink/) >=1.9
  * [samtools](https://github.com/samtools/samtools) >=1.9
  * [mosdepth](https://github.com/brentp/mosdepth)
  * [parallel](https://www.gnu.org/software/parallel/)
  * [mclust](https://cran.r-project.org/web/packages/mclust/index.html)

Once the installation of seGMM has been completed, you can run:
```
seGMM -h
```
to print a list of all command-line options. If these commands fail with an error, something goes wrong during the installation process.


## Parameters

|Parameter|Type| Description|Required|
|---|---|---|---|
|--vcf/-vcf|character|Input VCF file (Either multi-sample or single-sample data. If the sample size is < 10, please combine with reference data for prediction analysis). |``true``|
|--input/-i|character| Input file contains sampleid and directory of bam/cram files. A text file contain two columns with no header and is split by space. The first column is the sample ID which matches the sample ID in the input vcf file. The order of the sample ID in the input file and the order of the sample ID in the VCF file can be inconsistent. **An example file has been provided in the test fold.**|``true``|
|--alignment_format/-a|character| Alignment format type for the input data.**Optional is {BAM, CRAM}.**|``true``|
|--reference_fasta/-R|character| Reference genome for **CRAM** support (if CRAM is used). [default: '']|``true``|
|--chromosome/-c|character|Sex chromosomes are used to collect features. **Optional is {xy,x,y}. If --reference is used, you can no longer use this parameter**|``false``|
|--type/-t|character|Sequencing type. **Optional is {TGS, WES,WGS}.** Note that if your **don't provide additional reference data, you must use --type.** If the data type is WGS or WES, seGMM will automatically calculate all 5 features, otherwise if your **sequencing type is TGS you have to choose which sex chromosome you want to use (--chromosome/-c) and tell seGMM the SRY gene is included or not (--SRY/-s)**|``false``|
|--output/-o|character|Prefix of output directory.|``true``|
|--genome/-g|character|Genome version. **Default is hg19. Option is {hg19,hg38}**.|``false``|                        
|--SRY/-s|boolean|If **True**, seGMM will calculate the mean depth of SRY gene. **Option is {True,False}**. |``false``|
|--reference_additional/-r|character|The path of the additional reference file contain features. We have provided two additional files (**1000G.WES.txt and 1000G.WGS.txt in the reference folder**). If **--reference is used, seGMM will automatically calculate the same features in the reference file. The file (tab split) must contain at least two features, and the column names must be: sampleid,XH,Xmap,Ymap,XYratio,SRY. The ordering of the columns is arbitrary, except for the first instance, which must be the sample name** |``false``|
|--uncertain_threshold/-u|numeric|The threshold for detecting outliers in GMM model. **Default is 0.1. The range of threshold is 0-1.**|``false``|
|--num_threshold/-n|numeric|Number of additional threads to use. Default is 1.|``false``|
|--Qulity/-q|numeric|Mapping quality threshold of reads to count. Default is 30.|``false``|

## Usage examples
```shell
## For WES data (CRAM). Using 20 cores
seGMM -vcf input.vcf -i cram.file -R GRCh38.fa -g hg38 -a CRAM -t WES -n 20 -o outputdir

## For TGS data
# The gene panel contains only genes located on the X chromosome
seGMM -vcf input.vcf -i bam.file -a BAM -t TGS -o outputdir -c x -s False

# The gene panel contains genes located on the X and Y chromosome, but don't contain SRY
seGMM -vcf input.vcf -i bam.file -a BAM -t TGS -o outputdir -c xy -s False

## With an additional reference file. Note the header of referenc file must like: sampleid,XH,Xmap,Ymap,XYraio,SRY. 
## And seGMM will automatically calculated the same features in the reference file. 
## We have provided two additinal files (1000G.WES.txt and 1000G.WGS.txt in reference folder).
seGMM -vcf test.vcf -i cram.file -R GRCh38.fa -g hg38 -a CRAM -t WES -r 1000G.WES.txt -o outputdir

```

## Test for seGMM
We have provide two reference file named ``1000G.WES.txt and 1000G.WGS.txt in reference folder``. Users can download these files and integrate with your own vcf and bam files(WES or WGS sequencing) to test the utility of seGMM. In addition, you can download test data from exon-targetted sequencing for 1000 genes from the 1000 Genomes Project in the ``test folder`` or in Google Drive (https://drive.google.com/drive/folders/1OrPD8t7CFg7ytdZb7EHVmXnNWCFPA1oj?usp=sharing). After download the file, you should make a ``bam.list`` file which contain sample ID and the full path of bam files. Then you can run
```shell
seGMM -vcf test.vcf -i Target.bam.list -t TGS -a BAM -o output -c x -s False
```
If everything goes well, you will see:
```
seGMM -vcf test.vcf -i Target.bam.list -t TGS -a BAM -o output -c x -s False
*********************************************************************
* seGMM
* Version 1.3.0
* (C) 2021-2026 Sihan Liu
* Research Institute of Rare Disease / West china hospital
* GNU General Public License v3
*********************************************************************

Beginning to generate features at Tue Nov  1 17:17:09 2022
>> Collected feature of X chromosome heterozygosity
    Finish generating features of X chromosome heterozygosity at Tue Nov  1 17:17:13 2022

>> Collected feature of X mapping rate
    Finish generating features of X mapping rate at Tue Nov  1 17:17:17 2022

>> Combine features into a single file

>> Running sample classification based on GMM model
WARNING: ignoring the environment value of R_HOME
[1] "There are 0 outliers samples based on prediction uncertainty"
character(0)
outliers
FALSE
   10

Analysis complete for seGMM at Tue Nov  1 17:17:18 2022
Total time elapsed: 9.49s

*********************************************************************
* Thanks for using seGMM!
* Report bugs to liusihan@wchscu.cn
* seGMM homepage: https://github.com/liusihan/seGMM
*********************************************************************
```

## Citation
If you use seGMM, please cite our paper (thanks!):
> Liu S, Zeng Y, Wang C, Zhang Q, Chen M, Wang X, Wang L, Lu Y, Guo H, Bu F. seGMM: A New Tool for Gender Determination From Massively Parallel Sequencing Data. Front Genet. 2022 Mar 3;13:850804. doi: 10.3389/fgene.2022.850804.

## License
This project is licensed under GNU GPL v3.

## Authors
Sihan Liu (West china hospital)


