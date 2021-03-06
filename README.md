# seGMM
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

## Background
Computational tools have been developed to infer sex for genotype array, WES and WGS data such as plink, seXY and XYalign. Plink calculated F coefficient with X chromosome heterozygosity to infer sex for genotype array data. seXY considered both X chromosome heterozygosity and Y chromosome missingness to infer sex in genotype array data by logistic regression. XYalign extract read count mapped to sex chromosomes and calculated the ratio of X and Y counts to infer sex in WES and WGS data. However, evaluation the accuracy of these methods in targeted gene panel data is not yet fully and improvements in sex inference from gene panel data are warranted.In addition, PLINK, seXY, and XYalign could not report sex chromosome abnormality. 

## Description
`seGMM` applies unsupervised learning to determine the individual gender based on integrated information of the X and Y chromosomes from TGS, WES, or WGS data, providing the classification of six sex chromosome karyotypes (XX, XY, XYY, XXY, XXX, and X).<br>

Importantly, in clinical practice, individual patient is usually sequenced to get a molecular diagnosis. Hence, seGMM permits users to provide an additional reference data, by combining the features from reference data and testing data, seGMM can ensure the accuracy for clinical application. 

## Installation
### Robust install
In order to install the software and dependencies, we recommend using a dedicated conda environment and `seGMM` is available on conda (installation time: ~ 15min).
First you will need the `Conda` Python distribution and package manager. 

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

Once the installation of seGMM has completed, you can run:
```
seGMM -h
```
to print a list of all command-line options. If these commands fail with an error, then something as gone wrong during the installation process.


## Parameters

|Parameter|Type| Description|Required|
|---|---|---|---|
|--input/-i|character|The input vcf file which contains all the samples you want to infer gender. |``true``|
|--bam/-b|character|"Bam.file". A text file which contain two columns with no header and split by space. The first column is the sample ID which match the sample ID in the input vcf file. The order of the sample ID in bam.file and the order of the sample ID in the VCF file can be inconsistent. **An example file has been provided in the test fold.**|``true``|
|--chromosome/-c|character|Sex chromosome to use collect features. **Optional is {xy,x,y}. If --reference is used, you can no longer use this parameter**|``false``|
|--type/-t|character|Sequencing type. Note that if your **don't provide an additional reference data, you must use --type.** If the data type is WGS or WES, seGMM will automatic calculated all 5 features, otherwise if your **sequencing type is TGS you have to choice which sex chromosome you want to use (--chromosome/-c) and tell seGMM the SRY gene is included or not (--SRY/-s)**|``false``|
|--output/-o|character|Prefix of output directory.|``true``|
|--genome/-g|character|Genome version. **Default is hg19. Option is {hg19,hg38}**.|``false``|                        
|--SRY/-s|boolean|If **True**, seGMM will calculate the mean depth of SRY gene. **Option is {True,False}**. |``false``|
|--reference/-r|character|The path of additional reference file contain features. We have provided two additinal files (**1000G.WES.txt and 1000G.WGS.txt in reference folder**). If **--reference is used, seGMM will automatically calculated the same features in the reference file. The file (tab split) must contain at least two features, and the column names must be: sampleid,XH,Xmap,Ymap,XYratio,SRY. The ordering of the columns is arbitrary, except for the first instance, which must be the sample name** |``false``|
|--uncertain_threshold/-u|numeric|The threshold for detecting outliers in GMM model. **Default is 0.1. The range of threshold is 0-1.**|``false``|
|--num_threshold/-n|numeric|Number of additional threads to use. Default is 1.|``false``|
|--Qulity/-q|numeric|Mapping quality threshold of reads to count. Default is 30.|``false``|

## Usage examples
```shell
## For WES and WGS data. Using 20 cores
seGMM -i input.vcf -b bam.file -t WES -o outputdir -n 20

## For TGS data
# The gene panel contains only genes located on the X chromosome
seGMM -i input.vcf -b bam.file -t TGS -o outputdir -c x -s False

# The gene panel contains genes located on the X and Y chromosome, but don't contain SRY
seGMM -i input.vcf -b bam.file -t TGS -o outputdir -c xy -s False

## With an additional reference file. Note the header of referenc file must like: sampleid,XH,Xmap,Ymap,XYraio,SRY. 
## And seGMM will automatically calculated the same features in the reference file. 
## We have provided two additinal files (1000G.WES.txt and 1000G.WGS.txt in reference folder).
seGMM -i input.vcf -b bam.file -r reference.txt -o outputdir

```

## Test for seGMM
We have provide two reference file named ``1000G.WES.txt and 1000G.WGS.txt in reference folder``. Users can download these files and integrate with your own vcf and bam files(WES or WGS sequencing) to test the utility of seGMM. In addition, you can download test data from exon-targetted sequencing for 1000 genes from the 1000 Genomes Project in the ``test folder``. After download the file, you should make a ``bam.list`` file which contain sample ID and the full path of bam files. Then you can run
```shell
seGMM -i test.vcf -b bam.list -t TGS -c xy -s False -o seGMM_test
```
If everything goes well, you will see:
```
seGMM -i test.vcf -b bam.list -t TGS -c xy -s False -o seGMM_test
*********************************************************************
* seGMM
* Version 1.2.1
* (C) 2021-2026 Sihan Liu
* Research Institute of Rare disease / West china hospital
* GNU General Public License v3
*********************************************************************

Warning, the output file is not exist, seGMM creates the output folder of seGMM_test first!
Beginning generate features at Thu Nov 25 14:10:23 2021
>> Collected feature of X chromosome heterozygosity
    Finish generate features of X chromosome heterozygosity at Thu Nov 25 14:10:24 2021

>> Collected feature of X mapping rate
    Finish generate features of X mapping rate at Thu Nov 25 14:10:37 2021

>> Collected feature of Y mapping rate
    Finish generate features of Y mapping rate at Thu Nov 25 14:10:38 2021

>> Combine features into a single file

>> Running sample classfication based on GMM model
WARNING: ignoring environment value of R_HOME
[1] "There are 0 outliers samples based on prediction uncertainty"
character(0)
outliers
FALSE
   10

Analysis complete for seGMM at Thu Nov 25 14:10:38 2021
Total time elapsed: 15.59s

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


