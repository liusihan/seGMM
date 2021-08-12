# seGMM
seGMM is a tool that determines the gender of a sample from the called genotype data integrated BAM files and jointly considers information on the X and Y chromosomes in diverse genotype data including `panel data`. seGMM apply `Gaussian Mixture Model (GMM)` clustering to classify the samples into two clusters. The workflow of seGMM is shown in `Figule 1`.<br>
### If you use echolocatoR, please cite our preprint (thanks!):
>Sihan Liu (2021) seGMM: a new tool to infer sex from massively parallel sequencing data. bioRxiv
## Usage
''
sh seGMM.sh vcf.gz bam.list output
''
