#!/bin/sh

function usage() {
        echo "USAGE:"
        echo "sh seGMM.sh [-h] [-i <filename>] [-b <Bamlist>] [-c <Chromosome(x/y/b)>] [-s <y/n>] [-v] [-o <filename>]"
        echo "    -h help"
        echo "    -v genome version(default: hg19. If your vcf data is mapping to hg38, please use this parameter with no value!)"
        echo "    -i input vcf file"
        echo "    -b file contain sampleid and directory of bam file (no header)"
        echo "    -c Sex chromosome to use collect features. optional x,y,b"
        echo "    -s Including SRY gene or not. optional y,n"
        echo "    -o Prefix of output directory"
        exit -1
}

version=19

while getopts "hi:b:c:s:o:" opt
do
    case $opt in
            i)
                vcf=$OPTARG
                if [ ! -f $vcf ]
                then
                    echo "Error: the input vcf file $vcf doesn't exist!"
                    exit
                fi
                ;;
            b)
                bam=$OPTARG
                if [ ! -f $bam ]
                then
                    echo "Error: the bamlist file $bam doesn't exist!"
                    exit
                fi
                ;;
            c)
                chr=$OPTARG
                ;;
            s)
                use_SRY=$OPTARG
                ;;
            o)
                outputdir=$OPTARG
                if [ ! -d "$outputdir" ]
                then
                    echo "The output path `dirname $outputdir` doesn't exist!"
                    echo "Create the output path!"
                    mkdir "$outputdir"
                fi
                Outputdir=$(readlink -f $outputdir)
                ;;
            h)
                usage
                exit
                ;;
            v)
                version=38
                ;;
            :)
                echo "the option -$OPTARG require an arguement"
                usage
                exit 1
                ;;
            ?)
                usage
                exit 1
                ;;
    esac
done

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ [$chr = "x"] ]; then
    #Collect feature 1: X het rate
    plink --vcf $vcf --make-bed  --out $Outputdir/plink
    plink --bfile $Outputdir/plink --chr X --recode  --out $Outputdir/plink.X
    cat $Outputdir/plink.X.ped | awk 'BEGIN{OFS="\t"}{het=0;hom=0;missing=0;all=0;for(i=1;i<=((NF-6)/2);i++){if($(i*2+5)!=$(2*i+6)){het++;all++}else if($(i*2+5)==0){missing++;all++}else{hom++;all++}};if(all!=missing){print $1,het/(all-missing)}else{print $1,0};het=0;hom=0;missing=0;all=0}' | sort -k 1 >$Outputdir/XH.txt
        #Collect feature of X mapping rate
    if [ ! -d "$Outputdir/Reads_stat" ]; then
        mkdir $Outputdir/Reads_stat
    fi
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} X \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.X.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.total.stat
    cat $bam | awk '{print $1}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.X.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Xmap.txt
    paste $Outputdir/XH.txt $Outputdir/Xmap.txt | cut -f1,2,4 | awk 'BEGIN{OFS="\t";print "sampleid","XH","Xmap"}{print $0}' >$Outputdir/feature.txt
elif [ [$chr = "y"] ]; then
    if [ ! -d "$Outputdir/Reads_stat" ]; then
        mkdir $Outputdir/Reads_stat
    fi
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} Y \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.Y.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.total.stat
    cat $bam | awk '{print $1}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.Y.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Ymap.txt
    cat $Outputdir/Ymap.txt | awk 'BEGIN{OFS="\t";print "sampleid","Ymap"}{print $0}' >$Outputdir/feature.txt
    if [ [$use_SRY != "y"] ]; then
        if [ ! -d "$Outputdir/SRY" ]; then
            mkdir $Outputdir/SRY
        fi
        cat $bam | awk 'NR>=2{print $1,$4}' | while read id dir; do mosdepth -t 4 -Q 30 -b $script_dir/data/SRY_hg"$version".bed -n $Outputdir/SRY/$id $dir; done
        cat $bam | awk '{print $1}' | while read id; do cat $Outputdir/SRY/$id.mosdepth.summary.txt | tail -n 1 | awk 'BEGIN{OFS="\t"}{print "'$i'",$4}'; done | sort -k 1 >$Outputdir/SRY.txt
        paste $Outputdir/Ymap.txt $Outputdir/SRY.txt | cut -f1,2,4 | awk 'BEGIN{OFS="\t";print "sampleid","Ymap","SRY"}{print $0}' >$Outputdir/feature.txt
    fi
elif [ [$chr = "b"] ]; then
    plink --vcf $vcf --make-bed  --out $Outputdir/plink
    plink --bfile $Outputdir/plink --chr X --recode  --out $Outputdir/plink.X
    cat $Outputdir/plink.X.ped | awk 'BEGIN{OFS="\t"}{het=0;hom=0;missing=0;all=0;for(i=1;i<=((NF-6)/2);i++){if($(i*2+5)!=$(2*i+6)){het++;all++}else if($(i*2+5)==0){missing++;all++}else{hom++;all++}};if(all!=missing){print $1,het/(all-missing)}else{print $1,0};het=0;hom=0;missing=0;all=0}' | sort -k 1 >$Outputdir/XH.txt
        #Collect feature of X mapping rate
    if [ ! -d "$Outputdir/Reads_stat" ]; then
        mkdir $Outputdir/Reads_stat
    fi
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.total.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} X \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.X.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} Y \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.Y.stat
    cat $bam | awk '{print $1}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.X.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Xmap.txt
    cat $bam | awk '{print $1}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.Y.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Ymap.txt
    paste $Outputdir/XH.txt $Outputdir/Xmap.txt $Outputdir/Ymap.txt | cut -f1,2,4,6 | awk 'BEGIN{OFS="\t";print "sampleid","XH","Xmap",Ymap","XYratio"}{print $1,$2,$3,$2/$3,$4}' >$Outputdir/feature.txt
    if [ [$use_SRY != "y"] ]; then
        if [ ! -d "$Outputdir/SRY" ]; then
            mkdir $Outputdir/SRY
        fi
        cat $bam | awk 'NR>=2{print $1,$4}' | while read id dir; do mosdepth -t 4 -Q 30 -b $script_dir/data/SRY_hg"$version".bed -n $Outputdir/SRY/$id $dir; done
        cat $bam | awk '{print $1}' | while read id; do cat $Outputdir/SRY/$id.mosdepth.summary.txt | tail -n 1 | awk 'BEGIN{OFS="\t"}{print "'$i'",$4}'; done | sort -k 1 >$Outputdir/SRY.txt
        paste $Outputdir/XH.txt $Outputdir/Xmap.txt $Outputdir/Ymap.txt $Outputdir/SRY.txt | cut -f1,2,4,6,8 | awk 'BEGIN{OFS="\t";print "sampleid","XH","Xmap",Ymap","XYratio","SRY"}{print $1,$2,$3,$2/$3,$4,$5}' >$Outputdir/feature.txt
    fi
else
    echo "Error, please use the following command for running seGMM!"
    usage
    exit
fi

Rscript $script_dir/script/seGMM.r $Outputdir/feature.txt $Outputdir/
