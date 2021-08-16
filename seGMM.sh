#!/bin/sh
# MAKE SURE index.r RNA_processing.r and FDR.r ARE IN YOUR PATH
# FILL IN THESE PATHS

function usage() {
        echo "USAGE:"
        echo "sh seGMM.sh [-h] [-i <filename>] [-b <Bamlist>] [-c <Chromosome(x/y/b)>] [-s <y/n>] [-o <filename>]"
        echo "    -h help"
        echo "    -i input vcf file"
        echo "    -b file contain sampleid and directory of bam file (no header)"
        echo "    -c Sex chromosome to use collect features. optional x,y,b"
        echo "    -s Including SRY gene or not. optional y,n"
        echo "    -o Prefix of output directory"
        exit -1
}

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

if [ $chr = "x" ]; then
    #Collect feature 1: X het rate
    plink --vcf $vcf --recode  --out $Outputdir/plink
    cat $Outputdir/plink.ped | awk 'BEGIN{OFS="\t"}{het=0;hom=0;missing=0;all=0;for(i=1;i<=((NF-6)/2);i++){if($(i*2+5)!=$(2*i+6)){het++;all++}else if($(i*2+5)==0){missing++;all++}else{hom++;all++}};if(all!=missing){print $1,het/(all-missing)}else{print $1,0};het=0;hom=0;missing=0;all=0}' | sort -k 1 >$Outputdir/XH.txt
        #Collect feature of X mapping rate
    if [ ! -d "$Outputdir/Reads_stat" ]; then
        mkdir $Outputdir/Reads_stat
    fi
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} X \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.X.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.total.stat
    ls $Outputdir/Reads_stat/*.X.stat | awk '{split($1,a,".");split(a[1],b,"/");print b[length(b)]}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.X.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Xmap.txt
    paste $Outputdir/XH.txt $Outputdir/Xmap.txt | cut -f1,2,4 | awk 'BEGIN{OFS="\t";print "sampleid","XH","Xmap"}{print $0}' >$Outputdir/feature.txt
elif [ $chr = "y" ]; then
    if [ ! -d "$Outputdir/Reads_stat" ]; then
        mkdir $Outputdir/Reads_stat
    fi
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} Y \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.Y.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.total.stat
    ls $Outputdir/Reads_stat/*.Y.stat | awk '{split($1,a,".");split(a[1],b,"/");print b[length(b)]}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.Y.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Ymap.txt
    cat $Outputdir/Ymap.txt | awk 'BEGIN{OFS="\t";print "sampleid","Ymap"}{print $0}' >$Outputdir/feature.txt
    if [ $use_SRY != "y" ]; then
        if [ ! -d "$Outputdir/SRY" ]; then
            mkdir $Outputdir/SRY
        fi
        cat $bam | awk 'NR>=2{print $1,$4}' | while read id dir; do mosdepth -t 4 -Q 30 -b $basename/SRY.exon.bed -n $Outputdir/SRY/$id $dir; done
        ls $Outputdir/SRY/*.mosdepth.summary.txt | awk '{split($1,a,".");split(a[1],b,"/");print b[length(b)]}' | while read id; do cat $Outputdir/SRY/$id.mosdepth.summary.txt | tail -n 1 | awk 'BEGIN{OFS="\t"}{print "'$i'",$4}'; done | sort -k 1 >$Outputdir/SRY.txt
        paste $Outputdir/Ymap.txt $Outputdir/SRY.txt | cut -f1,2,4 | awk 'BEGIN{OFS="\t";print "sampleid","Ymap","SRY"}{print $0}' >$Outputdir/feature.txt
    fi
elif [ $chr = "b" ]; then
    plink --vcf $VCF --recode  --out $Outputdir/plink
    cat $Outputdir/plink.ped | awk 'BEGIN{OFS="\t"}{het=0;hom=0;missing=0;all=0;for(i=1;i<=((NF-6)/2);i++){if($(i*2+5)!=$(2*i+6)){het++;all++}else if($(i*2+5)==0){missing++;all++}else{hom++;all++}};if(all!=missing){print $1,het/(all-missing)}else{print $1,0};het=0;hom=0;missing=0;all=0}' | sort -k 1 >$Outputdir/XH.txt
        #Collect feature of X mapping rate
    if [ ! -d "$Outputdir/Reads_stat" ]; then
        mkdir $Outputdir/Reads_stat
    fi
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.total.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} X \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.X.stat
    cat $bam | awk '{print $1"\n"$2}' | parallel -j 10 --max-args 2 samtools view -bh -q 30 {2} Y \| samtools flagstat - \>$Outputdir/Reads_stat\/{1}.Y.stat
    ls $Outputdir/Reads_stat/*.X.stat | awk '{split($1,a,".");split(a[1],b,"/");print b[length(b)]}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.X.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Xmap.txt
    ls $Outputdir/Reads_stat/*.Y.stat | awk '{split($1,a,".");split(a[1],b,"/");print b[length(b)]}' | while read id; do paste -d' ' $Outputdir/Reads_stat/$id.Y.stat $Outputdir/Reads_stat/$id.total.stat | sed -n '5p' | cut -d ' ' -f1,8,15| sed 's/^/'$id' /'; done | awk 'BEGIN{OFS="\t"}{print $1,$2/$3}' | sort -k 1 >$Outputdir/Ymap.txt
    paste $Outputdir/XH.txt $Outputdir/Xmap.txt $Outputdir/Ymap.txt | cut -f1,2,4,6 | awk 'BEGIN{OFS="\t";print "sampleid","XH","Xmap",Ymap","XYratio"}{print $1,$2,$3,$2/$3,$4}' >$Outputdir/feature.txt
    if [ $use_SRY != "y" ]; then
        if [ ! -d "$Outputdir/SRY" ]; then
            mkdir $Outputdir/SRY
        fi
        cat $bam | awk 'NR>=2{print $1,$4}' | while read id dir; do mosdepth -t 4 -Q 30 -b $basename/SRY.exon.bed -n $Outputdir/SRY/$id $dir; done
        ls $Outputdir/SRY/*.mosdepth.summary.txt | awk '{split($1,a,".");split(a[1],b,"/");print b[length(b)]}' | while read id; do cat $Outputdir/SRY/$id.mosdepth.summary.txt | tail -n 1 | awk 'BEGIN{OFS="\t"}{print "'$i'",$4}'; done | sort -k 1 >$Outputdir/SRY.txt
        paste $Outputdir/XH.txt $Outputdir/Xmap.txt $Outputdir/Ymap.txt $Outputdir/SRY.txt | cut -f1,2,4,6,8 | awk 'BEGIN{OFS="\t";print "sampleid","XH","Xmap",Ymap","XYratio","SRY"}{print $1,$2,$3,$2/$3,$4,$5}' >$Outputdir/feature.txt
    fi
else
    echo "Please select the sex chromosome you want to use! Optional value is x, y and b(both)."
    usage
    exit
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
Rscript $script_dir/seGMM.r $Outputdir/feature.txt $Outputdir/