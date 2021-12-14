#!/usr/bin/env python
############################################################################
# Created Time: 2021-12-13
# File Name: seGMM
# Last Change:.
# Description: a new tool to infer sex from massively parallel sequencing data
############################################################################

import time, os, re, sys, argparse, traceback,subprocess
from functools import reduce
from pathlib import Path

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

def runcmd(command):
    try:
        return_info = subprocess.run(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if return_info.returncode != 0:
            print("When running: {cmd}".format(cmd=command))
            print("An error was occured, please check the parameters!")
            sys.exit()
    except Exception as e:
        print("When running: {cmd}".format(cmd=command))
        print("An error was occured, please check the parameters!")
        sys.exit()

def collect_XH(input_vcf,outdir):
    print(">> Collected feature of X chromosome heterozygosity")
    cmd = "plink --vcf " + input_vcf + " --make-bed --out " + outdir+"/plink"
    runcmd(cmd)
    cmd = "plink --bfile " + outdir+"/plink"+ " --chr X --recode --out " + outdir+"/plink.X"
    runcmd(cmd)
    pedfile=outdir+"/plink.X.ped"
    outfile=outdir+"/XH.txt"
    if os.path.isfile(outfile):
        cmd="rm "+outfile
        runcmd(cmd)
    XH_list=[]
    with open(outfile,'a+')as xh:
        with open(pedfile) as f:
            line=f.readline()
            while line:
                lines=line.rstrip("\n").split(" ")
                het=0
                hom=0
                missing=0
                all_num=0
                for i in range(0,int((len(lines)-6)/2)):
                    if lines[i*2+4]!=lines[2*i+5]:
                        het+=1
                        all_num+=1
                    elif lines[i*2+4]==0:
                        missing+=1
                        all_num+=1
                    else:
                        hom+=1
                        all_num+=1
                if all_num!=missing:
                    xh.write(str(lines[0])+'\t'+str(het/(all_num-missing))+'\n')
                else:
                    xh.write(str(lines[0])+'\t'+float(0)+'\n')
                line=f.readline()
        f.close()
    xh.close()
    cmd="cat "+ outdir +"/XH.txt "+"| sort -n -k 1 >"+ outdir+"/XH.sorted.txt"
    runcmd(cmd)
    print('    Finish generate features of X chromosome heterozygosity at {T} \n'.format(T=time.ctime()))

def collect_Xmap(bamfile,quality,num_threshold,outdir):
    print(">> Collected feature of X mapping rate")
    if not os.path.exists(outdir+"/"+"Read_stat"):
        os.makedirs(outdir+"/"+"Read_stat")
    cmd="cat " + bamfile + ''' | awk '{print $1"\\n"$2}' | parallel -j ''' + num_threshold +" --max-args 2 samtools view -@ 10 -bh -q " + quality + " {2} X \| samtools flagstat -@ 10 - \>" +outdir+"\/"+"Read_stat"+ "\/{1}.X.stat"
    runcmd(cmd)
    cmd="cat " + bamfile + ''' | awk '{print $1}' | while read id; do cat ''' + outdir+"/"+"Read_stat/$id.X.stat" + ''' | awk 'BEGIN{OFS="\\t"}NR==9{print "'$id'",$1}' ;done >'''+outdir+"/Read_stat/Xmap.txt"
    runcmd(cmd)
    if not os.path.isfile(outdir+"/total.txt"):
        collect_total(bamfile,quality,num_threshold,outdir)
        cmd = "paste " + outdir+"/Read_stat/Xmap.txt " + outdir+"/total.txt" + ''' | awk 'BEGIN{OFS="\\t"}{print $1,$2/$4}' >''' + outdir+"/Xmap.txt"
        runcmd(cmd)
    else:
        cmd = "paste " + outdir+"/Read_stat/Xmap.txt " + outdir+"/total.txt" + ''' | awk 'BEGIN{OFS="\\t"}{print $1,$2/$4}' >''' + outdir+"/Xmap.txt"
        runcmd(cmd)
    print('    Finish generate features of X mapping rate at {T} \n'.format(T=time.ctime()))

def collect_Ymap(bamfile,quality,num_threshold,outdir):
    print(">> Collected feature of Y mapping rate")
    if not os.path.exists(outdir+"/"+"Read_stat"):
        os.makedirs(outdir+"/"+"Read_stat")
    cmd="cat " + bamfile + ''' | awk '{print $1"\\n"$2}' | parallel -j ''' + num_threshold +" --max-args 2 samtools view -@ 10 -bh -q " + quality + " {2} Y \| samtools flagstat -@ 10 - \>" +outdir+"/"+"Read_stat"+ "\/{1}.Y.stat"
    runcmd(cmd)
    cmd="cat " + bamfile + ''' | awk '{print $1}' | while read id; do cat ''' + outdir+"/"+"Read_stat/$id.Y.stat" + ''' | awk 'BEGIN{OFS="\\t"}NR==9{print "'$id'",$1}' ;done >'''+outdir+"/Read_stat/Ymap.txt"
    runcmd(cmd)
    if not os.path.isfile(outdir+"/total.txt"):
        collect_total(bamfile,quality,num_threshold,outdir)
        cmd = "paste " + outdir+"/Read_stat/Ymap.txt " + outdir+"/total.txt" + ''' | awk 'BEGIN{OFS="\\t"}{print $1,$2/$4}' >''' + outdir+"/Ymap.txt"
        runcmd(cmd)
    else:
        cmd = "paste " + outdir+"/Read_stat/Ymap.txt " + outdir+"/total.txt" + ''' | awk 'BEGIN{OFS="\\t"}{print $1,$2/$4}' >''' + outdir+"/Ymap.txt"
        runcmd(cmd)
    print('    Finish generate features of Y mapping rate at {T} \n'.format(T=time.ctime()))

def collect_total(bamfile,quality,num_threshold,outdir):
    if not os.path.exists(outdir+"/"+"Read_stat"):
        os.makedirs(outdir+"/"+"Read_stat")
    cmd="cat " + bamfile + ''' | awk '{print $1"\\n"$2}' | parallel -j ''' + num_threshold +" --max-args 2 samtools view -@ 10 -bh -q " + quality + " {2} \| samtools flagstat -@ 10 - \>" +outdir+"/"+"Read_stat"+ "\/{1}.total.stat"
    runcmd(cmd)
    cmd="cat " + bamfile + ''' | awk '{print $1}' | while read id; do cat ''' + outdir+"/"+"Read_stat/$id.total.stat" + ''' | awk 'BEGIN{OFS="\\t"}NR==9{print "'$id'",$1}' ;done >'''+outdir+"/total.txt"
    runcmd(cmd)

def collect_SRY(bamfile,quality,genome_version,outdir):
    print(">> Collected feature of mean depth of SRY gene")
    if not os.path.exists(outdir+"/"+"SRY"):
        os.makedirs(outdir+"/"+"SRY")
    cmd="cat " + bamfile + ''' | awk '{print $1,$2}' | while read id dir; do mosdepth -t 4 -Q ''' + quality + " -b " + str(Path(__file__).absolute().parent)+"/data/SRY_"+ genome_version + ".bed -n "+ outdir+"/"+"SRY/$id $dir; done"
    runcmd(cmd)
    cmd="cat "+bamfile+''' | awk '{print $1}' | while read id; do cat '''+outdir+"/SRY/$id.mosdepth.summary.txt | "+'''awk 'BEGIN{OFS="\\t"}NR==2{out1=$4}NR==3{out2=$4}END{if(out1!=0){print "'$id'",out1}else if(out2!=0){print "'$id'",out2}else{print "'$id'",0}}'; done | sort -n -k 1 >''' + outdir+"/SRY.txt"
    runcmd(cmd)
    print('    Finish generate features of mean depth of SRY gene at {T} \n'.format(T=time.ctime()))


def with_reference(feature, input_vcf,bamfile,quality,num_threshold,genome_version,outdir):
    if feature == "XH":
        collect_XH(input_vcf,outdir)
    if feature == "Xmap":
        collect_Xmap(bamfile,quality,num_threshold,outdir)
    if feature == "Ymap":
        collect_Ymap(bamfile,quality,num_threshold,outdir)
    if feature == "SRY":
        collect_SRY(bamfile,quality,genome_version,outdir)
    if feature == "XYratio":
        Xmap = outdir+"/Xmap.txt"
        Ymap = outdir+"/Ymap.txt"
        if not os.path.isfile(Xmap):
            collect_Xmap(bamfile,quality,num_threshold,outdir)
        if not os.path.isfile(Ymap):
            collect_Ymap(bamfile,quality,num_threshold,outdir)
        XYratio = outdir+"/XYratio.txt"
        cmd = "paste "+Xmap+" " + Ymap + ''' | awk 'BEGIN{OFS="\\t"}{print $1,$2/$4}' >''' + XYratio
        runcmd(cmd)

def main():
    description = 'seGMM is a new tool to infer sex from massively parallel sequencing data. \n Written by Sihan Liu, liusihan@wchscu.cn. \n Please contact the authors for commercial use. \n Copyright (C) 2021 Institute of rare diseases\n============================================================================'
    __version__ = '1.2.1'
    header = "\n"
    header = "*********************************************************************\n"
    header += "* seGMM\n"
    header += "* Version {V}\n".format(V=__version__)
    header += "* (C) 2021-2026 Sihan Liu\n"
    header += "* Research Institute of Rare disease / West china hospital\n"
    header += "* GNU General Public License v3\n"
    header += "*********************************************************************\n"

    end = "*********************************************************************\n"
    end += "* Thanks for using seGMM!\n"
    end += "* Report bugs to liusihan@wchscu.cn\n"
    end += "* seGMM homepage: https://github.com/liusihan/seGMM\n"
    end += "*********************************************************************"

#Usage
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("--input","-i",required = True,help = "Input of VCF file.")
    parser.add_argument("--bam","-b",required = True,help = "File contain sampleid and directory of bam file (no header)")
    parser.add_argument("--chromosome","-c",required = False,help = "Sex chromosome to use collect features. ",choices=["xy","x","y"])
    parser.add_argument("--type","-t",required = False,help = "Sequencing type. Note that if your don't provide an additional reference data, you must use --type. If the data type is WGS or WES, seGMM will automatic calculated all 5 features, otherwise if your data type is TGS you have to choice which sex chromosome you want to use and tell seGMM the SRY gene is included or not!",choices=["TGS","WES","WGS"])
    parser.add_argument("--output","-o",required = True,help = "Prefix of output directory.")
    parser.add_argument("--genome","-g",required = False,help = "Genome version. Default is hg19. ",default='hg19',choices=["hg19","hg38"])
    parser.add_argument("--SRY","-s",required = False,help = "Including SRY gene or not.",choices=["True","False"])
    parser.add_argument("--reference","-r",required = False,help = "Reference file contain features.")
    parser.add_argument("--uncertain_threshold","-u",required = False,help = "The threshold for detecting outliers in GMM model. Default is 0.1. The range of threshold is 0-1!",default=0.1)
    parser.add_argument("--num_threshold","-n",required = False,help = "Number of additional threads to use. Default is 1.", default=1)
    parser.add_argument("--Quality","-q",required = False,help = "Mapping quality threshold of reads to count. Default is 30.", default=30)
    genome="hg19"
    uncertain_threshold="0.1"
    num_threshold="1"
    quality="30"
    args = parser.parse_args()
    print(header)
    try:
        if args.chromosome:
            chromosome=str(args.chromosome)
        if args.genome:
            genome=str(args.genome)
        if args.SRY:
            SRY=str(args.SRY)
        if args.uncertain_threshold:
            uncertain_threshold=str(args.uncertain_threshold)
        if args.num_threshold:
            num_threshold=str(args.num_threshold)
        if args.Quality:
            quality=str(args.Quality)
        if not os.path.isfile(args.input):
            sys.exit('Error, the input vcf file is not exist, please check that you have provided the correct path!')
        if not os.path.isfile(args.bam):
            sys.exit('Error, the input bam file is not exist, please check that you have provided the correct path!')
        if not os.path.exists(args.output):
            os.makedirs(args.output)
            print("Warning, the output file is not exist, seGMM creates the output folder of {s} first!".format(s=args.output))
        if os.path.isfile(args.input) and os.path.isfile(args.bam):
            if args.reference is None:
                if args.type == "WES" or args.type == "WGS":
                    print('Beginning generate features at {T}'.format(T=time.ctime()))
                    start_time = time.time()
                    collect_XH(args.input,args.output)
                    collect_Xmap(args.bam,quality,num_threshold,args.output)
                    collect_Ymap(args.bam,quality,num_threshold,args.output)
                    collect_SRY(args.bam,quality,genome,args.output)
                    print(">> Combine features into a single file\n")
                    cmd="paste "+args.output+"/XH.sorted.txt "+args.output+"/Xmap.txt "+args.output+"/Ymap.txt "+args.output+"/SRY.txt" +''' | awk 'BEGIN{OFS="\\t";print "sampleid","XH","Xmap","Ymap","XYratio","SRY"}{print $1,$2,$4,$6,$4/$6,$NF}' > '''+args.output+"/feature.txt"
                    runcmd(cmd)
                    print(">> Running sample classfication based on GMM model")
                    cmd="Rscript "+ str(Path(__file__).absolute().parent)+"/script/seGMM.r " +args.output+"/feature.txt "+str(uncertain_threshold)+" "+args.output
                    subprocess.run(cmd,shell=True)
                elif args.type == "TGS":
                    if args.chromosome=="x" and args.SRY!="True":
                        print('Beginning generate features at {T}'.format(T=time.ctime()))
                        start_time = time.time()
                        collect_XH(args.input,args.output)
                        collect_Xmap(args.bam,quality,num_threshold,args.output)
                        print(">> Combine features into a single file\n")
                        cmd="paste "+args.output+"/XH.sorted.txt "+args.output+"/Xmap.txt"+''' | awk 'BEGIN{OFS="\\t";print "sampleid","XH","Xmap"}{print $1,$2,$4}' > '''+args.output+"/feature.txt"
                        runcmd(cmd)
                        print(">> Running sample classfication based on GMM model")
                        cmd="Rscript "+ str(Path(__file__).absolute().parent)+"/script/seGMM.r " +args.output+"/feature.txt "+str(uncertain_threshold)+" "+args.output
                        subprocess.run(cmd,shell=True)
                    elif args.chromosome=="y" and args.SRY=="True":
                        print('Beginning generate features at {T}'.format(T=time.ctime()))
                        start_time = time.time()
                        collect_Ymap(args.bam,quality,num_threshold,args.output)
                        collect_SRY(args.bam,quality,genome,args.output)
                        print(">> Combine features into a single file\n")
                        cmd="paste "+args.output+"/Ymap.txt "+args.output+"/SRY.txt" +''' | awk 'BEGIN{OFS="\\t";print "sampleid","Ymap","SRY"}{print $1,$2,$NF}' >'''+args.output+"/feature.txt"
                        runcmd(cmd)
                        print(">> Running sample classfication based on GMM model")
                        cmd="Rscript "+ str(Path(__file__).absolute().parent)+"/script/seGMM.r " +args.output+"/feature.txt "+str(uncertain_threshold)+" "+args.output
                        subprocess.run(cmd,shell=True)
                    elif args.chromosome=="y" and args.SRY=="False":
                        sys.exit('Error, at least 2 features are required by GMM model')
                    elif args.chromosome=="xy" and args.SRY=="False":
                        print('Beginning generate features at {T}'.format(T=time.ctime()))
                        start_time = time.time()
                        collect_XH(args.input,args.output)
                        collect_Xmap(args.bam,quality,num_threshold,args.output)
                        collect_Ymap(args.bam,quality,num_threshold,args.output)
                        print(">> Combine features into a single file\n")
                        cmd="paste "+args.output+"/XH.sorted.txt "+args.output+"/Xmap.txt "+args.output+"/Ymap.txt "+''' | awk 'BEGIN{OFS="\\t";print "sampleid","XH","Xmap","Ymap","XYratio"}{print $1,$2,$4,$6,$4/$6}' >'''+args.output+"/feature.txt"
                        runcmd(cmd)
                        print(">> Running sample classfication based on GMM model")
                        cmd="Rscript "+ str(Path(__file__).absolute().parent)+"/script/seGMM.r " +args.output+"/feature.txt "+str(uncertain_threshold)+" "+args.output
                        subprocess.run(cmd,shell=True)
                    elif args.chromosome=="xy" and args.SRY=="True":
                        print('Beginning generate features at {T}'.format(T=time.ctime()))
                        start_time = time.time()
                        collect_XH(args.input,args.output)
                        collect_Xmap(args.bam,quality,num_threshold,args.output)
                        collect_Ymap(args.bam,quality,num_threshold,args.output)
                        collect_SRY(args.bam,quality,genome,args.output)
                        print(">> Combine features into a single file\n")
                        cmd="paste "+args.output+"/XH.sorted.txt "+args.output+"/Xmap.txt "+args.output+"/Ymap.txt "+args.output+"/SRY.txt" +''' | awk 'BEGIN{OFS="\\t";print "sampleid","XH","Xmap","Ymap","XYratio","SRY"}{print $1,$2,$4,$6,$4/$6,$NF}' >'''+args.output+"/feature.txt"
                        runcmd(cmd)
                        print(">> Running sample classfication based on GMM model")
                        cmd="Rscript "+ str(Path(__file__).absolute().parent)+"/script/seGMM.r " +args.output+"/feature.txt "+str(uncertain_threshold)+" "+args.output
                        subprocess.run(cmd,shell=True)
                    else:
                        print("Please note that you must choices the sex chromosome you want to use and tell seGMM the SRY gene is included in your reference data or not!")
                        sys.exit()
                else:
                    print("Please select the sequencing method for you data!")
                    sys.exit()
            else:
                if not os.path.isfile(args.reference):
                    print("The reference data is not exist!")
                    sys.exit()
                else:
                    ref = os.path.abspath(args.reference)
                    if args.chromosome is None:
                        order = ["XH","Xmap","Ymap","XYratio","SRY"]
                        features = ["sampleid"]
                        with open(ref) as f:
                            line=f.readline()
                            header=line.rstrip("\n").split("\t")
                        f.close()
                        if header[0]!="sampleid":
                            sys.exit("Error, the first column for reference file must is sampleid!")
                        for i in range(1,len(header)):
                            if header[i] not in order:
                                sys.exit("Error. The header of reference data is wrong. "
                                    "Please make sure the header of reference data is: sampleid, XH, Xmap, Yamp, XYratio, SRY")
                            else:
                                features.append(header[i])
                        if len(features)<=1 :
                            sys.exit(
                                "Error. At least two of features required to be "
                                "included within reference file if runnning --refenrence.")
                        else:
                            print('Beginning generate features at {T}'.format(T=time.ctime()))
                            start_time = time.time()
                            cmd = "cp " + ref + " " + args.output+"/feature.txt"
                            runcmd(cmd)
                            idx = 0
                            feature_combine = []
                            SRY = 0
                            XYratio = 0
                            for i in range(1,len(features)):
                                feature = args.output+"/"+features[i]+".txt"
                                if features[i] == "SRY":
                                    SRY = 1
                                    SRY_index = i
                                if features[i] == "XYratio":
                                    XYratio = 1
                                    XYratio_index = i
                                if not os.path.isfile(feature):
                                    with_reference(features[i], args.input,args.bam,quality,num_threshold,genome,args.output)
                                if idx==0:
                                    with open(feature) as f:
                                        line=f.readline()
                                        while line:
                                            lines=line.rstrip("\n").split("\t")
                                            feature_combine.append(lines)
                                            line = f.readline()
                                    f.close()
                                    idx=1
                                else:
                                    line_num=0
                                    with open(feature) as f:
                                        line=f.readline()
                                        while line:
                                            lines=line.rstrip("\n").split("\t")
                                            feature_combine[line_num].append(lines[1])
                                            line_num += 1
                                            line=f.readline()
                                    f.close()
                            feature_file = args.output+"/feature.txt"
                            with open(feature_file,'a+')as xh:
                                for i in range(0,len(feature_combine)):
                                    if XYratio == 1 and SRY == 1:
                                        feature_combine[i][SRY_index] = str(float(feature_combine[i][SRY_index])/float(feature_combine[i][XYratio_index]))
                                    xh.write('\t'.join(feature_combine[i])+'\n')
                            xh.close()
                        print(">> Running sample classfication based on GMM model")
                        cmd="Rscript "+ str(Path(__file__).absolute().parent)+"/script/seGMM.r " +args.output+"/feature.txt "+str(uncertain_threshold)+" "+args.output
                        subprocess.run(cmd,shell=True)
                    else:
                        print("Error, the --chromosome paremeter is not useful with an additional file!")
                        sys.exit()
        print('\nAnalysis complete for seGMM at {T}'.format(T=time.ctime()))
        time_elapsed = round(time.time()-start_time,2)
        print('Total time elapsed: {T} \n'.format(T=sec_to_str(time_elapsed)))
        print(end)
    except Exception:
        print("Error, please read the protocol of seGMM")
        raise

if __name__ == '__main__':
    main()
