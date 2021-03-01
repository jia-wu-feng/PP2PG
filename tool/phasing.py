#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
import re
import os
import sys,getopt

def getfile():
    g1_g2_snpfile=""
    samtools_th="1"
    genome1_alignment=""
    genome2_alignment=""
    genome1_output_prefix=""
    genome2_output_prefix=""
    phase_type=""
    help_str='''
phasing.py --pt=<phase_type> --g1=<genome1_alignment> --g2=<genome2_alignment> 
           --snp=<snpfile> --gop1=<genome1_output_prefix> --gop2=<genome2_output_prefix> [-h] 
    
    =====================================================
    This is a phasing tool based on two parental genomes. 
    =====================================================

    --pt=<phase_type> 
    Type: Isoseq  : Iso-Seq (PacBio Isoform Sequence) 
          RNAseq  : RNA-Seq for paired-end reads 
          BSseq   : BS-Seq (Bisulfite Sequencing) for paired-end reads
          Riboseq : Ribo-seq (Ribosome Profiling).


    --g1=<genome1_alignment>
    The bam/sam file that the reads aligned to first genome.

    The detailed methods were as follows:
    1) The type of Iso-seq reads was used minimap2 software for alignment.
     eg. minimap2 -ax splice -uf --secondary=no -C5 -O6,24 -B4 --MD ref query > results.sam

    2) The type of RNA-seq paired end reads was utilized hisat2 software for alignment.
     eg. hisat2 --dta -k 1 -x genome_index -1 query_1 -2 query_2 -S out.sam

    3) The type of BS-Seq reads was used Bismark software for alignment.
     eg. bismark --bowtie2 -genome genome_index -1 query_1 -2 query_2
     Note: Our pipeline only considers all unique mapped results in the current genomes in BS-Seq. Ambiguous reads are regarded 
           as unmaped reads.

    4) The type of Ribo-seq reads was used Tophat2(Bowtie2) software for alignment.
     eg. tophat2 -N 2 -I 50000 -G gff_file/gtf_file -o prefix genome_index query

    --g2=<genome2_alignment> 
    The bam/sam file that the reads aligned to second genome.

    --snp=<snpfile> 
    The SNP file was obtained by mummer program which used first genome as reference genome and second genome as query genome.
     eg. nucmer --mum --maxgap=500 --mincluster=100 --prefix=g1_g2 g1_genome.fa g2_genome.fa
         delta-filter -1 -q -r g1_g2.delta > g1_g2.filter.delta
         show-snps -C -H -I -T -r -l g1_g2.filter.delta > g1_g2.snp

    --st=<INT>
    The default is 1. Number of threads using samtools sort

    --gop1=<genome1_output_prefix> 
    The prefix file for first genome.

    --gop2=<genome2_output_prefix>
    The prefix file for second genome.

    -h help

    ========
    Require: 
    ========

    pysam     Version: 0.15.2
    samtools  Version: 1.6
    
    ===========
    Output file
    ===========
    Two genome-based alignment results file, containing 2 folders. One is the phasing results of g1 genome alignment, 
    which includes: g1-synteny reads, g2-synteny reads, unknown reads, and reads that only mapped to the g1 genome (g1-only reads); 
    the other is the phasing results of g2 genome alignment, which includes: g1-synteny reads, g2-synteny reads, unknown reads, 
    and reads that only mapped to the g2 genome (g2-only reads).
    g1_genome/g1reads.bam --- g1-synteny reads mapped to g1 genome
    g1_genome/g2reads.bam --- g2-synteny reads mapped to g1 genome
    g1_genome/unknowreads.bam --- unknown reads mapped to g1 genome
    g1_genome/g1onlyreads.bam --- reads that only mapped to the g1 genome (g1-only reads)
    g2_genome/g1reads.bam --- g1-synteny reads mapped to g2 genome
    g2_genome/g2reads.bam --- g2-synteny reads mapped to g2 genome
    g2_genome/unknowreads.bam --- unknown reads mapped to g2 genome
    g2_genome/g2onlyreads.bam --- reads that only mapped to the g2 genome (g2-only reads)

    g1_g2_Phasing_Report.txt --- Final Phasing Reads Report

    More information
    https://github.com/jia-wu-feng/PP2PG 
'''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"h",["g1=","g2=","snp=","gop1=","gop2=","pt=","st="])
    except getopt.GetoptError:
        print(help_str)
        sys.exit(2)
    allargv=[]
    for opt, arg in opts:
        if opt == '-h':
            print(help_str)
            sys.exit()
        elif opt == "--pt":
            if arg not in ("Isoseq","RNAseq","BSseq","Riboseq"):
                print("error:-pt")
                print(help_str)
                sys.exit()
            else:
                phase_type=arg
                allargv.append("--pt")
        elif opt =="--g1":
            genome1_alignment=arg
            allargv.append("--g1")
        elif opt =="--g2":
            genome2_alignment=arg
            allargv.append("--g2")
        elif opt =="--snp":
            g1_g2_snpfile=arg
            allargv.append("--snp")
        elif opt =="--st":
            samtools_th=arg
            allargv.append("--st")
        elif opt =="--gop1":
            genome1_output_prefix=arg
            allargv.append("--gop1")
        elif opt =="--gop2":
            genome2_output_prefix=arg
            allargv.append("--gop2")
    if "--pt" not in allargv:
        print("error: Missing parameters --pt")
        print(help_str)
        sys.exit(2)
    if "--g1" not in allargv:
        print("error: Missing parameters --g1")
        print(help_str)
        sys.exit(2)
    if "--g2" not in allargv:
        print("error: Missing parameters --g2")
        print(help_str)
        sys.exit(2)
    if "--snp" not in allargv:
        print("error: Missing parameters --snp")
        print(help_str)
        sys.exit(2)
    if "--gop1" not in allargv:
        print("error: Missing parameters --gop1")
        print(help_str)
        sys.exit(2)
    if "--gop2" not in allargv:
        print("error: Missing parameters --gop2")
        print(help_str)
        sys.exit(2)
    return genome1_alignment,genome2_alignment,g1_g2_snpfile,genome1_output_prefix,genome2_output_prefix,phase_type,samtools_th

def snpcall(r):
    allpair=r.get_aligned_pairs(with_seq=True)
    query=r.query_sequence
    snp_list=[]
    for pair in allpair:
        if pair[0]!=None and pair[1]!=None and pair[2]!=None:
            if pair[2].upper()!=query[pair[0]].upper():
                snp_list.append([pair[0],pair[1],pair[2].upper(),query[pair[0]].upper()])
    return snp_list

def makefile(gop1,gop2):
    if not os.path.exists(gop1):
        os.mkdir(gop1)
    else:
        print("Warning:"+gop1+" exit.\n")
    if not os.path.exists(gop2):
        os.mkdir(gop2)
    else:
        print("Warning:"+gop2+" exit.\n")
    return

def exitfile(filename):
    if os.path.exists(filename):
        try:
            os.remove(filename)
        except Exception as e:
            print(e)
            sys.exit(2)
        else:
            print("Warning:"+filename+" file exit, the result will be overwritten.\n")
    return

def trans_nut(a):
    if a=="T":
        b="A"
    elif a=="A":
        b="T"
    elif a=="C":
        b="G"
    elif a=="G":
        b="C"
    return b

def readsnpfile(g1_s):
    g1snp={}
    g2snp={}
    snpfile=open(g1_s,"r")
    snplines=snpfile.readlines()
    for snpline in snplines:
        op=snpline.split()
        g1snp.update({op[10]+'_'+op[0]:{"g1":op[1],"g2":op[2],"other":op[11]+'_'+op[3]}})
        if op[9]=="-1":
            g2snp.update({op[11]+'_'+op[3]:{"g1":trans_nut(op[1]),"g2":trans_nut(op[2]),"other":op[10]+'_'+op[0]}})
        else:
            g2snp.update({op[11]+'_'+op[3]:{"g1":op[1],"g2":op[2],"other":op[10]+'_'+op[0]}})
    return g1snp,g2snp

def meth_snp(zssnp,mhsnp,zs,mh):
    y=False
    if (mhsnp=="A" and zssnp=="T") or (mhsnp=="T" and zssnp=="A"):
        if mhsnp==mh and zssnp==zs:
            y=True
        else:
            y=False
    elif mhsnp=="A" and zssnp=="C":
        if mh=="A" and (zs=="C" or zs=="T"):
            y=True
        else:
            y=False
    elif mhsnp=="T" and zssnp=="G":
        if mh=="T" and (zs=="G" or zs=="A"):
            y=True
        else:
            y=False
    elif mhsnp=="C" and zssnp=="A":
        if (mh=="T" or mh=="C") and zs=="A":
            y=True
        else:
            y=False
    elif mhsnp=="C" and zssnp=="G":
        if (mh=="T" or mh=="C") and (zs=="G" or zs=="A"):
            y=True
        else:
            y=False
    elif mhsnp=="G" and zssnp=="T":
        if (mh=="G" or mh=="A") and zs=="T":
            y=True
        else:
            y=False
    elif mhsnp=="G" and zssnp=="C":
        if (mh=="G" or mh=="A") and (zs=="T" or zs=="C"):
            y=True
        else:
            y=False
    else:
        y=False
    return y

def isoreadsphase(samlist,g1snp,g2snp,gop1,gop2):
    reads_type=""
    g1tag=0
    g2tag=0
    unknowtag=0
    supportfile={"g1tag_support":["g1"],"g1tag_unknow_support":["g1_un"],"g2tag_support":["g2"],"g2tag_unknow_support":["g2_un"]}
    if len(samlist["g1"])==1 and len(samlist["g2"])==1:
        if samlist["g1"][0].flag==4 and samlist["g2"][0].flag==4:
            reads_type="discarded"
        elif samlist["g1"][0].flag==4 and samlist["g2"][0].flag!=4:
            reads_type="g2only"
        elif samlist["g1"][0].flag!=4 and samlist["g2"][0].flag==4:
            reads_type="g1only"
        else:
            snp_po=snpcall(samlist["g1"][0])
            for snp in snp_po:
                if samlist["g1"][0].reference_name+"_"+str(snp[1]+1) in g1snp:
                    snplist=g1snp[samlist["g1"][0].reference_name+"_"+str(snp[1]+1)]
                    if snp[3]==snplist["g2"]:
                        g2tag+=1
                        supportfile["g2tag_support"].append(gop1+"_"+samlist["g1"][0].reference_name+"_"+str(int(snp[1])+1))
                    else:
                        unknowtag+=1
                        supportfile["g2tag_unknow_support"].append(gop1+"_"+samlist["g1"][0].reference_name+"_"+str(int(snp[1])+1))
            snp_po=snpcall(samlist["g2"][0])
            for snp in snp_po:
                if samlist["g2"][0].reference_name+"_"+str(snp[1]+1) in g2snp:
                    snplist=g2snp[samlist["g2"][0].reference_name+"_"+str(snp[1]+1)]
                    if snp[3]==snplist["g1"]:
                        g1tag+=1
                        supportfile["g1tag_support"].append(gop1+"_"+samlist["g2"][0].reference_name+"_"+str(int(snp[1])+1))
                    else:
                        unknowtag+=1
                        supportfile["g1tag_unknow_support"].append(gop1+"_"+samlist["g2"][0].reference_name+"_"+str(int(snp[1])+1))
            if g1tag > g2tag :
                reads_type="g1"
            elif g1tag < g2tag:
                reads_type="g2"
            else:
                reads_type="unk"
    elif len(samlist["g1"])>1 and len(samlist["g2"])>1:
        reads_type="unk"
    elif len(samlist["g1"])==0 and len(samlist["g2"])>=1:
        reads_type="g2only"
    elif len(samlist["g1"])>=1 and len(samlist["g2"])==0:
        reads_type="g1only"
    else:
        reads_type="unk"
    return reads_type,supportfile


def rnapairreadsphase(samlist,g1snp,g2snp,gop1,gop2):
    reads_type=""
    g1tag=0
    g2tag=0
    unknowtag=0
    supportfile={"g1tag_support":["g1"],"g1tag_unknow_support":["g1_un"],"g2tag_support":["g2"],"g2tag_unknow_support":["g2_un"]}
    if len(samlist["g1_mate1"])==1 and len(samlist["g1_mate2"])==1 and len(samlist["g2_mate1"])==1 and len(samlist["g2_mate2"])==1:
        if (not samlist["g1_mate1"][0].is_unmapped) and ( not samlist["g1_mate2"][0].is_unmapped) and samlist["g2_mate1"][0].is_unmapped and samlist["g2_mate2"][0].is_unmapped:
            reads_type="g1only"
        elif samlist["g1_mate1"][0].is_unmapped and samlist["g1_mate2"][0].is_unmapped and (not samlist["g2_mate1"][0].is_unmapped) and (not samlist["g2_mate2"][0].is_unmapped):
            reads_type="g2only"
        elif samlist["g1_mate1"][0].is_unmapped and samlist["g1_mate2"][0].is_unmapped and samlist["g2_mate1"][0].is_unmapped and samlist["g2_mate2"][0].is_unmapped:
            reads_type="discarded"
        else:
            if not samlist["g1_mate1"][0].is_unmapped:
                snp_po=snpcall(samlist["g1_mate1"][0])
                for snp in snp_po:
                    if samlist["g1_mate1"][0].reference_name+"_"+str(snp[1]+1) in g1snp:
                        snplist=g1snp[samlist["g1_mate1"][0].reference_name+"_"+str(snp[1]+1)]
                        if snp[3]==snplist["g2"]:
                            g2tag+=1
                            supportfile["g2tag_support"].append(gop1+"_"+samlist["g1_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g2tag_unknow_support"].append(gop1+"_"+samlist["g1_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
            if not samlist["g1_mate2"][0].is_unmapped:
                snp_po=snpcall(samlist["g1_mate2"][0])
                for snp in snp_po:
                    if samlist["g1_mate2"][0].reference_name+"_"+str(snp[1]+1) in g1snp:
                        snplist=g1snp[samlist["g1_mate2"][0].reference_name+"_"+str(snp[1]+1)]
                        if snp[3]==snplist["g2"]:
                            g2tag+=1
                            supportfile["g2tag_support"].append(gop1+"_"+samlist["g1_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g2tag_unknow_support"].append(gop1+"_"+samlist["g1_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
            if not samlist["g2_mate1"][0].is_unmapped:
                snp_po=snpcall(samlist["g2_mate1"][0])
                for snp in snp_po:
                    if samlist["g2_mate1"][0].reference_name+"_"+str(snp[1]+1) in g2snp:
                        snplist=g2snp[samlist["g2_mate1"][0].reference_name+"_"+str(snp[1]+1)]
                        if snp[3]==snplist["g1"]:
                            g1tag+=1
                            supportfile["g1tag_support"].append(gop1+"_"+samlist["g2_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g1tag_unknow_support"].append(gop1+"_"+samlist["g2_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
            if not samlist["g2_mate2"][0].is_unmapped:
                snp_po=snpcall(samlist["g2_mate2"][0])
                for snp in snp_po:
                    if samlist["g2_mate2"][0].reference_name+"_"+str(snp[1]+1) in g2snp:
                        snplist=g2snp[samlist["g2_mate2"][0].reference_name+"_"+str(snp[1]+1)]
                        if snp[3]==snplist["g1"]:
                            g1tag+=1
                            supportfile["g1tag_support"].append(gop1+"_"+samlist["g2_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g1tag_unknow_support"].append(gop1+"_"+samlist["g2_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
            if g1tag > g2tag :
                reads_type="g1"
            elif g1tag < g2tag:
                reads_type="g2"
            else:
                reads_type="unk"
    elif len(samlist["g1_mate1"])==0 and len(samlist["g1_mate2"])==0 and len(samlist["g2_mate1"])>=1 and len(samlist["g2_mate2"])>=1:
        reads_type="g2only"
    elif len(samlist["g1_mate1"])>=1 and len(samlist["g1_mate2"])>=1 and len(samlist["g2_mate1"])==0 and len(samlist["g2_mate2"])==0:
        reads_type="g1only"
    elif len(samlist["g1_mate1"])>=1 and len(samlist["g1_mate2"])==0 and len(samlist["g2_mate1"])>=1 and len(samlist["g2_mate2"])==0:
        reads_type="single"
    elif len(samlist["g1_mate1"])==0 and len(samlist["g1_mate2"])>=1 and len(samlist["g2_mate1"])==0 and len(samlist["g2_mate2"])>=1:
        reads_type="single"
    else:
        reads_type="unk"
    return reads_type,supportfile


def wgbsreadsphase(samlist,g1snp,g2snp,gop1,gop2):
    reads_type=""
    g1tag=0
    g2tag=0
    unknowtag=0
    supportfile={"g1tag_support":["g1"],"g1tag_unknow_support":["g1_un"],"g2tag_support":["g2"],"g2tag_unknow_support":["g2_un"]}
    if len(samlist["g1_mate1"])==1 and len(samlist["g1_mate2"])==1 and len(samlist["g2_mate1"])==1 and len(samlist["g2_mate2"])==1:
        if samlist["g1_mate1"][0].get_tag("XM",with_value_type = True)[1]=="C" or \
        samlist["g1_mate2"][0].get_tag("XM",with_value_type = True)[1]=="C" or \
        samlist["g2_mate1"][0].get_tag("XM",with_value_type = True)[1]=="C" or\
        samlist["g2_mate2"][0].get_tag("XM",with_value_type = True)[1]=="C":
            reads_type="discarded"
        elif (not samlist["g1_mate1"][0].is_unmapped) and (not samlist["g1_mate2"][0].is_unmapped) and samlist["g2_mate1"][0].is_unmapped and samlist["g2_mate2"][0].is_unmapped:
            reads_type="g1only"
        elif samlist["g1_mate1"][0].is_unmapped and samlist["g1_mate2"][0].is_unmapped and (not samlist["g2_mate1"][0].is_unmapped ) and ( not samlist["g2_mate2"][0].is_unmapped) :
            reads_type="g2only"
        elif samlist["g1_mate1"][0].is_unmapped and samlist["g1_mate2"][0].is_unmapped and samlist["g2_mate1"][0].is_unmapped and samlist["g2_mate2"][0].is_unmapped:
            reads_type="discarded"
        else:
            if not samlist["g1_mate1"][0].is_unmapped:
                snp_po=snpcall(samlist["g1_mate1"][0])
                for snp in snp_po:
                    if samlist["g1_mate1"][0].reference_name+"_"+str(snp[1]+1) in g1snp:
                        snplist=g1snp[samlist["g1_mate1"][0].reference_name+"_"+str(snp[1]+1)]
                        if meth_snp(snplist["g1"],snplist["g2"],snp[2],snp[3]):
                            g2tag+=1
                            supportfile["g2tag_support"].append(gop1+"_"+samlist["g1_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g2tag_unknow_support"].append(gop1+"_"+samlist["g1_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
            if not samlist["g1_mate2"][0].is_unmapped:
                snp_po=snpcall(samlist["g1_mate2"][0])
                for snp in snp_po:
                    if samlist["g1_mate2"][0].reference_name+"_"+str(snp[1]+1) in g1snp:
                        snplist=g1snp[samlist["g1_mate2"][0].reference_name+"_"+str(snp[1]+1)]
                        if meth_snp(snplist["g1"],snplist["g2"],snp[2],snp[3]):
                            g2tag+=1
                            supportfile["g2tag_support"].append(gop1+"_"+samlist["g1_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g2tag_unknow_support"].append(gop1+"_"+samlist["g1_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
            if not samlist["g2_mate1"][0].is_unmapped:
                snp_po=snpcall(samlist["g2_mate1"][0])
                for snp in snp_po:
                    if samlist["g2_mate1"][0].reference_name+"_"+str(snp[1]+1) in g2snp:
                        snplist=g2snp[samlist["g2_mate1"][0].reference_name+"_"+str(snp[1]+1)]
                        if meth_snp(snplist["g1"],snplist["g2"],snp[3],snp[2]):
                            g1tag+=1
                            supportfile["g1tag_support"].append(gop1+"_"+samlist["g2_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g1tag_unknow_support"].append(gop1+"_"+samlist["g2_mate1"][0].reference_name+"_"+str(int(snp[1])+1))
            if not samlist["g2_mate2"][0].is_unmapped:
                snp_po=snpcall(samlist["g2_mate2"][0])
                for snp in snp_po:
                    if samlist["g2_mate2"][0].reference_name+"_"+str(snp[1]+1) in g2snp:
                        snplist=g2snp[samlist["g2_mate2"][0].reference_name+"_"+str(snp[1]+1)]
                        if meth_snp(snplist["g1"],snplist["g2"],snp[3],snp[2]):
                            g1tag+=1
                            supportfile["g1tag_support"].append(gop1+"_"+samlist["g2_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
                        else:
                            unknowtag+=1
                            supportfile["g1tag_unknow_support"].append(gop1+"_"+samlist["g2_mate2"][0].reference_name+"_"+str(int(snp[1])+1))
            if g1tag > g2tag :
                reads_type="g1"
            elif g1tag < g2tag:
                reads_type="g2"
            else:
                reads_type="unk"
    elif len(samlist["g1_mate1"])==0 and len(samlist["g1_mate2"])==0 and len(samlist["g2_mate1"])>=1 and len(samlist["g2_mate2"])>=1:
        reads_type="g2only"
    elif len(samlist["g1_mate1"])>=1 and len(samlist["g1_mate2"])>=1 and len(samlist["g2_mate1"])==0 and len(samlist["g2_mate2"])==0:
        reads_type="g1only"
    elif len(samlist["g1_mate1"])>=1 and len(samlist["g1_mate2"])==0 and len(samlist["g2_mate1"])>=1 and len(samlist["g2_mate2"])==0:
        reads_type="single"
    elif len(samlist["g1_mate1"])==0 and len(samlist["g1_mate2"])>=1 and len(samlist["g2_mate1"])==0 and len(samlist["g2_mate2"])>=1:
        reads_type="single"
    else:
        reads_type="unk"
    return reads_type,supportfile


def readswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport):
    if readt=="unk":
        for rp in samlist["g1"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_unknowreadsfile.write(rp)
            snpsupport.write(rp.query_name+"\t"+str(",".join(support["g1tag_support"]))+"\t"+"\t"+str(",".join(support["g1tag_unknow_support"]))+"\t"+\
                             str(",".join(support["g2tag_support"]))+"\t"+str(",".join(support["g2tag_unknow_support"]))+"\n")
        for rp in samlist["g2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_unknowreadsfile.write(rp)
    elif readt=="g1":
        for rp in samlist["g1"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_g1readsfile.write(rp)
            snpsupport.write(rp.query_name+"\t"+str(",".join(support["g1tag_support"]))+"\t"+"\t"+str(",".join(support["g1tag_unknow_support"]))+"\t"+\
                            str(",".join(support["g2tag_support"]))+"\t"+str(",".join(support["g2tag_unknow_support"]))+"\n")
        for rp in samlist["g2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_g1readsfile.write(rp)
    elif readt=="g2":
        for rp in samlist["g1"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_g2readsfile.write(rp)
            snpsupport.write(rp.query_name+"\t"+str(",".join(support["g1tag_support"]))+"\t"+"\t"+str(",".join(support["g1tag_unknow_support"]))+"\t"+\
                            str(",".join(support["g2tag_support"]))+"\t"+str(",".join(support["g2tag_unknow_support"]))+"\n")
        for rp in samlist["g2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_g2readsfile.write(rp)
    elif readt=="g1only":
        for rp in samlist["g1"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_g1onlyreadsfile.write(rp)
    elif readt=="g2only":
        for rp in samlist["g2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_g2onlyreadsfile.write(rp)
    return

def pairsreadswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport):
    if readt=="unk" or readt=="single":
        for rp in samlist["g1_mate1"]+samlist["g1_mate2"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_unknowreadsfile.write(rp)
        for rp in samlist["g2_mate1"]+samlist["g2_mate2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_unknowreadsfile.write(rp)
    elif readt=="g1":
        for rp in samlist["g1_mate1"]+samlist["g1_mate2"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_g1readsfile.write(rp)
            snpsupport.write(rp.query_name+"\t"+str(",".join(support["g1tag_support"]))+"\t"+"\t"+str(",".join(support["g1tag_unknow_support"]))+"\t"+\
                            str(",".join(support["g2tag_support"]))+"\t"+str(",".join(support["g2tag_unknow_support"]))+"\n")
        for rp in samlist["g2_mate1"]+samlist["g2_mate2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_g1readsfile.write(rp)
    elif readt=="g2":
        for rp in samlist["g1_mate1"]+samlist["g1_mate2"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_g2readsfile.write(rp)
            snpsupport.write(rp.query_name+"\t"+str(",".join(support["g1tag_support"]))+"\t"+"\t"+str(",".join(support["g1tag_unknow_support"]))+"\t"+\
                            str(",".join(support["g2tag_support"]))+"\t"+str(",".join(support["g2tag_unknow_support"]))+"\n")
        for rp in samlist["g2_mate1"]+samlist["g2_mate2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_g2readsfile.write(rp)
    elif readt=="g1only":
        for rp in samlist["g1_mate1"]+samlist["g1_mate2"]:
            rp.query_name=rp.query_name.replace("_g1","")
            g1_g1onlyreadsfile.write(rp)
    elif readt=="g2only":
        for rp in samlist["g2_mate1"]+samlist["g2_mate2"]:
            rp.query_name=rp.query_name.replace("_g2","")
            g2_g2onlyreadsfile.write(rp)
    return

def riboreadsphase(samlist,g1snp,g2snp,gop1,gop2):
    reads_type=""
    g1tag=0
    g2tag=0
    unknowtag=0
    supportfile={"g1tag_support":["g1"],"g1tag_unknow_support":["g1_un"],"g2tag_support":["g2"],"g2tag_unknow_support":["g2_un"]}
    if len(samlist["g1"])==1 and len(samlist["g2"])==1:
        if samlist["g1"][0].is_unmapped and samlist["g2"][0].is_unmapped:
            reads_type="discarded"
        elif samlist["g1"][0].is_unmapped and (not samlist["g2"][0].is_unmapped):
            reads_type="g2only"
        elif (not samlist["g1"][0].is_unmapped ) and samlist["g2"][0].is_unmapped:
            reads_type="g1only"
        else:
            snp_po=snpcall(samlist["g1"][0])
            for snp in snp_po:
                if samlist["g1"][0].reference_name+"_"+str(snp[1]+1) in g1snp:
                    snplist=g1snp[samlist["g1"][0].reference_name+"_"+str(snp[1]+1)]
                    if snp[3]==snplist["g2"]:
                        g2tag+=1
                        supportfile["g2tag_support"].append(gop1+"_"+samlist["g1"][0].reference_name+"_"+str(int(snp[1])+1))
                    else:
                        unknowtag+=1
                        supportfile["g2tag_unknow_support"].append(gop1+"_"+samlist["g1"][0].reference_name+"_"+str(int(snp[1])+1))
            snp_po=snpcall(samlist["g2"][0])
            for snp in snp_po:
                if samlist["g2"][0].reference_name+"_"+str(snp[1]+1) in g2snp:
                    snplist=g2snp[samlist["g2"][0].reference_name+"_"+str(snp[1]+1)]
                    if snp[3]==snplist["g1"]:
                        g1tag+=1
                        supportfile["g1tag_support"].append(gop1+"_"+samlist["g2"][0].reference_name+"_"+str(int(snp[1])+1))
                    else:
                        unknowtag+=1
                        supportfile["g1tag_unknow_support"].append(gop1+"_"+samlist["g2"][0].reference_name+"_"+str(int(snp[1])+1))
            if g1tag > g2tag :
                reads_type="g1"
            elif g1tag < g2tag:
                reads_type="g2"
            else:
                reads_type="unk"
    elif len(samlist["g1"])>1 and len(samlist["g2"])>1:
        reads_type="unk"
    elif len(samlist["g1"])==0 and len(samlist["g2"])>=1:
        reads_type="g2only"
    elif len(samlist["g1"])>=1 and len(samlist["g2"])==0:
        reads_type="g1only"
    else:
        reads_type="unk"
    return reads_type,supportfile


def main():
    g1_al,g2_al,g1_s,gop1,gop2,pt,sth=getfile()
    makefile(gop1+"_genome",gop2+"_genome")
    g1snp,g2snp=readsnpfile(g1_s)

    g1samfile = pysam.AlignmentFile(g1_al,"r")
    g2samfile = pysam.AlignmentFile(g2_al,"r")

    exitfile("tmp1.bam")
    exitfile("tmp2.bam")
    exitfile("tmp3.bam")
    exitfile("tmp_sortname.bam")

    tmp1=pysam.AlignmentFile("tmp1.bam", "wb",template=g1samfile,threads=int(sth))

    for r in g1samfile:
        r.query_name=r.query_name+"_g1"
        tmp1.write(r)
    tmp1.close()

    tmp2=pysam.AlignmentFile("tmp2.bam", "wb",template=g2samfile,threads=int(sth))
    for r in g2samfile:
        r.query_name=r.query_name+"_g2"
        tmp2.write(r)
    tmp2.close()

    pysam.merge("-@",sth,"tmp3.bam","tmp1.bam","tmp2.bam")
    pysam.sort("-@",sth,"-n","-o","tmp_sortname.bam","tmp3.bam")

    tmp_all=pysam.AlignmentFile("tmp_sortname.bam","rb")

    tmp_bamheader=tmp_all.text
    g1_header=g1samfile.text
    g2_header=g2samfile.text

    snpsupport = open(gop1+"_"+gop2+"_"+"support.txt","w")

    g1_g1readsfile = pysam.AlignmentFile(gop1+"_genome/"+gop1+"reads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    g1_g2readsfile = pysam.AlignmentFile(gop1+"_genome/"+gop2+"reads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    g1_unknowreadsfile = pysam.AlignmentFile(gop1+"_genome/"+"unknownreads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    g1_g1onlyreadsfile = pysam.AlignmentFile(gop1+"_genome/"+gop1+"onlyreads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)

    g2_g1readsfile = pysam.AlignmentFile(gop2+"_genome/"+gop1+"reads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    g2_g2readsfile = pysam.AlignmentFile(gop2+"_genome/"+gop2+"reads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    g2_unknowreadsfile = pysam.AlignmentFile(gop2+"_genome/"+"unknownreads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    g2_g2onlyreadsfile = pysam.AlignmentFile(gop2+"_genome/"+gop2+"onlyreads.sam", "w", template=tmp_all,threads=int(sth),add_sam_header=False)
    
    stat_number={"g1":0,"g2":0,"g1only":0,"g2only":0,"unk":0,"single":0,"discarded":0}

    if pt=="Isoseq":
        r0_name=""
        samlist={"g1":[],"g2":[]}
        for r in tmp_all:
            if r0_name!=r.query_name.replace("_g1","").replace("_g2","") and r0_name!="":
                readt,support=isoreadsphase(samlist,g1snp,g2snp,gop1,gop2)
                readswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
                stat_number[readt]+=1
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
                samlist={"g1":[],"g2":[]}
            elif r0_name=="":
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
            if r.query_name.find("_g1")!=-1:
                samlist["g1"].append(r)
            elif r.query_name.find("_g2")!=-1:
                samlist["g2"].append(r)
        readt,support=isoreadsphase(samlist,g1snp,g2snp,gop1,gop2)
        stat_number[readt]+=1
        readswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
    elif pt=="RNAseq":
        r0_name=""
        samlist={"g1_mate1":[],"g1_mate2":[],"g2_mate1":[],"g2_mate2":[]}
        for r in tmp_all:
            if r0_name!=r.query_name.replace("_g1","").replace("_g2","") and r0_name!="":
                readt,support=rnapairreadsphase(samlist,g1snp,g2snp,gop1,gop2)
                if readt=="single":
                    stat_number[readt]+=1
                else:
                    stat_number[readt]+=2
                pairsreadswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
                samlist={"g1_mate1":[],"g1_mate2":[],"g2_mate1":[],"g2_mate2":[]}
            elif r0_name=="":
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
            if r.query_name.find("_g1")!=-1:
                if r.is_read1:
                    samlist["g1_mate1"].append(r)
                elif r.is_read2:
                    samlist["g1_mate2"].append(r)
            elif r.query_name.find("_g2")!=-1:
                if r.is_read1:
                    samlist["g2_mate1"].append(r)
                elif r.is_read2:
                    samlist["g2_mate2"].append(r)
        readt,support=rnapairreadsphase(samlist,g1snp,g2snp,gop1,gop2)
        if readt=="single":
            stat_number[readt]+=1
        else:
            stat_number[readt]+=2
        pairsreadswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
    elif pt=="BSseq":
        r0_name=""
        samlist={"g1_mate1":[],"g1_mate2":[],"g2_mate1":[],"g2_mate2":[]}
        for r in tmp_all:
            if r0_name!=r.query_name.replace("_g1","").replace("_g2","") and r0_name!="":
                readt,support=wgbsreadsphase(samlist,g1snp,g2snp,gop1,gop2)
                if readt=="single":
                    stat_number[readt]+=1
                else:
                    stat_number[readt]+=2
                pairsreadswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
                samlist={"g1_mate1":[],"g1_mate2":[],"g2_mate1":[],"g2_mate2":[]}
            elif r0_name=="":
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
            if r.query_name.find("_g1")!=-1:
                if r.is_read1:
                    samlist["g1_mate1"].append(r)
                elif r.is_read2:
                    samlist["g1_mate2"].append(r)
            elif r.query_name.find("_g2")!=-1:
                if r.is_read1:
                    samlist["g2_mate1"].append(r)
                elif r.is_read2:
                    samlist["g2_mate2"].append(r)
        readt,support=wgbsreadsphase(samlist,g1snp,g2snp,gop1,gop2)
        if readt=="single":
            stat_number[readt]+=1
        else:
            stat_number[readt]+=2
        pairsreadswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
    elif pt=="Riboseq":
        r0_name=""
        samlist={"g1":[],"g2":[]}
        for r in tmp_all:
            if r0_name!=r.query_name.replace("_g1","").replace("_g2","") and r0_name!="":
                readt,support=riboreadsphase(samlist,g1snp,g2snp,gop1,gop2)
                readswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)
                stat_number[readt]+=1
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
                samlist={"g1":[],"g2":[]}
            elif r0_name=="":
                r0_name=r.query_name.replace("_g1","").replace("_g2","")
            if r.query_name.find("_g1")!=-1:
                samlist["g1"].append(r)
            elif r.query_name.find("_g2")!=-1:
                samlist["g2"].append(r)
        readt,support=riboreadsphase(samlist,g1snp,g2snp,gop1,gop2)
        stat_number[readt]+=1
        readswrite(samlist,readt,support,g1_g1readsfile,g1_g2readsfile,g1_unknowreadsfile,g1_g1onlyreadsfile,g2_g1readsfile,g2_g2readsfile,g2_unknowreadsfile,g2_g2onlyreadsfile,snpsupport)

    g1_g1readsfile.close()
    g1_g2readsfile.close()
    g1_unknowreadsfile.close()
    g1_g1onlyreadsfile.close()
    g2_g1readsfile.close()
    g2_g2readsfile.close()
    g2_unknowreadsfile.close()
    g2_g2onlyreadsfile.close()
    snpsupport.close()

    g1file=[gop1+"_genome/"+gop1+"reads",gop1+"_genome/"+gop2+"reads",gop1+"_genome/unknownreads",gop1+"_genome/"+gop1+"onlyreads"]
    g2file=[gop2+"_genome/"+gop1+"reads",gop2+"_genome/"+gop2+"reads",gop2+"_genome/unknownreads",gop2+"_genome/"+gop2+"onlyreads"]
    PG='@PG\tID:phasing.py\tPN:PP2PG\tVN:1.0.0\tCL:'+' '.join(sys.argv)+'\n'

    tmp4=open("tmp4.header","w")
    tmp4.write(g1_header+PG)
    tmp4.close()
    tmp5=open("tmp5.header","w")
    tmp5.write(g2_header+PG)
    tmp5.close()

    for f in g1file:
        exitfile(f+".bam")
        exitfile(f+".tmp.sam")
        try:
            os.system("cat tmp4.header "+f+".sam >"+f+".tmp.sam")
            print("Excuting: samtools view -@ "+sth+" -b -S "+f+".tmp.sam -o "+f+".bam\n")
            os.system("samtools view -@ "+sth+" -b -S "+f+".tmp.sam -o "+f+".bam")
            os.system("rm "+f+".sam "+f+".tmp.sam")
        except Exception as e:
            print(e)
            sys.exit(2)
        else:
            print("*** "+f+".sam converted bam file successfully. ***\n")
    for f in g2file:
        exitfile(f+".bam")
        exitfile(f+".tmp.sam")
        try:
            os.system("cat tmp5.header "+f+".sam >"+f+".tmp.sam")
            print("Excuting: samtools view -@ "+sth+" -b -S "+f+".tmp.sam -o "+f+".bam\n")
            os.system("samtools view -@ "+sth+" -b -S "+f+".tmp.sam -o "+f+".bam")
            os.system("rm "+f+".sam "+f+".tmp.sam")
        except Exception as e:
            print(e)
            sys.exit(2)
        else:
            print("*** "+f+".sam converted bam file successfully. ***\n")
    try:
        os.remove("tmp1.bam")
        os.remove("tmp2.bam")
        os.remove("tmp3.bam")
        os.remove("tmp_sortname.bam")
        os.remove("tmp4.header")
        os.remove("tmp5.header")
        os.remove(gop1+"_"+gop2+"_support.txt")
    except Exception as e:
        print(e)
        sys.exit(2)
    else:
        print("*** Deleted tmp file successfully. ***\n")

    g1samfile.close()
    g2samfile.close()
    tmp_all.close()
    g1reads_stats=stat_number["g1"]
    g2reads_stats=stat_number["g2"]
    unknowreads_stats=stat_number["unk"]+stat_number["single"]
    g1onlyreads_stats=stat_number["g1only"]
    g2onlyreads_stats=stat_number["g2only"]
    datatype={"Isoseq":"Iso-Seq (PacBio Isoform Sequence).",\
              "RNAseq":"RNA-Seq for paired-end reads.",\
              "BSseq":"BS-Seq (Bisulfite Sequencing) for paired-end reads." ,\
              "Riboseq":"Ribo-seq (Ribosome Profiling)."}
    head="""
==========================
Final Phasing Reads Report
==========================

"""+"Phasing Data Type:\t"+datatype[pt]+"""

Note: When calculating the Separation Rate, the reads which are unmapped in two parental genomes are discarded.

Type of reads:\tRead Counts
"""
    report=gop1+"-synteny reads:\t"+str(g1reads_stats)+"\n"+\
           gop2+"-synteny reads:\t "+str(g2reads_stats)+"\n"+\
           "Unknown reads:\t"+str(unknowreads_stats)+"\n"+\
           gop1+"-only reads:\t"+str(g1onlyreads_stats)+"\n"+\
           gop2+"-only reads:\t"+str(g2onlyreads_stats)+"\n\n"+\
           "Separation Rate:\t"+str(round((float(g1reads_stats)+float(g2reads_stats))/(float(g1reads_stats)+float(g2reads_stats)+float(g1onlyreads_stats)+float(g2onlyreads_stats)+float(unknowreads_stats))*100,2))+"%\n"
    file_report=open(gop1+"_"+gop2+"_Phasing_Report.txt","w")
    file_report.write(head+report)
    file_report.close()
    print(head+report)
    print("""
=================
Phasing Finished!
=================
""")

if __name__ == '__main__':
    main()
