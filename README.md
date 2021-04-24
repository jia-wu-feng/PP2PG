# PP2PG

PP2PG(Phasing pipeline based on two parental genomes)


Require: 
  
  1) MUMmer     Version:4.0.0beta2 (https://github.com/mummer4/mummer)
  2) Bismark    Version:v0.22.1    (https://github.com/FelixKrueger/Bismark)
  3) minimap2   Version:2.17       (https://github.com/lh3/minimap2)
  4) HISAT2     Version:2.1.0      (http://ccb.jhu.edu/software/hisat2/manual.shtml)
  5) pysam      Version:0.15.2     (https://pypi.org/project/pysam/)     
  6) SAMtools   Version:1.6        (http://samtools.sourceforge.net/)  

Installation:

       Note: Please confirm the required software before installation
       git clone git@github.com:jia-wu-feng/PP2PG.git
       cd PP2PG
       chmod +X tool/phasing.py 

![image1](https://github.com/jia-wu-feng/PP2PG/blob/master/img/img1.png)

Fig.1 The phasing pipeline based on two parental genomes (PP2PG). It mainly consists of four steps: (i) Mapping, (ii) SNP Calling, (iii) Phasing and (iv) Classification.

1.Mapping:

1) The type of Iso-Seq reads was used minimap2 software for alignment:

       minimap2 -ax splice -uf --secondary=no -C5 -O6,24 -B4 --MD ref query > results.sam      

2) The type of RNA-Seq paired end reads was utilized hisat2 software for alignment:

       hisat2 --dta -k 1 -x genome_index -1 query_1 -2 query_2 -S out.sam

3) The type of BS-Seq reads was used Bismark software for alignment:
      
       bismark --bowtie2 -genome genome_index -1 query_1 -2 query_2
    Note: Our pipeline only considers all unique mapped results in the current genomes in BS-Seq. Ambiguous reads are regarded 
    as unmaped reads. 

4) The type of Ribo-seq reads was used Tophat2(Bowtie2) software for alignment:

       tophat2 -N 2 -I 50000 -G gff_file/gtf_file -o prefix genome_index query


2.SNP Calling:

    nucmer --mum --maxgap=500 --mincluster=100 --prefix=g1_g2 g1_genome.fa g1_genome.fa
    
    delta-filter -1 -q -r g1_g2.delta > g1_g2.filter.delta
    
    show-snps -C -H -I -T -r -l g1_g2.filter.delta > g1_g2.snp


3.Scoring and Classification:

    phasing.py --pt=<phase_type> --g1=<genome1_alignment> --g2=<genome2_alignment> 
               --snp=<snpfile> --gop1=<genome1_output_prefix> --gop2=<genome2_output_prefix>

4.Output file:

Two genome-based alignment results file, containing 2 folders. One is the phasing results of g1 genome alignment, which includes: g1-synteny reads, g2-synteny reads, unknown reads, and reads that only mapped to the g1 genome (g1-only reads); the other is the phasing results of g2 genome alignment, which includes: g1-synteny reads, g2-synteny reads, unknown reads, and reads that only mapped to the g2 genome (g2-only reads).

    g1_genome/g1reads.bam --- g1-synteny reads mapped to g1 genome
    g1_genome/g2reads.bam --- g2-synteny reads mapped to g1 genome
    g1_genome/unknowreads.bam --- unknown reads mapped to g1 genome
    g1_genome/g1onlyreads.bam --- reads that only mapped to the g1 genome (g1-only reads)
    
    g2_genome/g1reads.bam --- g1-synteny reads mapped to g2 genome
    g2_genome/g2reads.bam --- g2-synteny reads mapped to g2 genome
    g2_genome/unknowreads.bam --- unknown reads mapped to g2 genome
    g2_genome/g2onlyreads.bam --- reads that only mapped to the g2 genome (g2-only reads)


g1_g2_Phasing_Report.txt --- Final Phasing Reads Report

Example:

    ==========================
    Final Phasing Reads Report
    ==========================
    
    Phasing Data Type:      Iso-Seq (PacBio Isoform Sequence).
    
    Note: When calculating the Separation Rate, the reads which are unmapped in two parental genomes are discarded.
    
    Type of reads:  Read Counts
    MH63-synteny reads:	129
    ZS97-synteny reads:	129
    Unknown reads:	47
    MH63-only reads:	2
    ZS97-only reads:	3
    
    Separation Rate:        83.23%

5.More analysis:

Iso-Seq:

Identification of alternative splicing

(1) Samtools was used to convert the phasing bam files into sequence files. 

(2) Then alignPacBio.py (TAPIS) was applied for correction, and run_tapis.py was applied for assembling the corrected transcripts (Abdel-Ghany et al. 2016). 

(3) SUPPA2 was utilized to identify the alternative splicing (AS) from the transcripts which corrected and assembled by TAPIS (Trincado et al. 2018). 

Identification of allele-specific expressed (ASE) genes

(1) HTSeq (Anders, Pyl, and Huber 2015) was employed for counting the number of full-length transcripts. 

(2) To identify the gene bias, binom.test in R language was applied for testing whether the ratio of alleles deviated from 0.5 ( total reads≥5,qvalue≤0.05 ).

RNA-Seq:

Identification of ASE genes and trans-regulation

(1) Samtools was used to sort bam files by name.

(2) HTSeq (Anders, Pyl, and Huber 2015) was used to count reads for every sample after phasing. 

(3) DESeq2 was performed to identify the differential expression (reads ≥5, qvalue ≤0.05 and |log2FoldChange| > 1).

Construction the allele co-expression network of SY63

(1) StringTie (Pertea et al. 2015) was utilized to calculate the expression of two parental alleles in the same genome. 

(2) The two files of allele labels were modified and merged these two files into one file. 

(3) An allele co-expression network was constructed based on Weighted Gene Co-expression Network Analysis (WGCNA) (Langfelder and Horvath 2008).

BS-Seq

Identification of allele-specific DNA methylation (ASM)

(1) samtools was used to sort bam files by name after phasing. 

(2) The bismark_methylation_extractor was utilized to recall methylations based on the MH63 genome with the parameters (-p --gzip --bedGraph --comprehensive --cytosine_report -CX) (Krueger and Andrews 2011). 

(3) Cytosine report can be used to obtain each C methylation and unmethylation situation. Users can calculate the differential methylation sites and regions.


Reference

Abdel-Ghany, S. E., M. Hamilton, J. L. Jacobi, P. Ngam, N. Devitt, F. Schilkey, A. Ben-Hur, and A. S. Reddy. 2016. 'A survey of the sorghum transcriptome using single-molecule long reads', Nat Commun, 7: 11706.

Anders, S., P. T. Pyl, and W. Huber. 2015. 'HTSeq--a Python framework to work with high-throughput sequencing data', Bioinformatics, 31: 166-9.

Krueger, F., and S. R. Andrews. 2011. 'Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications', Bioinformatics, 27: 1571-2.

Langfelder, P., and S. Horvath. 2008. 'WGCNA: an R package for weighted correlation network analysis', BMC Bioinformatics, 9: 559.

Pertea, M., G. M. Pertea, C. M. Antonescu, T. C. Chang, J. T. Mendell, and S. L. Salzberg. 2015. 'StringTie enables improved reconstruction of a transcriptome from RNA-seq reads', Nat Biotechnol, 33: 290-5.

Trincado, J. L., J. C. Entizne, G. Hysenaj, B. Singh, M. Skalic, D. J. Elliott, and E. Eyras. 2018. 'SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions', Genome Biol, 19: 40.
