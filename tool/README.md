phasing.py --pt=<phase_type> --g1=<genome1_alignment> --g2=<genome2_alignment> 
           --snp=<snpfile> --gop1=<genome1_output_prefix> --gop2=<genome2_output_prefix> [-h] 

    This is a phasing tool based on two genomes. 

    --pt=<phase_type> 
    Type: Isoseq : Iso-Seq (PacBio Isoform Sequence) 
          RNAseq : RNA-Seq for paired-end reads 
          BSseq  : BS-Seq (Bisulfite Sequencing) for paired-end reads

    --g1=<genome1_alignment>
    The bam/sam file that the reads aligned to first genome.

    The detailed methods were as follows:
    1) The type of Iso-Seq reads was used minimap2 software for alignment.
     eg. minimap2 -ax splice -uf --secondary=no -C5 -O6,24 -B4 --MD ref query > results.sam      

    2) The type of RNA-Seq paired end reads was utilized hisat2 software for alignment.
     eg. hisat2 --dta -k 1 -x genome_index -1 query_1 -2 query_2 -S out.sam

    3) The type of BS-Seq reads was used Bismark software for alignment.
     eg. bismark --bowtie2 -genome genome_index -1 query_1 -2 query_2
    Note: Our pipeline only considers all unique mapped results in the current genomes in BS-Seq. Ambiguous reads are regarded 
    as unmaped reads.

    --g2=<genome2_alignment> 
    The bam/sam file that the reads aligned to second genome.

    --snp=<snpfile> 
    The SNP file was obtained by mummer program which used first genome as reference genome and second genome as query genome.
     eg. nucmer --mum --maxgap=500 --mincluster=100 --prefix=g1_g2 g1_genome.fa g1_genome.fa
         delta-filter -1 -q -r g1_g2.delta > g1_g2.filter.delta
         show-snps -C -H -I -T -r -l g1_g2.filter.delta > g1_g2.snp

    --st=<INT>
    The default is 1. Number of threads using samtools sort

    --gop1=<genome1_output_prefix> 
    The prefix file for first genome.

    --gop2=<genome2_output_prefix>
    The prefix file for second genome.

    -h help

Require: 

    pysam     Version: 0.15.2
    samtools  Version: 1.6
    
Output file:

    Two genome-based alignment results file, containing 2 folders. One is the phasing results of g1 genome alignment, 
    which includes: g1-synteny reads, g2-synteny reads, unknown reads, and reads that only mapped to the g1 genome (g1-only reads); 
    the other is the phasing results of g2 genome alignment, which includes: g1-synteny reads, g2-synteny reads, 
    unknown reads, and reads that only mapped to the g2 genome (g2-only reads).
    g1_genome/g1reads.bam --- g1-synteny reads mapped to g1 genome
    g1_genome/g2reads.bam --- g2-synteny reads mapped to g1 genome
    g1_genome/unknowreads.bam --- unknown reads mapped to g1 genome
    g1_genome/g1onlyreads.bam --- reads that only mapped to the g1 genome (g1-only reads)
    g2_genome/g1reads.bam --- g1-synteny reads mapped to g2 genome
    g2_genome/g2reads.bam --- g2-synteny reads mapped to g2 genome
    g2_genome/unknowreads.bam --- unknown reads mapped to g2 genome
    g2_genome/g2onlyreads.bam --- reads that only mapped to the g2 genome (g2-only reads)

    g1_g2_Phasing_Report.txt --- Final Phasing Reads Report
