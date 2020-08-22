# PP2PG
Phasing pipeline based on two parental genomes

Require: 
MUMmer     Version:4.0.0beta2 (https://github.com/mummer4/mummer)
Bismark    Version:v0.22.1    (https://github.com/FelixKrueger/Bismark)
minimap2   Version:2.17       (https://github.com/lh3/minimap2)
HISAT2     Version:2.1.0      (http://ccb.jhu.edu/software/hisat2/manual.shtml)
pysam      Version:0.15.2     (https://pypi.org/project/pysam/)     
SAMtools   Version:1.6        (http://samtools.sourceforge.net/)  


1.Mapping:

1) The type of Iso-seq reads was used minimap2 software for alignment:

       minimap2 -ax splice -uf --secondary=no -C5 -O6,24 -B4 --MD ref query > results.sam      

2) The type of RNA-seq paired end reads was utilized hisat2 software for alignment:

       hisat2 --dta -k 1 -x genome_index -1 query_1 -2 query_2 -S out.sam

3) The type of WGBS reads was used Bismark software for alignment:
      
       bismark --bowtie2 -genome genome_index -1 query_1 -2 query_2
    Note: Our pipeline only considers all unique mapped results in the current genomes in BS-Seq. Ambiguous reads are regarded 
    as unmaped reads. 



2.SNP Calling:

    nucmer --mum --maxgap=500 --mincluster=100 --prefix=g1_g2 g1_genome.fa g1_genome.fa
    
    delta-filter -1 -q -r g1_g2.delta > g1_g2.filter.delta
    
    show-snps -C -H -I -T -r -l g1_g2.filter.delta > g1_g2.snp


3.Scoring of Phasing and Classification Output:

    phasing.py --pt=<phase_type> --g1=<genome1_alignment> --g2=<genome2_alignment> 
               --snp=<snpfile> --gop1=<genome1_output_prefix> --gop2=<genome2_output_prefix>


