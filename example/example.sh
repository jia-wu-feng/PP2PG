###Unzip files
gzip -d MH63.fa.gz ZS97.fa.gz Isoseq.fa.gz BSseq_1.fq.gz BSseq_2.fq.gz RNAseq_1.fq.gz RNAseq_2.fq.gz

###Iso-Seq phasing
mkdir Isoseq
cd Isoseq

######Calling allelic SNPs by MUMmer
nucmer --mum --maxgap=500 --mincluster=100 --prefix=MH63_ZS97 ../MH63.fa ../ZS97.fa
delta-filter -1 -q -r MH63_ZS97.delta > MH63_ZS97.filter.delta
show-snps -C -H -I -T -r -l MH63_ZS97.filter.delta > MH63_ZS97.snp

######Mapping to the two parental genomes for alignment by minimap2
minimap2 -ax splice -uf --secondary=no -C5 -O6,24 -B4 --MD ../MH63.fa ../Isoseq.fa > MH63.sam
minimap2 -ax splice -uf --secondary=no -C5 -O6,24 -B4 --MD ../ZS97.fa ../Isoseq.fa > ZS97.sam

######Scoring and Classification
../../tool/phasing.py --pt=Isoseq --g1=MH63.sam --g2=ZS97.sam --snp=MH63_ZS97.snp --gop1=MH63 --gop2=ZS97

cd ../
###BS-Seq phasing
mkdir BSseq
cd BSseq

######Calling allelic SNPs by MUMmer
nucmer --mum --maxgap=500 --mincluster=100 --prefix=MH63_ZS97 ../MH63.fa ../ZS97.fa
delta-filter -1 -q -r MH63_ZS97.delta > MH63_ZS97.filter.delta
show-snps -C -H -I -T -r -l MH63_ZS97.filter.delta > MH63_ZS97.snp

#########Build index by bismark
mkdir genome
cd genome
mkdir MH63 ZS97
cp ../../MH63.fa MH63/. 
cp ../../ZS97.fa ZS97/.
bismark_genome_preparation --bowtie2 MH63
bismark_genome_preparation --bowtie2 ZS97

######Mapping to the two parental genomes for alignment by bismark
cd ../
bismark --bowtie2 -B MH63 -genome genome/MH63 -1 ../BSseq_1.fq -2 ../BSseq_2.fq
bismark --bowtie2 -B ZS97 -genome genome/ZS97 -1 ../BSseq_1.fq -2 ../BSseq_2.fq

######Scoring and Classification
../../tool/phasing.py --pt=BSseq --g1=MH63_pe.bam --g2=ZS97_pe.bam --snp=MH63_ZS97.snp --gop1=MH63 --gop2=ZS97

cd ../
###RNA-seq phasing
mkdir RNAseq
cd RNAseq

######Calling allelic SNPs by MUMmer
nucmer --mum --maxgap=500 --mincluster=100 --prefix=MH63_ZS97 ../MH63.fa ../ZS97.fa
delta-filter -1 -q -r MH63_ZS97.delta > MH63_ZS97.filter.delta
show-snps -C -H -I -T -r -l MH63_ZS97.filter.delta > MH63_ZS97.snp

######Build index by hisat2
mkdir genome
cd genome
mkdir MH63 ZS97
cp ../../MH63.fa MH63/. 
cp ../../ZS97.fa ZS97/.
cd MH63/
hisat2-build MH63.fa MH63
cd ../ZS97
hisat2-build ZS97.fa ZS97

######Mapping to the two parental genomes for alignment by hisat2
cd ../../
hisat2 --dta -k 1 -x genome/MH63/MH63 -1 ../RNAseq_1.fq -2 ../RNAseq_2.fq -S RNA.MH63.sam
hisat2 --dta -k 1 -x genome/ZS97/ZS97 -1 ../RNAseq_1.fq -2 ../RNAseq_2.fq -S RNA.ZS97.sam

######Scoring and Classification
../../tool/phasing.py --pt=RNAseq --g1=RNA.MH63.sam --g2=RNA.ZS97.sam --snp=MH63_ZS97.snp --gop1=MH63 --gop2=ZS97
