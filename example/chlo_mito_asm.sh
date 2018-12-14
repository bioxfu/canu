#### Install NGMLR (https://github.com/philres/ngmlr/releases)
# cd ~
# wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
# tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz

#### Assemble chloroplast and mitochondrion
source activate gmatic

cd ref 
# combine references for mapping
cat Gtrilobum_chloroplast.fasta Gtrilobum_mitochondrion.fasta > Gtrilobum_chlo_mito.fasta
# build blastdb for BLAST search
makeblastdb -in Gtrilobum_chloroplast.fasta -dbtype nucl
makeblastdb -in Gtrilobum_mitochondrion.fasta -dbtype nucl
cd ..

mkdir chlo_mito
# mapping trimmed reads to chlo and mito genome
zcat Gossypium_trilobum-trim/Gossypium_trilobum.trimmedReads.fasta.gz | ~/ngmlr-0.2.7/ngmlr -t 20 -r ref/Gtrilobum_chlo_mito.fasta | samtools view -Shub -F 4 | samtools sort - -o chlo_mito/Gossypium_trilobum.trimmedReads.chlo.mito.bam

cd chlo_mito
# extract reads sequences from bam file
samtools index Gossypium_trilobum.trimmedReads.chlo.mito.bam
samtools fasta Gossypium_trilobum.trimmedReads.chlo.mito.bam > Gossypium_trilobum.trimmedReads.chlo.mito.fa

# get the reads length
samtools faidx Gossypium_trilobum.trimmedReads.chlo.mito.fa
cut -f1,2 Gossypium_trilobum.trimmedReads.chlo.mito.fa.fai|sort > Gossypium_trilobum.trimmedReads.chlo.mito.reads.length

# mapping the aligned reads to chlo and mito genome again using BLAST
blastn -num_threads 30 -query Gossypium_trilobum.trimmedReads.chlo.mito.fa -db ../ref/Gtrilobum_chloroplast.fasta   -outfmt 6 -evalue 1e-200 > Gossypium_trilobum.trimmedReads.chlo.blast.out
blastn -num_threads 30 -query Gossypium_trilobum.trimmedReads.chlo.mito.fa -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > Gossypium_trilobum.trimmedReads.mito.blast.out

# find the length of maximum aligned region for each read
cat Gossypium_trilobum.trimmedReads.chlo.blast.out|cut -f1,7,8|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o max > Gossypium_trilobum.trimmedReads.chlo.blast.out.maxaln.length
cat Gossypium_trilobum.trimmedReads.mito.blast.out|cut -f1,7,8|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o max > Gossypium_trilobum.trimmedReads.mito.blast.out.maxaln.length

cd ..
# classify the aligned reads into different categories based on the proportion of reads that can be aligned on the reference genome
Rscript script/chlo_mito_aln_prop.R

# extract the sequences based on the read iDs using faSomeRecords
# https://bioconda.github.io/recipes/ucsc-fasomerecords/README.html
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_100 read_id_chlo_percent_100.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_90  read_id_chlo_percent_90.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_80  read_id_chlo_percent_80.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_70  read_id_chlo_percent_70.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_60  read_id_chlo_percent_60.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_50  read_id_chlo_percent_50.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_100 read_id_mito_percent_100.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_90  read_id_mito_percent_90.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_80  read_id_mito_percent_80.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_70  read_id_mito_percent_70.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_60  read_id_mito_percent_60.fa
#faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_50  read_id_mito_percent_50.fa

# assemble the mito and chlo genomes
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_100 -d chlo_percent_100 -pacbio-corrected read_id_chlo_percent_100.fa genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_90  -d chlo_percent_90  -pacbio-corrected read_id_chlo_percent_90.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_80  -d chlo_percent_80  -pacbio-corrected read_id_chlo_percent_80.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_70  -d chlo_percent_70  -pacbio-corrected read_id_chlo_percent_70.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_60  -d chlo_percent_60  -pacbio-corrected read_id_chlo_percent_60.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_50  -d chlo_percent_50  -pacbio-corrected read_id_chlo_percent_50.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_100 -d mito_percent_100 -pacbio-corrected read_id_mito_percent_100.fa genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_90  -d mito_percent_90  -pacbio-corrected read_id_mito_percent_90.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_80  -d mito_percent_80  -pacbio-corrected read_id_mito_percent_80.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_70  -d mito_percent_70  -pacbio-corrected read_id_mito_percent_70.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_60  -d mito_percent_60  -pacbio-corrected read_id_mito_percent_60.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
#~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_50  -d mito_percent_50  -pacbio-corrected read_id_mito_percent_50.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60

# calculate the length of contigs
cp mito_percent_100/mito_percent_100.contigs.fasta .; samtools faidx mito_percent_100.contigs.fasta
#cp mito_percent_90/mito_percent_90.contigs.fasta .; samtools faidx mito_percent_90.contigs.fasta
#cp mito_percent_80/mito_percent_80.contigs.fasta .; samtools faidx mito_percent_80.contigs.fasta
#cp mito_percent_70/mito_percent_70.contigs.fasta .; samtools faidx mito_percent_70.contigs.fasta
#cp mito_percent_60/mito_percent_60.contigs.fasta .; samtools faidx mito_percent_60.contigs.fasta
#cp mito_percent_50/mito_percent_50.contigs.fasta .; samtools faidx mito_percent_50.contigs.fasta
cp chlo_percent_100/chlo_percent_100.contigs.fasta .; samtools faidx chlo_percent_100.contigs.fasta
#cp chlo_percent_90/chlo_percent_90.contigs.fasta .; samtools faidx chlo_percent_90.contigs.fasta
#cp chlo_percent_80/chlo_percent_80.contigs.fasta .; samtools faidx chlo_percent_80.contigs.fasta
#cp chlo_percent_70/chlo_percent_70.contigs.fasta .; samtools faidx chlo_percent_70.contigs.fasta
#cp chlo_percent_60/chlo_percent_60.contigs.fasta .; samtools faidx chlo_percent_60.contigs.fasta
#cp chlo_percent_50/chlo_percent_50.contigs.fasta .; samtools faidx chlo_percent_50.contigs.fasta
cut -f1,2 mito_percent_100.contigs.fasta.fai|sort -nr -k2 > mito_percent_100.contigs_length
#cut -f1,2 mito_percent_90.contigs.fasta.fai |sort -nr -k2 > mito_percent_90.contigs_length
#cut -f1,2 mito_percent_80.contigs.fasta.fai |sort -nr -k2 > mito_percent_80.contigs_length
#cut -f1,2 mito_percent_70.contigs.fasta.fai |sort -nr -k2 > mito_percent_70.contigs_length
#cut -f1,2 mito_percent_60.contigs.fasta.fai |sort -nr -k2 > mito_percent_60.contigs_length
#cut -f1,2 mito_percent_50.contigs.fasta.fai |sort -nr -k2 > mito_percent_50.contigs_length
cut -f1,2 chlo_percent_100.contigs.fasta.fai|sort -nr -k2 > chlo_percent_100.contigs_length 
#cut -f1,2 chlo_percent_90.contigs.fasta.fai |sort -nr -k2 > chlo_percent_90.contigs_length 
#cut -f1,2 chlo_percent_80.contigs.fasta.fai |sort -nr -k2 > chlo_percent_80.contigs_length 
#cut -f1,2 chlo_percent_70.contigs.fasta.fai |sort -nr -k2 > chlo_percent_70.contigs_length 
#cut -f1,2 chlo_percent_60.contigs.fasta.fai |sort -nr -k2 > chlo_percent_60.contigs_length 
#cut -f1,2 chlo_percent_50.contigs.fasta.fai |sort -nr -k2 > chlo_percent_50.contigs_length 

# BLAST the contigs to reference genome
blastn -num_threads 30 -query mito_percent_100.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_100.contigs.blast.out
#blastn -num_threads 30 -query mito_percent_90.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_90.contigs.blast.out
#blastn -num_threads 30 -query mito_percent_80.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_80.contigs.blast.out
#blastn -num_threads 30 -query mito_percent_70.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_70.contigs.blast.out
#blastn -num_threads 30 -query mito_percent_60.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_60.contigs.blast.out
#blastn -num_threads 30 -query mito_percent_50.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_50.contigs.blast.out
blastn -num_threads 30 -query chlo_percent_100.contigs.fasta -db ../ref/Gtrilobum_chloroplast.fasta -outfmt 6 -evalue 1e-200 > chlo_percent_100.contigs.blast.out
#blastn -num_threads 30 -query chlo_percent_90.contigs.fasta -db ../ref/Gtrilobum_chloroplast.fasta -outfmt 6 -evalue 1e-200 > chlo_percent_90.contigs.blast.out
#blastn -num_threads 30 -query chlo_percent_80.contigs.fasta -db ../ref/Gtrilobum_chloroplast.fasta -outfmt 6 -evalue 1e-200 > chlo_percent_80.contigs.blast.out
#blastn -num_threads 30 -query chlo_percent_70.contigs.fasta -db ../ref/Gtrilobum_chloroplast.fasta -outfmt 6 -evalue 1e-200 > chlo_percent_70.contigs.blast.out
#blastn -num_threads 30 -query chlo_percent_60.contigs.fasta -db ../ref/Gtrilobum_chloroplast.fasta -outfmt 6 -evalue 1e-200 > chlo_percent_60.contigs.blast.out
#blastn -num_threads 30 -query chlo_percent_50.contigs.fasta -db ../ref/Gtrilobum_chloroplast.fasta -outfmt 6 -evalue 1e-200 > chlo_percent_50.contigs.blast.out

# visualize the BLAST results
Rscript ../script/BLAST_viewer.R mito_percent_100.contigs.blast.out mito_percent_100.contigs_length 10000
#Rscript ../script/BLAST_viewer.R mito_percent_90.contigs.blast.out mito_percent_90.contigs_length 10000
#Rscript ../script/BLAST_viewer.R mito_percent_80.contigs.blast.out mito_percent_80.contigs_length 10000
#Rscript ../script/BLAST_viewer.R mito_percent_70.contigs.blast.out mito_percent_70.contigs_length 10000
#Rscript ../script/BLAST_viewer.R mito_percent_60.contigs.blast.out mito_percent_60.contigs_length 10000
#Rscript ../script/BLAST_viewer.R mito_percent_50.contigs.blast.out mito_percent_50.contigs_length 10000
Rscript ../script/BLAST_viewer.R chlo_percent_100.contigs.blast.out chlo_percent_100.contigs_length 10000
#Rscript ../script/BLAST_viewer.R chlo_percent_90.contigs.blast.out chlo_percent_90.contigs_length 10000
#Rscript ../script/BLAST_viewer.R chlo_percent_80.contigs.blast.out chlo_percent_80.contigs_length 10000
#Rscript ../script/BLAST_viewer.R chlo_percent_70.contigs.blast.out chlo_percent_70.contigs_length 10000
#Rscript ../script/BLAST_viewer.R chlo_percent_60.contigs.blast.out chlo_percent_60.contigs_length 10000
#Rscript ../script/BLAST_viewer.R chlo_percent_50.contigs.blast.out chlo_percent_50.contigs_length 10000

# extract reads for nuclear genome assembly
cat read_id_chlo_percent_100 read_id_mito_percent_100 > read_id_chlo_and_mito_percent_100
faSomeRecords -exclude ../Gossypium_trilobum-trim/Gossypium_trilobum.trimmedReads.fasta.gz read_id_chlo_and_mito_percent_100 Gossypium_trilobum_nucl.trimmedReads.fasta
gzip Gossypium_trilobum_nucl.trimmedReads.fasta
cd ..

# assemble the nuclear genome
NAME=Gossypium_trilobum_nucl
GENOME_SIZE=851m
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p ${NAME} -d ${NAME}-erate-0.045 -pacbio-corrected chlo_mito/${NAME}.trimmedReads.fasta.gz genomeSize=${GENOME_SIZE} correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60

