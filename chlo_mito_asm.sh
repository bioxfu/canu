#### Install NGMLR (https://github.com/philres/ngmlr/releases)
# cd ~
# wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
# tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz

#### Assemble chloroplast and mitochondrion
source activate gmatic
mkdir chlo_mito
cd ref 
cat Gtrilobum_chloroplast.fasta Gtrilobum_mitochondrion.fasta > Gtrilobum_chlo_mito.fasta
makeblastdb -in Gtrilobum_chloroplast.fasta -dbtype nucl
makeblastdb -in Gtrilobum_mitochondrion.fasta -dbtype nucl
cd ..

zcat Gossypium_trilobum-trim/Gossypium_trilobum.trimmedReads.fasta.gz | ~/ngmlr-0.2.7/ngmlr -t 20 -r ref/Gtrilobum_chlo_mito.fasta | samtools view -Shub -F 4 | samtools sort - -o chlo_mito/Gossypium_trilobum.trimmedReads.chlo.mito.bam

cd chlo_mito
samtools index Gossypium_trilobum.trimmedReads.chlo.mito.bam
samtools fasta Gossypium_trilobum.trimmedReads.chlo.mito.bam > Gossypium_trilobum.trimmedReads.chlo.mito.fa
samtools faidx Gossypium_trilobum.trimmedReads.chlo.mito.fa

cut -f1,2 Gossypium_trilobum.trimmedReads.chlo.mito.fa.fai|sort > Gossypium_trilobum.trimmedReads.chlo.mito.reads.length

blastn -num_threads 30 -query Gossypium_trilobum.trimmedReads.chlo.mito.fa -db ../ref/Gtrilobum_chloroplast.fasta   -outfmt 6 -evalue 1e-200 > Gossypium_trilobum.trimmedReads.chlo.blast.out
blastn -num_threads 30 -query Gossypium_trilobum.trimmedReads.chlo.mito.fa -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > Gossypium_trilobum.trimmedReads.mito.blast.out

#cat Gossypium_trilobum.trimmedReads.chlo.blast.out|cut -f1,7,8|sort -k1,1 -k2,2n -k3,3n|mergeBed|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o sum > Gossypium_trilobum.trimmedReads.chlo.blast.out.aln.length
#cat Gossypium_trilobum.trimmedReads.mito.blast.out|cut -f1,7,8|sort -k1,1 -k2,2n -k3,3n|mergeBed|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o sum > Gossypium_trilobum.trimmedReads.mito.blast.out.aln.length
cat Gossypium_trilobum.trimmedReads.chlo.blast.out|cut -f1,7,8|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o max > Gossypium_trilobum.trimmedReads.chlo.blast.out.maxaln.length
cat Gossypium_trilobum.trimmedReads.mito.blast.out|cut -f1,7,8|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o max > Gossypium_trilobum.trimmedReads.mito.blast.out.maxaln.length

cd ..
Rscript script/chlo_mito_aln_prop.R

faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_100 read_id_chlo_percent_100.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_90  read_id_chlo_percent_90.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_80  read_id_chlo_percent_80.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_70  read_id_chlo_percent_70.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_60  read_id_chlo_percent_60.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_chlo_percent_50  read_id_chlo_percent_50.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_100 read_id_mito_percent_100.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_90  read_id_mito_percent_90.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_80  read_id_mito_percent_80.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_70  read_id_mito_percent_70.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_60  read_id_mito_percent_60.fa
faSomeRecords Gossypium_trilobum.trimmedReads.chlo.mito.fa read_id_mito_percent_50  read_id_mito_percent_50.fa

~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_100 -d chlo_percent_100 -pacbio-corrected read_id_chlo_percent_100.fa genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_90  -d chlo_percent_90  -pacbio-corrected read_id_chlo_percent_90.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_80  -d chlo_percent_80  -pacbio-corrected read_id_chlo_percent_80.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_70  -d chlo_percent_70  -pacbio-corrected read_id_chlo_percent_70.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_60  -d chlo_percent_60  -pacbio-corrected read_id_chlo_percent_60.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p chlo_percent_50  -d chlo_percent_50  -pacbio-corrected read_id_chlo_percent_50.fa  genomeSize=160109 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60

~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_100 -d mito_percent_100 -pacbio-corrected read_id_mito_percent_100.fa genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_90  -d mito_percent_90  -pacbio-corrected read_id_mito_percent_90.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_80  -d mito_percent_80  -pacbio-corrected read_id_mito_percent_80.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_70  -d mito_percent_70  -pacbio-corrected read_id_mito_percent_70.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_60  -d mito_percent_60  -pacbio-corrected read_id_mito_percent_60.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p mito_percent_50  -d mito_percent_50  -pacbio-corrected read_id_mito_percent_50.fa  genomeSize=644460 correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60

scp xfu@10.41.25.100:/cluster/home/xfu/Project/Gossypium_trilobum/canu/chlo_mito/mito_percent_100/mito_percent_100.contigs.fasta .; samtools faidx mito_percent_100.contigs.fasta
scp xfu@10.41.25.100:/cluster/home/xfu/Project/Gossypium_trilobum/canu/chlo_mito/mito_percent_90/mito_percent_90.contigs.fasta .; samtools faidx mito_percent_90.contigs.fasta
scp xfu@10.41.25.100:/cluster/home/xfu/Project/Gossypium_trilobum/canu/chlo_mito/mito_percent_80/mito_percent_80.contigs.fasta .; samtools faidx mito_percent_80.contigs.fasta
scp xfu@10.41.25.100:/cluster/home/xfu/Project/Gossypium_trilobum/canu/chlo_mito/mito_percent_70/mito_percent_70.contigs.fasta .; samtools faidx mito_percent_70.contigs.fasta
scp xfu@10.41.25.100:/cluster/home/xfu/Project/Gossypium_trilobum/canu/chlo_mito/mito_percent_60/mito_percent_60.contigs.fasta .; samtools faidx mito_percent_60.contigs.fasta
scp xfu@10.41.25.100:/cluster/home/xfu/Project/Gossypium_trilobum/canu/chlo_mito/mito_percent_50/mito_percent_50.contigs.fasta .; samtools faidx mito_percent_50.contigs.fasta

cut -f1,2 mito_percent_100.contigs.fasta.fai|sort -nr -k2 > mito_percent_100.contigs_length
cut -f1,2 mito_percent_90.contigs.fasta.fai |sort -nr -k2 > mito_percent_90.contigs_length
cut -f1,2 mito_percent_80.contigs.fasta.fai |sort -nr -k2 > mito_percent_80.contigs_length
cut -f1,2 mito_percent_70.contigs.fasta.fai |sort -nr -k2 > mito_percent_70.contigs_length
cut -f1,2 mito_percent_60.contigs.fasta.fai |sort -nr -k2 > mito_percent_60.contigs_length
cut -f1,2 mito_percent_50.contigs.fasta.fai |sort -nr -k2 > mito_percent_50.contigs_length

#blastn -num_threads 30 -query mito_percent_100.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 7 -evalue 1e-200 
blastn -num_threads 30 -query mito_percent_100.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_100.contigs.blast.out
blastn -num_threads 30 -query mito_percent_90.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_90.contigs.blast.out
blastn -num_threads 30 -query mito_percent_80.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_80.contigs.blast.out
blastn -num_threads 30 -query mito_percent_70.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_70.contigs.blast.out
blastn -num_threads 30 -query mito_percent_60.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_60.contigs.blast.out
blastn -num_threads 30 -query mito_percent_50.contigs.fasta -db ../ref/Gtrilobum_mitochondrion.fasta -outfmt 6 -evalue 1e-200 > mito_percent_50.contigs.blast.out

Rscript ../script/BLAST_viewer.R mito_percent_100.contigs.blast.out mito_percent_100.contigs_length
Rscript ../script/BLAST_viewer.R mito_percent_90.contigs.blast.out mito_percent_90.contigs_length
Rscript ../script/BLAST_viewer.R mito_percent_80.contigs.blast.out mito_percent_80.contigs_length
Rscript ../script/BLAST_viewer.R mito_percent_70.contigs.blast.out mito_percent_70.contigs_length
Rscript ../script/BLAST_viewer.R mito_percent_60.contigs.blast.out mito_percent_60.contigs_length
Rscript ../script/BLAST_viewer.R mito_percent_50.contigs.blast.out mito_percent_50.contigs_length

