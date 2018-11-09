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

cat Gossypium_trilobum.trimmedReads.chlo.blast.out|cut -f1,7,8|sort -k1,1 -k2,2n -k3,3n|mergeBed|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o sum > Gossypium_trilobum.trimmedReads.chlo.blast.out.aln.length
cat Gossypium_trilobum.trimmedReads.mito.blast.out|cut -f1,7,8|sort -k1,1 -k2,2n -k3,3n|mergeBed|awk '{print $1"\t"$3-$2+1}'|sort -k1,1 -k2,2n|groupBy -g 1 -c 2 -o sum > Gossypium_trilobum.trimmedReads.mito.blast.out.aln.length

cd ..

script/chlo_mito_aln_prop.R

