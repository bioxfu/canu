### Canu Assembly Workflow
#### Install Canu
```
## https://github.com/marbl/canu/releases
cd ~
wget https://github.com/marbl/canu/releases/download/v1.7.1/canu-1.7.1.Linux-amd64.tar.xz
xz -dc canu-1.7.1.*.tar.xz |tar -xf -
```

#### Initiate the project
```
cp example/init.sh .
# edit and run init.sh file according to your project
. init.sh
```

#### Prepare raw data
```
## Convert BAM files to fasta file
## The subreads.bam contains all the usable data. 
## https://pacbiofileformats.readthedocs.io/en/3.0/Primer.html
## http://seqanswers.com/forums/showpost.php?p=202171&postcount=9

nohup bash -c "find $DATA_DIR/*/*/*subreads.bam | xargs -I {} samtools fasta {} | gzip -c > raw/$NAME.subreads.fasta.gz" &
```

#### Run Canu on HPC
```
## head4 IP (10.41.25.99)
## correct
~/canu-1.7.1/Linux-amd64/bin/canu -correct -p ${NAME} -d ${NAME}-trim -pacbio-raw raw/${NAME}.subreads.fasta.gz genomeSize=${GENOME_SIZE} gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100
## trim
~/canu-1.7.1/Linux-amd64/bin/canu -trim -p ${NAME} -d ${NAME}-trim -pacbio-corrected ${NAME}-trim/${NAME}.correctedReads.fasta.gz genomeSize=${GENOME_SIZE} gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100
## assemble
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p ${NAME} -d ${NAME}-erate-0.045 -pacbio-corrected ${NAME}-trim/${NAME}.trimmedReads.fasta.gz genomeSize=${GENOME_SIZE} correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60

## check the status of your job
squeue -l
```

#### Evaluate the assembly
```
module add quast/4.6.1
quast.py ${NAME}.contigs.fasta 
```
