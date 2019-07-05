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
## login head4 IP (10.41.25.99)
ssh head4

## you can run it in three steps:
## step1: correct
~/canu-1.7.1/Linux-amd64/bin/canu -correct -p ${NAME} -d ${NAME}-trim-corOutCoverage40 -pacbio-raw raw/${NAME}.subreads.fastq.gz genomeSize=${GENOME_SIZE} gnuplotImageFormat=svg gridOptions="-t 1000:00:00 --mem-per-cpu=8g" gridEngine=slurm canuIterationMax=100 corOutCoverage=40

## step2: trim
~/canu-1.7.1/Linux-amd64/bin/canu -trim -p ${NAME} -d ${NAME}-trim-corOutCoverage40 -pacbio-corrected ${NAME}-trim-corOutCoverage40/${NAME}.correctedReads.fasta.gz genomeSize=${GENOME_SIZE} gnuplotImageFormat=svg gridOptions="-t 1000:00:00 --mem-per-cpu=8g" gridEngine=slurm canuIterationMax=100 corOutCoverage=40

## step3: assemble
~/canu-1.7.1/Linux-amd64/bin/canu -assemble -p ${NAME} -d ${NAME}-erate-0.045-corOutCoverage40 -pacbio-corrected ${NAME}-trim-corOutCoverage40/${NAME}.trimmedReads.fasta.gz genomeSize=${GENOME_SIZE} correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00 --mem-per-cpu=8g" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60 corOutCoverage=40 

## Or in one step
~/canu-1.7.1/Linux-amd64/bin/canu -p ${NAME} -d ${NAME}-erate-0.045-corOutCoverage40 -pacbio-raw raw/${NAME}.subreads.fasta.gz genomeSize=${GENOME_SIZE} correctedErrorRate=0.045 gnuplotImageFormat=svg gridOptions="-t 1000:00:00 --mem-per-cpu=8g" gridEngine=slurm canuIterationMax=100 redMemory=60 oeaMemory=60 corOutCoverage=40

## check the status of your job
squeue -l

## remove all of your job
scancel -u xfu
```

#### Evaluate the assembly
```
module add quast/4.6.1
quast.py ${NAME}.contigs.fasta 
```

#### Tip 1
When you see the following error in *canu.out*:
```
job correction/${NAME}.ovlStore.BUILDING/1008 FAILED
```
and you check the logs: *correction/${NAME}.ovlStore.BUILDING/logs/2-sort.59339_1008.out*

you would see the info like this:
```
Running job 1008 based on SLURM_ARRAY_TASK_ID=1008 and offset=0.

Attempting to increase maximum allowed processes and open files.
  Max processes per user limited to 514848, no increase possible.
  Changed max open files from 51200 to 65536 (max 65536).

Job 1008 is finished (remove './1008' or -force to try again).
```
How to solve this problem:

1. edit the file: *correction/${NAME}.ovlStore.BUILDING/scripts/2-sort.sh*

2. add **-force** to ovStoreSorter
```
$bin/ovStoreSorter -force
```
3. run Canu again.

if you still can't fix the problem, just remove the *${NAME}.ovlStore.BUILDING* and run it again.

#### Tip 2
If Bogart failed again and again, run it on node2 with 1024Gb memory (-M 1024)
```
ssh node2 
cd ${NAME}-erate-0.045/unitigging/4-unitigger

~/canu-1.7.1/Linux-amd64/bin/bogart -G ../${NAME}.gkpStore -O ../${NAME}.ovlStore -o ./${NAME} -gs 1400000000 -eg 0.045 -eM 0.045 -mo 500 -dg 6 -db 6 -dr 3 -ca 2100 -cp 200 -threads 16 -M 1024 -unassembled 2 0 1.0 0.5 3 > ./unitigger.err 2>&1

mv ./qq74.ctgStore ../qq74.ctgStore
mv ./qq74.utgStore ../qq74.utgStore
~/canu-1.7.1/Linux-amd64/bin/tgStoreDump -G ../qq74.gkpStore -T ../qq74.ctgStore 1 -sizes -s 1400000000 > ../qq74.ctgStore/seqDB.v001.sizes.txt

# run Canu again 
```

#### Tip 3
If meryl failed again and again, run it on node2
```
ssh node2 
cd Cq_PBL191261-trim-corOutCoverage40/trimming/0-mercounts

/cluster/home/xfu/canu-1.7.1/Linux-amd64/bin/meryl -B -C -L 2 -v -m 22 -threads 32 -memory 209715 -s ../Cq_PBL191261.gkpStore -o ./Cq_PBL191261.ms22.WORKING 

mv ./Cq_PBL191261.ms22.WORKING.mcdat ./Cq_PBL191261.ms22.mcdat
mv ./Cq_PBL191261.ms22.WORKING.mcidx ./Cq_PBL191261.ms22.mcidx

/cluster/home/xfu/canu-1.7.1/Linux-amd64/bin/meryl -Dh -s ./Cq_PBL191261.ms22 > ./Cq_PBL191261.ms22.histogram.WORKING 2> ./Cq_PBL191261.ms22.histogram.info 
mv -f ./Cq_PBL191261.ms22.histogram.WORKING ./Cq_PBL191261.ms22.histogram

/cluster/home/xfu/canu-1.7.1/Linux-amd64/bin/estimate-mer-threshold -h ./Cq_PBL191261.ms22.histogram -c 38 > ./Cq_PBL191261.ms22.estMerThresh.out.WORKING 2> ./Cq_PBL191261.ms22.estMerThresh.err
mv ./Cq_PBL191261.ms22.estMerThresh.out.WORKING ./Cq_PBL191261.ms22.estMerThresh.out
```
