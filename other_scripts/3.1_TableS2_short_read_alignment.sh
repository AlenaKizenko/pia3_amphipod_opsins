## example commands to process a folder (with *1.fastq.gz and *2.fastq.gz folders)
for sample in `ls /media/quartemary/rnaseq_reads/*/*1.fastq.gz`; do base=$(basename  $sample _R1.fastq.gz); dir=$(dirname $sample); bowtie2 -x Gpu_MWS_LWS -p 9 -1 $dir/${base}_R1.fastq.gz  -2 $dir/${base}_R2.fastq.gz  --no-unal -S $base.opsins.sam; done

for sample in `ls /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_Eve*/*R1.fastq.gz`; do base=$(basename  $sample _R1.fastq.gz); dir=$(dirname $sample); bowtie2 -x Gpu_MWS_LWS -p 9 -1 $dir/${base}_R1.fastq.gz  -2 $dir/${base}_R2.fastq.gz  --no-unal -S $base.opsins.sam; done
for sample in `ls /media/tertiary/transcriptome_raw/Eulimnogammarus_cyaneus/Sample_Ecy*/*R1.fastq.gz`; do base=$(basename  $sample _R1.fastq.gz); dir=$(dirname $sample); bowtie2 -x Gpu_MWS_LWS -p 9 -1 $dir/${base}_R1.fastq.gz  -2 $dir/${base}_R2.fastq.gz  --no-unal -S $base.opsins.sam; done
for sample in `ls /media/tertiary/transcriptome_raw/Gammarus_lacustris/Sample_Gla*/*R1.fastq.gz`; do base=$(basename  $sample _R1.fastq.gz); dir=$(dirname $sample); bowtie2 -x Gpu_MWS_LWS -p 9 -1 $dir/${base}_R1.fastq.gz  -2 $dir/${base}_R2.fastq.gz  --no-unal -S $base.opsins.sam; done

## an example command to merge alignment results for the same species
for file in Sample_Eve*sam; do samtools view -bS $file >$file.bam; done
samtools merge Eve2gpu.bam *bam

## if you accidentally ran verbose and don't want to rerun everything
#for file in Sample_Eve*sam; do tail -n +5 $file | samtools view -bS  >$file.bam; done
#samtools merge Eve2eve.bam *bam



