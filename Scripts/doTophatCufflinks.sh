#!/bin/bash
SAMPLE=$1
FASTQ_DIR=rawdata
GENOME=/mnt/sequence/genomes/mouse/mm10
ANNOTATION=/mnt/sequence/genomes/mouse/mm10-iGenomes.gtf
CRG_ANNOTATION=/mnt/sequence/xzheng/Genomes/mouse_refFlat.Genes

mkdir -p $SAMPLE
fastqc -o $SAMPLE $FASTQ_DIR/$SAMPLE.fastq
tophat -o $SAMPLE/$SAMPLE\_thout -G $ANNOTATION -p4 --no-coverage-search $GENOME $FASTQ_DIR/$SAMPLE.fastq
samtools index $SAMPLE/$SAMPLE\_thout/accepted_hits.bam
samtools view $SAMPLE/$SAMPLE\_thout/accepted_hits.bam | $CRG_FLOW $CRG_ANNOTATION $SAMPLE/$SAMPLE.tagcount
cufflinks -p8 -o $SAMPLE/$SAMPLE\_clout -G $ANNOTATION $SAMPLE/$SAMPLE\_thout/accepted_hits.bam