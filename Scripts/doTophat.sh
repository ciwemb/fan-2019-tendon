#!/bin/bash
SAMPLE=$1
FASTQ_DIR=rawdata
GENOME=/mnt/sequence/genomes/mouse/mm9
ANNOTATION=/mnt/sequence/xzheng/Genomes/mm9_refFlat.gtf
CRG_FLOW=/mnt/sequence/xzheng/bin/CRG-Flow
CRG_ANNOTATION=/mnt/sequence/xzheng/Genomes/mouse_refFlat.Genes

mkdir -p $SAMPLE
fastqc -o $SAMPLE $FASTQ_DIR/$SAMPLE.fastq
tophat -o $SAMPLE/$SAMPLE\_thout -G $ANNOTATION -p4 --no-coverage-search $GENOME $FASTQ_DIR/$SAMPLE.fastq
samtools index $SAMPLE/$SAMPLE\_thout/accepted_hits.bam
samtools view $SAMPLE/$SAMPLE\_thout/accepted_hits.bam | $CRG_FLOW $CRG_ANNOTATION $SAMPLE/$SAMPLE.tagcount


