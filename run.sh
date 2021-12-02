#!/usr/bin/env bash
## Set these variables to your liking.
READSFILE=/media/DataXFS/ngmlr_data/PacBio_CCS_10kb/3000_reads.fastq
REFFILE=/media/DataXFS/ngmlr_data/GRCh38/genome_assemblies_genome_fasta/ncbi-genomes-2021-11-16/GCA_000001405.15_GRCh38_genomic.fna.gz
THREADS=1
OPTIMIZE=pacbio

#!/bin/bash
for SUBSEGMENT in $@
do
    # These are automatically generated variables
    # SUBSEGMENT=$1
    DIR=.
    NICKNAME="$OPTIMIZE-$(basename $READSFILE | sed s/.fastq$//)"
    DATE=$(date --iso-8601=seconds)
    GPROFFILENAME=gprof.$NICKNAME.$SUBSEGMENT.txt
    OUTPUTFILENAME=$NICKNAME.result.$SUBSEGMENT.sam
    OUTPUTDIR=$DIR/docs/$NICKNAME/threads-$THREADS/

    mkdir -p $OUTPUTDIR


    start=`date +%s`
    $DIR/bin/ngmlr-0.2.8-debug/ngmlr-debug --bam-fix -x $OPTIMIZE -t $THREADS --subread-length $SUBSEGMENT -q $READSFILE -r $REFFILE -o $OUTPUTDIR/$OUTPUTFILENAME
    end=`date +%s`
    runtime=$((end-start))
    echo "+++++ Clock time: $runtime"
    echo "+++++ Clock time: $runtime" > $OUTPUTDIR/$GPROFFILENAME

    gprof $DIR/bin/ngmlr-0.2.8-debug/ngmlr-debug >> $OUTPUTDIR/$GPROFFILENAME    


    # convert to BAM
    samtools sort -O BAM -@11 -o $OUTPUTDIR/$OUTPUTFILENAME.bam $OUTPUTDIR/$OUTPUTFILENAME
    # create index for visualization in IGV
    samtools index $OUTPUTDIR/$OUTPUTFILENAME.bam
done


# diff -y --suppress-common-lines PacBio_3000_head.result.256-2021-11-26T16:15:58+09:00.sam.sort.txt PacBio_3000_head.result.1024-2021-11-26T16:30:40+09:00.sam.sort.txt | wc -l