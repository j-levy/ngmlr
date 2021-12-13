#!/usr/bin/env bash
## Set these variables to your liking.
READSFILE=/media/DataXFS/ngmlr_data/Genomes/c_elegans/pbsim_coverage10/sd_0004.fastq
REFFILE=/media/DataXFS/ngmlr_data/Genomes/c_elegans/refs/sd_0004_sv0.fasta
# READSFILE=/media/DataXFS/ngmlr_data/Genomes/Amanita_phalloides/pbsim_cov1/amanita_reads.fastq
# REFFILE=/media/DataXFS/ngmlr_data/Genomes/Amanita_phalloides/survivor/amanita0.0.fasta # survivor_amanita/amanita0.1smol.fasta
THREADS=12
# this means that at least SNIFFLES_MIN_COVERAGE reads must agree on an SV to be detected.
# you should use a dataset using a coverage large enough to be sure you detect your SVs.
# I don't know exactly what "large enough" is. In production they often go for a coverage of 40-50.
# this means that statistically DNA is covered around 50 times, this redundancy allows to be certain about mapping.
SNIFFLES_MIN_COVERAGE=3
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
    $DIR/bin/ngmlr-0.2.8/ngmlr \
        --bam-fix -x $OPTIMIZE -t $THREADS \
        --subread-length $SUBSEGMENT \
        -q $READSFILE -r $REFFILE \
        -o $OUTPUTDIR/$OUTPUTFILENAME #--verbose 2> $SUBSEGMENT.log
    end=`date +%s`
    runtime=$((end-start))
    echo "+++++ Clock time: $runtime"
    echo "+++++ Clock time: $runtime" > $OUTPUTDIR/$GPROFFILENAME

    gprof $DIR/bin/ngmlr-0.2.8/ngmlr >> $OUTPUTDIR/$GPROFFILENAME    

    # convert to BAM
    # create index for visualization in IGV
    # using samtools = 1.9
    $HOME/miniconda3/envs/sniffles/bin/samtools view -o $OUTPUTDIR/$OUTPUTFILENAME.bam $OUTPUTDIR/$OUTPUTFILENAME
    $HOME/miniconda3/envs/sniffles/bin/samtools sort -o $OUTPUTDIR/$OUTPUTFILENAME.sort.bam $OUTPUTDIR/$OUTPUTFILENAME.bam 
    $HOME/miniconda3/envs/sniffles/bin/samtools index $OUTPUTDIR/$OUTPUTFILENAME.sort.bam
    $HOME/miniconda3/envs/sniffles/bin/sniffles -s $SNIFFLES_MIN_COVERAGE -m $OUTPUTDIR/$OUTPUTFILENAME.sort.bam -v $OUTPUTDIR/$OUTPUTFILENAME.vcf
    # cat $OUTPUTDIR/$OUTPUTFILENAME.vcf
done

# diff -y --suppress-common-lines PacBio_3000_head.result.256-2021-11-26T16:15:58+09:00.sam.sort.txt PacBio_3000_head.result.1024-2021-11-26T16:30:40+09:00.sam.sort.txt | wc -l