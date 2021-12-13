quick notes on how to evaluate ngmlr

take data set (fasta)

## simulate pacbio reads

> main parameters : 
> - depth (min: 1, prod: 40-50, O(n))
> - length (min: ?, prod: ?, ?)

use pbsim to split the fasta into chromosomes and chromosomes reads.
each sd_00xx is a chromosome fasta reference (.ref) and simulated read (.fastq) with a given error model.
I always use P6C4 for no particular reason (it was suggested in pbsim2 tutorial)
running it :

pbsim --depth NUMBER --hmm_model PATH/TO/ERROR/MODEL PATH/TO/FASTA/FILE

pbsim --depth 15 --hmm_model ../../pbsim2/data/P6C4~stuff ../genome.fasta 

You can also adjust other stuff, including
  --length-mean        mean of length model (9000.0).
  --length-sd          standard deviation of length model (7000.0).
  --accuracy-mean      mean of accuracy model (0.85).


Note : depth (or coverage) indicates how many time DNA must be covered by sampling.
In production they often go for 50, meaning that if your chromosome reference is 40MB, your reads is 50x40MB = 2GB 
Yes, this produces *a huge amount* of data.

You now have sd_00xx.fastq (simulated reads for chromosome xx) and sd_00xx.ref (fasta reference for chromosome xx)

## simulating SV

> main parameters:
> - number of SV, type of SV, size of SV, frequency of SNP

create a parameter file for SURVIVOR

SURVIVOR simSV param.svr 

edit param.svr to your liking
If you're doing only 1 chromosome mapping, disable translocation since you need at least 2 chromosomes (ofc)

generate the altered reference fasta with:

SURVIVOR simSV PATH/TO/REFERENCE.fasta PATH/TO/param.svr FREQUENCY_OF_SNP(0..1) 0(always using simulated reads) OUTPUT_PREFIX

SURVIVOR simSV sd_00xx.ref param.svr 0.1 0 sd_00xx_sv.fasta

(note : all .ref are just like .fasta, simply they're split by chromosome)

this generates:

- sd_00xx_sv.fasta: reference with SVs
- sd_00xx_sv.vcf: list of SV added
- sd_00xx_sv.bed: reference of SV added (slightly different from VCF I guess?), will be used for evaluation.

## Run ngmlr

> main parameters:
> - kmer length (default: 13)
> - subsegment length (default: 256)


run the thing with ./run.sh after having adjusted variables inside, it's the simplest thing.

It generates an alignment file (SAM).

## Convert SAM to BAM

use samtools for conversion. These are usually reasonably fast.

    $HOME/miniconda3/envs/sniffles/bin/samtools view -o $OUTPUTDIR/$OUTPUTFILENAME.bam $OUTPUTDIR/$OUTPUTFILENAME
    $HOME/miniconda3/envs/sniffles/bin/samtools sort -o $OUTPUTDIR/$OUTPUTFILENAME.sort.bam $OUTPUTDIR/$OUTPUTFILENAME.bam 
    $HOME/miniconda3/envs/sniffles/bin/samtools index $OUTPUTDIR/$OUTPUTFILENAME.sort.bam

the index is needed to visualize it in the interactive genome viewer, `igv`
    
    
## Call for SV with sniffles

> main parameters
> - minimum coverage (def: 10), **must be inferior to depth if you want to see anything, of course**
> - lots of other stuff I haven't seen yet
>>> General:
    -s <int>,  --min_support <int>
        Minimum number of reads that support a SV. [10]
    --max_num_splits <int>
        Maximum number of splits per read to be still taken into account. [7]
    -d <int>,  --max_distance <int>
        Maximum distance to group SV together. [1000]
    -t <int>,  --threads <int>
        Number of threads to use. [3]
    -l <int>,  --min_length <int>
        Minimum length of SV to be reported. [30]
    -q <int>,  --minmapping_qual <int>
        Minimum Mapping Quality. [20]
    -n <int>,  --num_reads_report <int>
        Report up to N reads that support the SV in the vcf file. -1: report all. [0]
    -r <int>,  --min_seq_size <int>
        Discard read if non of its segment is larger then this. [2000]
    -z <int>,  --min_zmw <int>
        Discard SV that are not supported by at least x zmws. This applies only for PacBio recognizable reads. [0]
    --cs_string
        Enables the scan of CS string instead of Cigar and MD.  [false]



    $HOME/miniconda3/envs/sniffles/bin/sniffles -s $SNIFFLES_MIN_COVERAGE -m $OUTPUTDIR/$OUTPUTFILENAME.sort.bam -v $OUTPUTDIR/$OUTPUTFILENAME.vcf

This outputs a .VCF file. You can read it with a text editor to check which SV are detected.

Note: .VCF is the extension for some sort of digital contact card information. Most OS bind the VCF extension to some contact info management program (could be your email program too). Don't double-click, open manually in a text editor. (I disabled the .VCF file association on my computer to simplify my life)

## Evaluate sniffles output

SURVIVOR is used to evalute the output.

> main parameter:
> - max_allowed_distance (?) maximum distance between detected SV and truth-SV (created by our simulation, referenced in BED file)

SURVIVOR eval PATH/TO/GENERATED/VCF PATH/TO/REFERENCE/BED max_allowed_distance output_file

