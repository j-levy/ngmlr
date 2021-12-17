#/usr/bin/env python3

import argparse
import os
import re # regular expressions
import time
import json
import subprocess

utf8 =  'utf-8'

parser = argparse.ArgumentParser(description='Generate SV, generate reads, align them, call SV and compare the result.')


parser.add_argument("--align")
parser.add_argument("--bam_convert")
parser.add_argument("--sv_call")
parser.add_argument("--eval")

parser.add_argument("--read", dest='read', type=str)
parser.add_argument('--ref', dest='ref', type=str)

parser.add_argument('--simulatedepth', dest='depth', type=int, default=20)
parser.add_argument('--threads', dest='threads', type=int, default=os.cpu_count())
parser.add_argument('--snifflescoverage', dest='coverage', type=int, default=8)
parser.add_argument('--subsegment', dest='subsegment', type=int, default=256)
parser.add_argument('--optimize', dest='optimize', type=str, default='pacbio')
parser.add_argument('--nick', dest='nick', type=str, default='')
parser.add_argument('--min_sv_size', dest='min_sv_size', type=int, default=100)

args = parser.parse_args()

dataset = dict()


# map args to the dict. Easier to use
dataset['read'] = args.read
dataset['ref'] = args.ref
dataset['simulatedepth'] = args.depth
dataset['threads'] = args.threads
dataset['snifflescoverage'] = args.coverage
dataset['subsegment'] = args.subsegment
dataset['optimize'] = args.optimize
dataset['nick'] = args.nick
dataset['min_sv_size'] = args.min_sv_size

# add some pre-computed fields to the dict
dataset['nickname'] = f"{args.nick}-{os.path.basename(args.read)}-{os.path.basename(args.ref)}"
dataset['gprof_filename'] = f"gprof.{args.nick}-{args.subsegment}.txt"
dataset['output_filename'] = f"{args.nick}-{args.subsegment}.sam"
dataset['output_dir']= f"out/{args.nick}/threads-{args.threads}/"


def prepare_folder(dataset):
    os.mkdir(dataset['output_dir'])


def create_sv(dataset):
    print("not implemented yet")


def align(git_dir, dataset):
    start_time = time.time()
    os.system(f"""{git_dir}/bin/ngmlr-0.2.8/ngmlr \
        --bam-fix -x {dataset['optimize']} -t {dataset['threads']} \
        --subread-length {dataset['subsegment']} \
        -q {dataset['read']} -r {dataset['ref']} \
        -o {dataset['output_dir']}/{dataset['output_filename']}""")
    end_time= time.time()
    runtime = end_time - start_time
    print(f"+++++ Clock time: {runtime}")
    # echo "+++++ Clock time: $runtime" > $OUTPUTDIR/$GPROFFILENAME
    os.system(f"gprof {git_dir}/bin/ngmlr-0.2.8/ngmlr >> {dataset['output_dir']}/{dataset['gprof_filename']}")

def bam_convert(conda_bin, dataset):
    os.system(f"{conda_bin}/samtools view -o {dataset['output_dir']}/{dataset['output_filename']}.bam {dataset['output_dir']}/{dataset['output_filename']}")
    os.system(f"{conda_bin}/samtools sort -o  {dataset['output_dir']}/{dataset['output_filename']}.sort.bam  {dataset['output_dir']}/{dataset['output_filename']}.bam")
    # )rm {dataset['output_dir']}/{dataset['output_filename']}.bam")
    os.system(f"{conda_bin}/samtools index  {dataset['output_dir']}/{dataset['output_filename']}.sort.bam")

def sv_call(conda_bin, dataset):
    os.system(f"{conda_bin}/sniffles -s {dataset['coverage']} -m {dataset['output_dir']}/{dataset['output_filename']}.sort.bam -v {dataset['output_dir']}/{dataset['output_filename']}.vcf")

def eval(conda_bin, dataset):
    res = subprocess.run(f"{conda_bin}/SURVIVOR", 
        "eval",  
        f"{dataset['output_dir']}/{dataset['output_filename']}.vcf",
        f"{dataset['ref'].replace('.fasta', '.bed')}",
        f"{dataset['min_sv_size']}",
        f"/dev/null", stdout=subprocess.PIPE).stdout.decode(utf8)
    print(res)

def bash_run(command):
    return subprocess.run(command.split(" "), stdout=subprocess.PIPE).stdout.decode("utf-8")


git_dir = bash_run("git rev-parse --show-toplevel").replace("\n", "")
conda_bin = "$HOME/miniconda3/envs/sniffles/bin"

# just remapping output_dir
dataset['output_dir']= f"{git_dir}/docs/{args.nick}/threads-{args.threads}/"

prepare_folder(dataset)
align(git_dir, dataset)
bam_convert(conda_bin, dataset)
sv_call(conda_bin, dataset)
eval(conda_bin, dataset)